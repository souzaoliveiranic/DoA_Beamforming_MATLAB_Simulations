% =========================================================================
% DoA_KW_UCA_Pipeline_BER_EVM
%
% Valida o PIPELINE COMPLETO com modelo de PACOTE (preambulo + payload),
% como nos scripts Conjunto*:
%
%     X (M x K) = [ zeros | PREAMBULO (K_kw) | PAYLOAD (resto) ]
%
%   - PREAMBULO (forma de onda conhecida): estima DoA, b_hat, self-cal (Chat)
%     e a COVARIANCIA do Capon. Daqui saem os PESOS do beamformer.
%   - PAYLOAD: o beamforming (pesos do preambulo) e' aplicado e a BER/EVM
%     sao medidas (dados aleatorios do mesmo SOI).
%
%   Cenario: 1 SOI + 1 interferente co-canal, acoplamento mutuo no canal.
%   DoA p/ apontar: OS 4 METODOS (KW/DAS/Capon/MUSIC) comparados.
%   MUSIC com 1 FONTE. Beamformers: DAS e Capon (MPDR).
%
%   Estrategias (steering do beamformer, estimado no PREAMBULO):
%     - ideal      : sem acoplamento, a(phi_true)        (referencia)
%     - sem comp.  : dado acoplado,  a(phi_hat_metodo)
%     - inversao   : desacopla D=inv(Chat), a(phi_hat_metodo)
%     - manifold   : dado acoplado, steering medido b_hat (= Wiener c/ Capon)
%
%   Figuras (2x2: linhas=beamformer, colunas=estrategia que usa DoA):
%     (1) BER vs SNR (payload) : 4 metodos + manifold + ideal
%     (2) EVM vs SNR (payload) : idem
% =========================================================================

clear; clc; close all;
warning('off','MATLAB:singularMatrix'); warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix'); warning('off','estimate_C_circulant_uca:smallAlpha');

%% ---- Parametros ----
M = 8; fc = 500e6; c = 3e8; lambda = c/fc; r = 0.2*lambda; theta_sig_deg = 90;
fs = 288000; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

% --- estrutura do pacote ---
K_kw        = 2100;             % PREAMBULO (amostras p/ estimacao) -- "o que fazemos hoje"
payload_len = 9900;             % PAYLOAD (amostras p/ medir BER/EVM)
K           = K_kw + payload_len;
k0          = 1;                % preambulo comeca na amostra 1 (sem zeros antes)

ISR_dB   = -2;            % interferente-para-sinal (co-canal). Use -100 p/ desligar.
sep_deg  = 10;          % separacao angular SOI-interferente
interf_t = 'fsk2';      % interferente co-canal (mesma modulacao)

maxIter  = 6; u = 0.5;  % self-cal (KW)
range_SNR_dB = [-12 -9 -6 -3 0 3 6 9 12];
methods = {'KW','DAS','CAPON','MUSIC'};  nM = numel(methods);
bf_list = {'DAS','CAPON'};  nBF = numel(bf_list);
phi_grid_deg = -180:0.5:180;
beta_uca = 2*pi*(0:M-1).'/M;

% Monte Carlo (ajustavel; payload longo -> mantenha moderado)
n_angles = 100; n_trials = 10; n_real = n_angles*n_trials;
rng(2026,'twister');
gstep = phi_grid_deg(2)-phi_grid_deg(1);
phi_set = round((-180 + 360*rand(1,n_angles))/gstep)*gstep;

outDir = fullfile(pwd,'Pipeline BER-EVM Graphs'); if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- C_true ----
C_true = compute_Ctx_for_R(fc, M, r, 50);
% Normaliza C_true para PRESERVAR energia (rede de acoplamento passiva: nao
% cria potencia / nao melhora SNR). Escala s tal que a media angular de
% ||C a(phi)||^2 seja M = ||a||^2 (igual ao caso sem acoplamento). Sem isso,
% como ||C a||^2 pode superar M, o caso "Ideal" apareceria com BER MENOR que
% o "No-MC" (ganho ficticio do acoplamento).
gbar = 0; CtC = C_true'*C_true;
for ig = 1:numel(phi_grid_deg)
    ag = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig));
    gbar = gbar + real(ag'*CtC*ag);
end
gbar = gbar/numel(phi_grid_deg);            % media de ||C a||^2
s_norm = sqrt(M/gbar);  C_true = s_norm*C_true;
fprintf('C_true (raio=%.2f lambda, cond=%.1f, norm s=%.3f, <||Ca||^2>=%.2f->%d), ISR=%+d dB, sep=%d deg, K_kw=%d, payload=%d\n', ...
    r/lambda, cond(C_true), s_norm, gbar, M, ISR_dB, sep_deg, K_kw, payload_len);

%% ---- Acumuladores (BER linear, EVM rms linear) ----
nSNR = numel(range_SNR_dB);
sBER_id  = zeros(nBF,nSNR);    sEVM_id  = zeros(nBF,nSNR);    % No-MC  (sem acoplamento, a_true)
sBER_idl = zeros(nBF,nSNR);    sEVM_idl = zeros(nBF,nSNR);    % Ideal  (acoplado, C_true e angulo conhecidos: s=C_true*a_true)
sBER_kw  = zeros(nBF,nSNR);    sEVM_kw  = zeros(nBF,nSNR);    % BF-KW (beamforming c/ known-waveform: steering medido b_hat)
sBER_sc  = zeros(nBF,nM,nSNR); sEVM_sc  = zeros(nBF,nM,nSNR); % sem comp por metodo (No-Comp = Capon)
sBER_iv  = zeros(nBF,nM,nSNR); sEVM_iv  = zeros(nBF,nM,nSNR); % inversao por metodo
cnt = zeros(1,nSNR);

t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n== SNR=%+d dB ==\n', SNR_dB);
    for ia = 1:n_angles
        phi_sig = phi_set(ia);
        phi_int = wrapTo180(phi_sig + sep_deg);
        for itr = 1:n_trials
            [Xtot, y_ref, k0t, Xsig, Xint, Xn, bits_full, sym_full, bits_pre, ~] = ...
                utils.simulate_data_uca_v3(M, r, lambda, phi_sig, phi_int, theta_sig_deg, theta_sig_deg, ...
                SNR_dB, ISR_dB, K, K_kw, fs, Rs, sps, alpha, span, fd, interf_t, k0);
            q = y_ref(:);  bits_full = bits_full(:);  sym_full = sym_full(:);
            nsym_pre = numel(bits_pre);                  % nº de simbolos do preambulo
            win = k0t : k0t + K_kw - 1;                  % janela do preambulo
            act = k0t : K;                               % parte ativa (preambulo + payload)

            X_ideal_full   = Xtot;                       % = Xsig+Xint+Xn (sem acoplamento)
            X_coupled_full = C_true*(Xsig + Xint) + Xn;  % com acoplamento
            Xk_i = X_ideal_full(:, win);                 % preambulo ideal
            Xk_c = X_coupled_full(:, win);               % preambulo acoplado
            b_hat = Xk_c*conj(q)/(q'*q);                 % steering efetivo do SOI (manifold)
            a_true = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_sig);

            % --- self-cal (KW) no PREAMBULO -> Chat -> desacopla sinal todo ---
            c_fin = selfcal_kw(Xk_c, b_hat, q, M, r, lambda, theta_sig_deg, beta_uca, phi_grid_deg, u, maxIter);
            D = safe_inv(reconstruct_C(c_fin, M));
            Yd_full = D*X_coupled_full;                  % payload desacoplado tambem
            Yd_pre  = Yd_full(:, win);                   % preambulo desacoplado

            % --- DoA p/ apontar (por metodo) no PREAMBULO ---
            a_cm = cell(1,nM); a_dm = cell(1,nM);
            for m = 1:nM
                pc = doa_estimate(Xk_c,   methods{m}, q, r, lambda, beta_uca, phi_grid_deg);
                pd = doa_estimate(Yd_pre, methods{m}, q, r, lambda, beta_uca, phi_grid_deg);
                a_cm{m} = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, pc);
                a_dm{m} = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, pd);
            end

            % --- beamforming (pesos do preambulo) + BER/EVM no PAYLOAD ---
            for ib = 1:nBF
                bf = bf_list{ib};
                % No-MC: canal SEM acoplamento, steering verdadeiro
                [b,e] = bf_demod(X_ideal_full,   Xk_i,  a_true,        bf, act, nsym_pre, bits_full, sym_full, Rs, sps, alpha, span, fd);
                sBER_id(ib,iSNR)=sBER_id(ib,iSNR)+b;    sEVM_id(ib,iSNR)=sEVM_id(ib,iSNR)+e;
                % Ideal: canal acoplado, C_true E angulo conhecidos (s = C_true*a_true)
                [b,e] = bf_demod(X_coupled_full, Xk_c,  C_true*a_true, bf, act, nsym_pre, bits_full, sym_full, Rs, sps, alpha, span, fd);
                sBER_idl(ib,iSNR)=sBER_idl(ib,iSNR)+b;  sEVM_idl(ib,iSNR)=sEVM_idl(ib,iSNR)+e;
                % BF-KW: steering medido b_hat (known-waveform; sem conhecer C nem o angulo)
                [b,e] = bf_demod(X_coupled_full, Xk_c,  b_hat,         bf, act, nsym_pre, bits_full, sym_full, Rs, sps, alpha, span, fd);
                sBER_kw(ib,iSNR)=sBER_kw(ib,iSNR)+b;    sEVM_kw(ib,iSNR)=sEVM_kw(ib,iSNR)+e;
                for m = 1:nM
                    [b,e] = bf_demod(X_coupled_full, Xk_c,  a_cm{m}, bf, act, nsym_pre, bits_full, sym_full, Rs, sps, alpha, span, fd);
                    sBER_sc(ib,m,iSNR)=sBER_sc(ib,m,iSNR)+b;  sEVM_sc(ib,m,iSNR)=sEVM_sc(ib,m,iSNR)+e;
                    [b,e] = bf_demod(Yd_full,        Yd_pre, a_dm{m}, bf, act, nsym_pre, bits_full, sym_full, Rs, sps, alpha, span, fd);
                    sBER_iv(ib,m,iSNR)=sBER_iv(ib,m,iSNR)+b;  sEVM_iv(ib,m,iSNR)=sEVM_iv(ib,m,iSNR)+e;
                end
            end
            cnt(iSNR) = cnt(iSNR) + 1;
        end
        fprintf('  phi_s=%+6.1f  | %.0fs\n', phi_sig, toc(t0));
    end
end

%% ---- Consolidacao ----
cS = reshape(cnt,1,nSNR);  cS3 = reshape(cnt,1,1,nSNR);
BER_id  = sBER_id./cS;   EVM_id  = 20*log10(sEVM_id./cS);    % No-MC
BER_idl = sBER_idl./cS;  EVM_idl = 20*log10(sEVM_idl./cS);   % Ideal
BER_kw  = sBER_kw./cS;   EVM_kw  = 20*log10(sEVM_kw./cS);    % BF-KW
BER_sc  = sBER_sc./cS3;  EVM_sc  = 20*log10(sEVM_sc./cS3);
BER_iv  = sBER_iv./cS3;  EVM_iv  = 20*log10(sEVM_iv./cS3);

% No-Comp (sem compensacao) = sem-comp apontando com Capon
iCap = find(strcmp(methods,'CAPON'));
BER_nc = reshape(BER_sc(:,iCap,:), nBF, nSNR);  EVM_nc = reshape(EVM_sc(:,iCap,:), nBF, nSNR);

nsym_pay = floor(payload_len/sps);
flrB = @(x) max(x, 0.5/(max(cnt)*nsym_pay));     % piso ~ 1 erro no payload

colsM = lines(nM);  mkM = {'o-','s-','d-','^-'};
strat_name = {'sem comp.','inversao (D)'};

%% =================== FIG 1: BER (2x2: beamformer x estrategia) ===================
% Duas versoes:
%   (a) FRACAO    : escala LOG (potencia de 10), sc=1   -> fracao 0..1
%   (b) PORCENTAGEM: escala LINEAR 0..100 (ticks inteiros), sc=100 -> %
% Mesmos dados, so muda escala/limites/ticks do eixo y.
allB = flrB([BER_sc(:); BER_iv(:); BER_id(:); BER_idl(:); BER_kw(:)]);  ylB0 = [min(allB)*0.7, 1];
ber_modes = { 'frac', 1,   'BER (payload)',    'pipeline_BER_vs_snr.png',     'log'; ...
              'pct',  100, 'BER (%, payload)', 'pipeline_BER_pct_vs_snr.png', 'linear' };
for imode = 1:size(ber_modes,1)
    sc    = ber_modes{imode,2};
    ylab  = ber_modes{imode,3};
    fout  = ber_modes{imode,4};
    yscl  = ber_modes{imode,5};
    fig = figure('Color','w','Position',[40 50 1250 860]);
    for ib = 1:nBF
        for js = 1:2
            subplot(2,2,(ib-1)*2+js); hold on; grid on;
            if js==1, Bm = BER_sc; else, Bm = BER_iv; end
            for m = 1:nM
                plot(range_SNR_dB, sc*flrB(squeeze(Bm(ib,m,:))), mkM{m}, 'Color',colsM(m,:), ...
                    'LineWidth',1.7,'MarkerFaceColor',colsM(m,:),'MarkerSize',6,'DisplayName',methods{m});
            end
            % referencias: No-MC (sem acoplamento), Ideal (C_true e phi conhecidos), No-Comp (Capon sem comp)
            plot(range_SNR_dB, sc*flrB(BER_id(ib,:)),  'k--','LineWidth',1.6,'DisplayName','No-MC');
            plot(range_SNR_dB, sc*flrB(BER_idl(ib,:)), 'k-', 'LineWidth',2.4,'DisplayName','Ideal');
            plot(range_SNR_dB, sc*flrB(BER_nc(ib,:)),  ':', 'Color',[0.55 0.55 0.55],'LineWidth',2.2,'DisplayName','No-Comp (Capon)');
            plot(range_SNR_dB, sc*flrB(BER_kw(ib,:)),  '-.','Color',[0.60 0 0.70],'LineWidth',2.2,'DisplayName','BF-KW');
            set(gca,'YScale',yscl); xlabel('SNR (dB)'); ylabel(ylab); xticks(range_SNR_dB);
            if strcmp(yscl,'linear')
                ylim([0 55]); yticks(0:10:55);     % % com numeros inteiros 0..100
            else
                ylim([min(allB)*0.7, 1]);            % fracao em potencias de 10
            end
            title(sprintf('%s  |  %s', bf_list{ib}, strat_name{js}));
            if ib==1 && js==1, legend('Location','southwest'); end
        end
    end
    unit_str = ber_modes{imode,1}; if strcmp(unit_str,'pct'), unit_str='%'; else, unit_str='fracao'; end
    sgtitle(sprintf('BER (payload, %s) vs SNR por metodo de DoA  (raio=%.2f\\lambda, ISR=%+d dB, sep=%d^o)', unit_str, r/lambda, ISR_dB, sep_deg),'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir,fout),'Resolution',170);
    matlab2tikz(fullfile(outDir, regexprep(fout,'\.png$','.tex')), 'width','\figurewidth','height','\figureheight');
end

%% =================== FIG 2: EVM (2x2: beamformer x estrategia) ===================
allE = [EVM_sc(:); EVM_iv(:); EVM_id(:); EVM_idl(:); EVM_kw(:)];  ylE = [min(allE)-1, max(allE)+1];
fig = figure('Color','w','Position',[40 50 1250 860]);
for ib = 1:nBF
    for js = 1:2
        subplot(2,2,(ib-1)*2+js); hold on; grid on;
        if js==1, Em = EVM_sc; else, Em = EVM_iv; end
        for m = 1:nM
            plot(range_SNR_dB, squeeze(Em(ib,m,:)), mkM{m}, 'Color',colsM(m,:), ...
                'LineWidth',1.7,'MarkerFaceColor',colsM(m,:),'MarkerSize',6,'DisplayName',methods{m});
        end
        % referencias: No-MC, Ideal (C_true e phi conhecidos), No-Comp (Capon sem comp)
        plot(range_SNR_dB, EVM_id(ib,:),  'k--','LineWidth',1.6,'DisplayName','No-MC');
        plot(range_SNR_dB, EVM_idl(ib,:), 'k-', 'LineWidth',2.4,'DisplayName','Ideal');
        plot(range_SNR_dB, EVM_nc(ib,:),  ':', 'Color',[0.55 0.55 0.55],'LineWidth',2.2,'DisplayName','No-Comp (Capon)');
        plot(range_SNR_dB, EVM_kw(ib,:),  '-.','Color',[0.60 0 0.70],'LineWidth',2.2,'DisplayName','BF-KW');
        xlabel('SNR (dB)'); ylabel('EVM (dB, payload)'); xticks(range_SNR_dB); ylim(ylE);
        title(sprintf('%s  |  %s', bf_list{ib}, strat_name{js}));
        if ib==1 && js==1, legend('Location','northeast'); end
    end
end
sgtitle(sprintf('EVM (payload) vs SNR por metodo de DoA  (raio=%.2f\\lambda, ISR=%+d dB, sep=%d^o)', r/lambda, ISR_dB, sep_deg),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'pipeline_EVM_vs_snr.png'),'Resolution',170);
matlab2tikz(fullfile(outDir,'pipeline_EVM_vs_snr.tex'), 'width','\figurewidth','height','\figureheight');

%% ---- Resumo numerico (Capon) ----
ibC = find(strcmp(bf_list,'CAPON'));
fprintf('\n===== BER (%%) payload, Capon, INVERSAO por metodo, vs SNR =====\n');
fprintf('%-10s', 'SNR'); fprintf('%9d', range_SNR_dB); fprintf('\n');
for m=1:nM, fprintf('%-10s', methods{m}); fprintf('%9.3f', 100*squeeze(BER_iv(ibC,m,:))); fprintf('\n'); end
fprintf('%-10s', 'BF-KW');    fprintf('%9.3f', 100*BER_kw(ibC,:));  fprintf('\n');
fprintf('--- referencias ---\n');
fprintf('%-10s', 'No-MC');    fprintf('%9.3f', 100*BER_id(ibC,:));  fprintf('\n');
fprintf('%-10s', 'Ideal');    fprintf('%9.3f', 100*BER_idl(ibC,:)); fprintf('\n');
fprintf('%-10s', 'No-Comp');  fprintf('%9.3f', 100*BER_nc(ibC,:));  fprintf('\n');

%% =================== FIG 3: BEAMPATTERNS (cenario representativo) ===================
% UM unico pacote (angulos e SNR fixos). Mostra o PADRAO ESPACIAL EFETIVO de
% cada estrategia: resposta |w_eff^H h(phi)|^2 a uma fonte real chegando de phi
% NO CANAL REAL com acoplamento, h(phi)=C_true*a(phi). Normalizado ao maximo.
% Marca o SOI e o interferente -> da' p/ ver onde aponta o lobulo e se ha' nulo.
rng(7,'twister');                                   % cenario reprodutivel
phi_bp_sig = 40;  SNR_bp = 6;                       % AJUSTAVEL
phi_bp_int = 60; %wrapTo180(phi_bp_sig + sep_deg);
[Xtot_b, yb, k0b, Xs_b, Xi_b, Xn_b] = utils.simulate_data_uca_v3(M, r, lambda, ...
    phi_bp_sig, phi_bp_int, theta_sig_deg, theta_sig_deg, SNR_bp, ISR_dB, K, K_kw, ...
    fs, Rs, sps, alpha, span, fd, interf_t, k0);
qb = yb(:);  winb = k0b:k0b+K_kw-1;
Xc_full_b = C_true*(Xs_b + Xi_b) + Xn_b;            % preambulo acoplado
Xk_ib = Xtot_b(:,winb);  Xk_cb = Xc_full_b(:,winb);
b_hatb  = Xk_cb*conj(qb)/(qb'*qb);
a_trueb = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_bp_sig);
c_finb  = selfcal_kw(Xk_cb, b_hatb, qb, M, r, lambda, theta_sig_deg, beta_uca, phi_grid_deg, u, maxIter);
Db      = safe_inv(reconstruct_C(c_finb, M));
Ydpre_b = Db*Xk_cb;                                 % preambulo desacoplado
% steerings: No-Comp aponta com Capon (acoplado); inversao aponta com KW (desacoplado)
a_nc   = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, doa_estimate(Xk_cb,  'CAPON', qb, r, lambda, beta_uca, phi_grid_deg));
a_ivkw = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, doa_estimate(Ydpre_b,'KW',    qb, r, lambda, beta_uca, phi_grid_deg));

% dicionario de steering e manifold real (com acoplamento) por azimute
A_plot = zeros(M, numel(phi_grid_deg));
for ig = 1:numel(phi_grid_deg), A_plot(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig)); end
H_plot = C_true*A_plot;

fig = figure('Color','w','Position',[40 60 1300 560]);
for ib = 1:nBF
    bf = bf_list{ib};
    w_nomc = bf_weights(Xk_ib,   a_trueb,        bf);   % opera no canal SEM acoplamento
    w_idl  = bf_weights(Xk_cb,   C_true*a_trueb, bf);
    w_kw   = bf_weights(Xk_cb,   b_hatb,         bf);
    w_nc   = bf_weights(Xk_cb,   a_nc,           bf);
    w_iv   = bf_weights(Ydpre_b, a_ivkw,         bf);   % opera no canal DESACOPLADO
    subplot(1,nBF,ib); hold on; grid on;
    plot(phi_grid_deg, bp_db(w_nomc, A_plot),    'k--','LineWidth',1.4,'DisplayName','No-MC');
    plot(phi_grid_deg, bp_db(w_idl,  H_plot),    'k-', 'LineWidth',1.8,'DisplayName','Ideal');
    plot(phi_grid_deg, bp_db(w_kw,   H_plot),    '-.', 'Color',[0.60 0 0.70],'LineWidth',2.0,'DisplayName','BF-KW');
    plot(phi_grid_deg, bp_db(w_nc,   H_plot),    ':',  'Color',[0.55 0.55 0.55],'LineWidth',2.0,'DisplayName','No-Comp (Capon)');
    plot(phi_grid_deg, bp_db(Db'*w_iv, H_plot),  '-',  'Color',[0.10 0.25 0.85],'LineWidth',1.3,'DisplayName','inversao (KW)');
    xline(phi_bp_sig,'-','Color',[0 0.6 0],'LineWidth',1.6,'DisplayName','SOI');
    xline(phi_bp_int,'-','Color',[0.85 0 0],'LineWidth',1.6,'DisplayName','interferente');
    xlim([-180 180]); xticks(-180:45:180); ylim([-50 2]);
    xlabel('azimute \phi (graus)'); ylabel('ganho normalizado (dB)');
    title(bf); if ib==1, legend('Location','southwest','FontSize',8); end
end
sgtitle(sprintf('Beampatterns efetivos no canal real  (SOI=%d^o, interf=%d^o, SNR=%+d dB, ISR=%+d dB, raio=%.2f\\lambda)', ...
    phi_bp_sig, phi_bp_int, SNR_bp, ISR_dB, r/lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'pipeline_beampatterns.png'),'Resolution',170);
matlab2tikz(fullfile(outDir,'pipeline_beampatterns.tex'), 'width','\figurewidth','height','\figureheight');

fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function [ber, evm] = bf_demod(Y_full, Y_pre, s, bf, act, nsym_pre, bits_full, sym_full, Rs, sps, alpha, span, fd)
% Pesos formados no PREAMBULO (Y_pre): DAS=s/(s'*s) ou Capon/MPDR com a
% covariancia do preambulo. Aplica ao sinal todo (Y_full) e mede BER/EVM
% APENAS no PAYLOAD (simbolos apos os nsym_pre do preambulo).
    w = bf_weights(Y_pre, s, bf);
    y = w'*Y_full;                                            % aplica ao sinal todo (1 x K)
    [bits_hat, ~, ~, sym_rx] = utils.fsk2_demod(y(act), bits_full, Rs, sps, alpha, span, fd);
    i0 = nsym_pre + 1;                                        % primeiro simbolo do payload
    Lb = min(numel(bits_hat), numel(bits_full));
    if i0 <= Lb, ber = mean(bits_hat(i0:Lb) ~= bits_full(i0:Lb)); else, ber = 0.5; end
    Ls = min(numel(sym_rx), numel(sym_full));
    if i0 <= Ls, [evm,~] = utils.calc_evm_real(sym_rx(i0:Ls), sym_full(i0:Ls)); else, evm = 1.0; end
    if ~isfinite(ber), ber = 0.5; end
    if ~isfinite(evm), evm = 1.0; end
end

function w = bf_weights(Y_pre, s, bf)
% Pesos do beamformer formados no PREAMBULO Y_pre.
%   DAS  : casado          w = s/(s's)
%   CAPON: MPDR/Capon      w = R^{-1}s / (s' R^{-1} s),  R = cov. do preambulo
    M = size(Y_pre,1);
    if strcmpi(bf,'DAS')
        w = s/(s'*s + eps);
    else  % CAPON / MPDR (covariancia do preambulo, com diagonal loading)
        Kp = size(Y_pre,2);  R = (Y_pre*Y_pre')/Kp;  R = R + 1e-3*trace(R)/M*eye(M);
        ws = R\s;  w = ws/(s'*ws + eps);
    end
end

function PdB = bp_db(w, H)
% Resposta espacial |w^H h(phi)|^2 (dB), normalizada ao maximo (0 dB).
% H(:,k) = manifold da fonte que chega de phi_k (a(phi) ou C*a(phi)).
    p = abs(w'*H).^2;
    PdB = 10*log10(p/(max(p)+eps));
end

function c_fin = selfcal_kw(X_coupled, b_hat_orig, q, M, r, lambda, theta_deg, beta_uca, phi_grid_deg, u, maxIter)
    phi0 = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, phi_grid_deg);
    C_hat = estimate_C_circulant_uca(b_hat_orig, utils.steering_vec_uca(M,r,lambda,theta_deg,phi0), M);
    for it = 1:maxIter
        D = safe_inv(C_hat);  Y = D*X_coupled;
        phi = doa_estimate(Y, 'KW', q, r, lambda, beta_uca, phi_grid_deg);
        a_hat = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
        C_ls = estimate_C_circulant_uca(b_hat_orig, a_hat, M);
        C_hat = (1-u)*C_hat + u*C_ls;
    end
    c_fin = C_hat(1, 2:floor(M/2)+1).';
end

function C = reconstruct_C(c, M)
    K = floor(M/2); is_even = (mod(M,2)==0); c = c(:);
    if is_even, first_row = [1, c(1:K-1).', c(K), flip(c(1:K-1).')];
    else,       first_row = [1, c(1:K).', flip(c(1:K).')]; end
    C = zeros(M);
    for i = 1:M, C(i,:) = circshift(first_row, [0, i-1]); end
end

function D = safe_inv(C)
    M = size(C,1);
    ws = warning('off','MATLAB:singularMatrix'); wn = warning('off','MATLAB:nearlySingularMatrix');
    cleanupObj = onCleanup(@() warning([ws wn]));
    rc = rcond(C);
    if isfinite(rc) && rc > 1e-12, D = inv(C); if all(isfinite(D(:))), return; end, end
    mu = 1e-6*(norm(C,'fro')^2/M + eps);  D = (C'*C + mu*eye(M))\C';
    if ~all(isfinite(D(:))), D = eye(M); end
end

function phi_hat = doa_estimate(X, method, q, r, lambda, beta_uca, phi_grid_deg)
% DoA do SOI. KW usa a forma de onda conhecida (robusto ao interferente);
% DAS/Capon/MUSIC varrem o espectro assumindo 1 FONTE dominante.
    M = size(X,1);
    persistent A_dict_cache phi_cache r_cache
    if isempty(phi_cache) || ~isequal(phi_cache,phi_grid_deg) || ~isequal(r_cache,r)
        A_dict_cache = zeros(M,numel(phi_grid_deg));
        for ig=1:numel(phi_grid_deg), A_dict_cache(:,ig)=utils.steering_vec_uca(M,r,lambda,90,phi_grid_deg(ig)); end
        phi_cache = phi_grid_deg;  r_cache = r;
    end
    A_dict = A_dict_cache;
    switch upper(method)
        case 'KW'
            [~, phi] = doa_kw_uca(X, q(:).', r, lambda, beta_uca);  phi_hat = phi(1);
        case 'DAS'
            R = (X*X')/size(X,2);
            scores = real(sum(conj(A_dict).*(R*A_dict), 1));
            [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        case 'CAPON'
            R = (X*X')/size(X,2); R = R + 1e-6*trace(R)/M*eye(M); Rinv = R\eye(M);
            den = real(sum(conj(A_dict).*(Rinv*A_dict),1));
            [~,ip] = max(1./max(den,eps)); phi_hat = phi_grid_deg(ip);
        case 'MUSIC'
            R = (X*X')/size(X,2); R = (R+R')/2;
            [V, Dg] = eig(R); [~, idx] = sort(real(diag(Dg)), 'descend'); V = V(:, idx);
            En = V(:, 2:end);            % 1 FONTE (estima so' o SOI dominante)
            proj = En'*A_dict;
            den = real(sum(conj(proj).*proj, 1));
            [~,ip] = max(1./max(den,eps)); phi_hat = phi_grid_deg(ip);
        otherwise
            error('doa_estimate: metodo %s', method);
    end
end

function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
    c = 3e8; lambda = c/fc;
    mp = dipole; mp.Length = 0.5*lambda; mp.Width = 0.01*lambda;
    uca = circularArray; uca.Element = mp; uca.NumElements = M; uca.Radius = R;
    Sobj = sparameters(uca, fc);
    Z_matrix = s2z(Sobj.Parameters(:,:,1), Z0);
    denom = diag(Z_matrix) + Z0;
    Ctx = eye(M);
    for j = 1:M
        for i = 1:M
            if i ~= j, Ctx(i,j) = Z_matrix(i,j) / denom(j); end
        end
    end
end
