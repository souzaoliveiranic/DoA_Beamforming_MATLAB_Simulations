classdef utils
    methods (Static)
        function [X, q, r, bits, pam_rrc_tx, pam_rect, sym_tx, qn, taus_sig, taus_int] = simulate_fsk_data( ...
                M, d, lambda, theta_sig_deg, theta_int_deg, ...
                SNRdB, ISR_dB, N, fs, Rs, sps, alpha, span, fd)

            % Potências dos sinais
            sigma_s2 = 1;
            sigma_i2 = sigma_s2 * 10^(ISR_dB/10);
            sigma_n2 = sigma_s2 / 10^(SNRdB/10);

            % Velocidade da luz
            c = 3e8; % m/s
            fc = c / lambda;           % Hz (usaremos fc para up/down-conversion)

            %% Sinais baseband
            phi0 = 2*pi*rand;
            %q = gen_fm_known(N, fs, fm_msg, deltaf_max, phi0);  % N×1
            [q, bits, pam_rrc_tx, pam_rect, sym_tx] = utils.fsk2_mod(N/sps, Rs, sps, alpha, span, fd);
            % interferidor vai ser também um sinal fm mas com outro modulante
            phi0 = 2*pi*rand;
            %r = gen_fm_known(N, fs, 2*fm_msg*rand, 2*deltaf_max*rand, phi0);  % siN×1
            %r = zeros(N,1); % para testar
            % Interferente: RUÍDO NARROWBAND (mesma banda do FM)
            Bfm = 2*(fd + 0);                 % Carson (Hz)
            %r   = gen_nb_noise(N, fs, Bfm);                % potência ~1, centrado em 0 Hz
            r    = utils.gen_ce_nb_noise(N, fs, Bfm, 0);   % |r|=1 sempre (const-envelope)

            q = q ./ sqrt(mean(abs(q).^2));
            r = r ./ sqrt(mean(abs(r).^2));

            q = sqrt(sigma_s2) * q;
            r = sqrt(sigma_i2) * r;

            %% Atrasos de cada antena em segundos (referência no elemento 0)
            taus_sig = utils.element_delays(M, d, theta_sig_deg, c);   % M×1 [s]
            taus_int = utils.element_delays(M, d, theta_int_deg, c);   % M×1 [s]

            %% steering (broadside)
            a_sig = steering_vec(M, d, lambda, theta_sig_deg);
            a_int = steering_vec(M, d, lambda, theta_int_deg);

            % disp('Atrasos do sinal (s):'); disp(taus_sig.');
            % disp('Atrasos do interferidor (s):'); disp(taus_int.');

            %% Aplica atrasos no envelope + fase da portadora
            Xsig = zeros(M, N);
            Xint = zeros(M, N);
            for m = 1:M
                q_del = delayseq(q, taus_sig(m), fs);  % atraso em tempo (segundos)
                r_del = delayseq(r, taus_int(m), fs);
                phase_sig = exp(-1j*2*pi*fc*taus_sig(m)); % fase devido à portadora
                phase_int = exp(-1j*2*pi*fc*taus_int(m));
                Xsig(m,:) = (q_del(:).') * phase_sig;
                Xint(m,:) = (r_del(:).') * phase_int;
            end


            %% Ruído
            Xn = sqrt(sigma_n2/2) * (randn(M,N) + 1j*randn(M,N));
            %Xn = zeros(M, N);

            %% Sinal total
            X = Xsig + Xint + Xn;
            %qn = q' + Xn(1,:);
            qn = Xsig(1,:) + Xn(1,:);
        end


        function [X, q, v, Xsig, Xint, Xn, bits, pam_rrc_tx, pam_rect, sym_tx, qn, taus_sig, taus_int] = simulate_fsk_data_uca( ...
                M, r, lambda, phi_sig_deg, phi_int_deg, theta_sig_deg, theta_int_deg, ...
                SNRdB, ISR_dB, N, fs, Rs, sps, alpha, span, fd)
            % SIMULATE_FSK_DATA_UCA
            %
            % Simula o sinal recebido por um UCA (Uniform Circular Array) com M elementos,
            % considerando um sinal desejado e um interferidor, ambos modulados em FSK,
            % com ruído AWGN adicionado. Retorna a matriz de sinais do array e os vetores
            % de atraso e steering correspondentes.
            %
            % Entradas:
            %   M, r, lambda      : nº de elementos, raio do UCA, comprimento de onda
            %   phi_sig_deg       : azimute do sinal desejado (graus)
            %   phi_int_deg       : azimute do interferidor (graus)
            %   theta_sig_deg     : elevação do sinal desejado (graus)
            %   theta_int_deg     : elevação do interferidor (graus)
            %   SNRdB, ISR_dB     : SNR e interferer-to-signal ratio (dB)
            %   N, fs, Rs, sps    : nº de amostras, taxa de amostragem, símbolo, oversampling
            %   alpha, span, fd   : roll-off, span do RRC, desvio de frequência FSK
            %
            % Saídas:
            %   X        : M×N matriz de amostras do array (cada linha = antena)
            %   q, v     : sinais baseband desejado/interferente (N×1)
            %   Xsig, Xint, Xn: Sinais que copõem o sinal X
            %   bits     : bits do sinal desejado
            %   pam_rrc_tx, pam_rect, sym_tx : intermediários da modulação
            %   qn       : canal 1 com ruído (referência)
            %   taus_sig, taus_int : vetores M×1 de atrasos (s)
            %

            % --------------------------------------------------------
            % Potências dos sinais
            sigma_s2 = 1;
            sigma_i2 = sigma_s2 * 10^(ISR_dB/10);
            sigma_n2 = sigma_s2 / 10^(SNRdB/10);

            % Constantes
            c = 3e8;                % velocidade da luz (m/s)
            fc = c / lambda;        % frequência portadora (Hz)
            beta = 2*pi*(0:M-1)'/M; % ângulos de cada antena do UCA (rad)

            % --------------------------------------------------------
            % Geração dos sinais baseband
            [q, bits, pam_rrc_tx, pam_rect, sym_tx] = utils.fsk2_mod(N/sps, Rs, sps, alpha, span, fd);
            v = utils.gen_ce_nb_noise(N, fs, 2*fd, 0); % interferente CE narrowband

            % Normalização de potência
            q = q ./ sqrt(mean(abs(q).^2));
            v = v ./ sqrt(mean(abs(v).^2));
            q = sqrt(sigma_s2) * q;
            v = sqrt(sigma_i2) * v;

            % --------------------------------------------------------
            % Atrasos de chegada (UCA)
            taus_sig = utils.element_delays_uca(M, r, theta_sig_deg, phi_sig_deg, c);  % [s]
            taus_int = utils.element_delays_uca(M, r, theta_int_deg, phi_int_deg, c);  % [s]

            % Vetores de direção (UCA)
            a_sig = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_sig_deg);
            a_int = utils.steering_vec_uca(M, r, lambda, theta_int_deg, phi_int_deg);

            % --------------------------------------------------------
            % Aplicação de atrasos e fases
            Xsig = zeros(M, N);
            Xint = zeros(M, N);

            for m = 1:M
                % Atraso temporal no envelope
                q_del = delayseq(q, taus_sig(m), fs);
                v_del = delayseq(v, taus_int(m), fs);

                % Fase da portadora (projeção espacial)
                phase_sig = exp(-1j * 2*pi*fc * taus_sig(m));
                phase_int = exp(-1j * 2*pi*fc * taus_int(m));

                Xsig(m,:) = q_del(:).' * phase_sig;
                Xint(m,:) = v_del(:).' * phase_int;
            end

            % --------------------------------------------------------
            % Adiciona ruído AWGN
            Xn = sqrt(sigma_n2/2) * (randn(M,N) + 1j*randn(M,N));

            % --------------------------------------------------------
            % Sinal total
            X = Xsig + Xint + Xn;
            qn = Xsig(1,:) + Xn(1,:);

            % (opcional) mostrar atrasos
            % fprintf('Atrasos sinal: %s\n', num2str(taus_sig', '%.3e '));
            % fprintf('Atrasos interferente: %s\n', num2str(taus_int', '%.3e '));
        end


        function taus = element_delays(M, d, theta_deg, c)
            % Atrasos de chegada em cada elemento de uma ULA em segundos:
            % tau_m = m * d * sin(theta) / c, m=0..M-1
            m = (0:M-1).';
            theta = deg2rad(theta_deg);
            taus = (m .* d .* sin(theta)) / c;   % [s]
        end

        function taus = element_delays_uca(M, r, theta_deg, phi_deg, c)
            % ELEMENT_DELAYS_UCA  Atrasos de chegada em um UCA (Uniform Circular Array)
            %
            %   taus = element_delays_uca(M, r, theta_deg, phi_deg, c)
            %
            %   M          : número de elementos
            %   r          : raio do UCA (m)
            %   theta_deg  : ângulo de elevação (graus) [0°=z, 90°=plano XY]
            %   phi_deg    : ângulo de azimute (graus)
            %   c          : velocidade da luz (m/s)
            %
            %   Retorna vetor M×1 de atrasos relativos [s] (referência = elemento 1).
            %
            %   Fórmula: τ_m = - (r/c) * sinθ * cos(φ - β_m)

            beta_m = 2*pi*(0:M-1)'/M;  % posição angular de cada antena
            theta  = deg2rad(theta_deg);
            phi    = deg2rad(phi_deg);

            taus = - (r / c) * sin(theta) .* cos(phi - beta_m);
            taus = taus - min(taus);    % referencia atraso mínimo = 0 (opcional)
        end

        function a = steering_vec(M, d, lambda, theta_deg)
            % Vetor de direção ULA (campo distante, banda estreita), ângulo a partir da BROADSIDE.
            % a(m) = exp(-j*2π (m-1) d/λ * sin(theta)), m=1..M
            m = (0:M-1).';
            theta = deg2rad(theta_deg);
            a = exp(-1j*2*pi*(d/lambda)*m*sin(theta));
            a = a / norm(a);  % normalização numérica
        end

        function a = steering_vec_uca(M, r, lambda, theta_deg, phi_deg)
            % STEERING_VEC_UCA  Vetor de direção (steering vector) de um UCA
            %
            %   a = steering_vec_uca(M, r, lambda, theta_deg, phi_deg)
            %
            %   M          : número de elementos
            %   r          : raio do UCA (m)
            %   lambda     : comprimento de onda (m)
            %   theta_deg  : ângulo de elevação (graus)
            %   phi_deg    : ângulo de azimute (graus)
            %
            %   Retorna vetor M×1 complexo: a_m = exp(j*k*r*sinθ*cos(φ−β_m))

            beta_m = 2*pi*(0:M-1)/M;      % posição angular de cada antena
            theta  = deg2rad(theta_deg);
            phi    = deg2rad(phi_deg);
            k      = 2*pi/lambda;

            phi = double(phi(1));
            theta = double(theta(1));

            a = zeros(M,1);
            for m = 1:M
                a(m) = exp(1j * k * r * sin(theta) * cos(phi - beta_m(m)));
            end

            % a = exp(1j * k * r * sin(theta) .* cos(phi - beta_m));
            a = a / norm(a);
        end


        function q = gen_fm_known(N, fs, fm_msg, deltaf_max, phi0)
            % Gera FM baseband
            % m[n] = cos(2π fm_msg n/fs)          (|m[n]| <= 1)
            % θ[n] = θ[n-1] + 2π Δf_max * m[n] / fs
            % s[k]=exp{j*2π*kf​*∑m[j]} somatório de j=0 até k-1

            n  = (0:N-1).';
            m  = cos(2*pi*fm_msg*n/fs);   % sinal modulante
            theta = phi0 + 2*pi*deltaf_max*cumsum(m)/fs;   % integra frequência instantânea
            q = exp(1j*theta);                % FM

            % normaliza potência para ~1
            q = q ./ sqrt(mean(abs(q).^2));
        end


        function r = gen_nb_noise(N, fs, BW)
            % Gera ruído complexo narrowband centrado em 0 Hz, com largura ~BW (Hz)
            % Estratégia: branco complexo -> FIR LP com fc = BW/2 -> normaliza potência.

            fc = BW/2;

            % Ruído branco complexo
            w = (randn(N,1) + 1j*randn(N,1))/sqrt(2);

            % FIR passa-baixas (janela Hamming)
            % ordem ~ fs/BW para transição razoável (clamp para não exagerar)
            ord = max(32, round(8*fs/max(BW,1)));   % heurística simples
            ord = min(ord, 4096);                   % limite superior
            b = fir1(ord, (2*fc)/fs, 'low', hamming(ord+1));  % cutoff normalizado

            % Filtra (fftfilt é eficiente em ordens maiores; filtfilt se quiser fase zero)
            r = fftfilt(b, w);

            % equaliza ganho (descarta transiente inicial simples)
            r = r(ord+1:end);
            if numel(r) < N
                r(end+1:N) = 0;
            else
                r = r(1:N);
            end

            % Normaliza potência para ~1
            r = r ./ sqrt(mean(abs(r).^2) + eps);
        end


        function x = gen_ce_nb_noise(N, fs, BW, f0)
            % Gera "ruído" narrowband de ENVELOPE CONSTANTE (|x|=1), centrado em f0 (Hz).
            % BW ~ largura ocupada (aprox. Carson) do sinal resultante.
            % Estratégia:
            %   - Cria um processo de frequência instantânea f_i[n] = f0 + u[n],
            %     onde u[n] é ruído branco filtrado (LP) com banda ~BW/2
            %   - Modula FM

            % Ruído branco real -> LP para limitar a banda do desvio de frequência
            w   = randn(N,1);
            fcH = min(BW/2, 0.45*fs/2);                % corte "em Hz" do modulante
            ord = max(64, round(8*fs/max(BW,1)));      % ordem FIR (heurística)
            ord = min(ord, 4096);

            % FIR LP (janela Hamming); normalização para [0..1] usa (2*fc)/fs
            b = fir1(ord, (2*fcH)/fs, 'low', hamming(ord+1));
            u = fftfilt(b, w);                          % u ~ ruído LP

            % Ajusta o desvio RMS de frequência (heurístico): std(u) -> BW/4 (Hz)
            u = u / (std(u)+eps) * (BW/4);

            % Frequência instantânea: f_i[n] = f0 + u[n]
            fi = f0 + u;

            % Integra para obter fase
            phi = 2*pi * cumsum(fi) / fs;              % rad
            % opcional: fase inicial aleatória
            phi = phi + 2*pi*rand;

            % Sinal CE
            x = exp(1j*phi);

            % (opcional) normaliza potência — aqui |x|=1, então já está ok
            % x = x / sqrt(mean(abs(x).^2));
        end




        function [s_tx, bits, pam_rrc_tx, pam_rect, sym_tx] = fsk2_mod(Nsym, Rs, sps, alpha, span, fd)
            %FSK2_MOD  Modulador 2-FSK via PAM + RRC + FM complexa BB
            %% Bits -> símbolos PAM (+1/-1)
            fs = Rs * sps;      % taxa de amostragem
            bits       = randi([0 1], Nsym, 1);
            sym_tx     = 2*bits - 1;                % 0->-1, 1->+1

            %% Sinal PAM retangular (zero-order hold) para visualizar
            imp        = upsample(sym_tx, sps);      % impulsos espaçados de sps
            pam_rect   = filter(ones(sps,1), 1, imp);% retangular (PAM "cru")

            %% Filtro RRC (TX) e modelagem do pulso
            rrc = rcosdesign(alpha, span, sps, 'sqrt');  % Root Raised-Cosine
            pam_rrc_tx = filter(rrc, 1, imp);            % PAM modelado por RRC (TX)

            %% Normalização do sinal modulante p/ usar bem o Δf
            m = pam_rrc_tx;
            m = m ./ (max(abs(m)) + eps);   % |m|<=1 => fi = fc ± fd

            %% FM complexa em banda base: s[n] = exp(j*2π ∑(fd*m)/fs)
            phi = 2*pi*cumsum(fd*m)/fs;     % fase acumulada (sem portadora)
            s_tx = exp(1j*phi);             % FM complexa (|s|=1)                    % FM real

        end

        function [bits_hat, BER, pam_rx_mf, sym_rx] = fsk2_demod(s_rx, bits, Rs, sps, alpha, span, fd)
            %FSK2_DEMOD  Demodulador 2-FSK com filtro casado RRC
            %% Demodulação FM por diferença de fase
            fs = Rs * sps;      % taxa de amostragem
            dphi      = angle( s_rx(2:end) .* conj(s_rx(1:end-1)) );  % Δfase
            dphi = dphi(:);   % força coluna
            instfreq  = [dphi(1); dphi] * fs/(2*pi);                  % Hz
            m_hat     = instfreq / fd;                                % estimativa do modulante

            %% Filtro casado RRC (RX) + amostragem & decisão
            rrc = rcosdesign(alpha, span, sps, 'sqrt');  % Root Raised-Cosine
            pam_rx_mf = filter(rrc, 1, m_hat);              % matched filter (sqrtRC)
            gd        = span*sps/2;                         % atraso de grupo por filtro
            % Amostrar no instante ideal (atraso total ~ 2*gd)
            idx0      = 2*gd + 1;

            % Calcula o número máximo de símbolos que cabem no vetor
            Nsym_valid = floor((length(pam_rx_mf) - idx0) / sps);

            sym_rx    = pam_rx_mf(idx0 : sps : idx0 + (Nsym_valid-1)*sps);

            % Decisão dura e BER
            decisions = sign(sym_rx);
            decisions(decisions==0) = 1;
            bits_hat  = decisions > 0;

            bits_ref = bits(1:length(bits_hat));
            BER      = mean(bits_hat ~= bits_ref);

            fprintf('BER: %.2f%% \n', BER*100);

            %% Plot 4×1 dos sinais citados
            nsym_view = 100;                         % quantos símbolos mostrar
            Nview     = nsym_view*sps;
            t1 = (0:Nview-1)/fs * 1e3;              % tempo em ms (apenas para janela)

            % figure('Name','2-FSK via PAM/RRC + FM complexa BB (Δf=4.8 kHz)','Color','w');
            % subplot(4,1,1);
            % plot(t1, pam_rect(1:Nview), 'LineWidth', 1);
            % grid on; xlabel('Tempo (ms)'); ylabel('Amp.');
            % title('1) Sinal PAM retangular (pré-RRC)');
            %
            % subplot(4,1,2);
            % plot(t1, pam_rrc_tx(1:Nview), 'LineWidth', 1);
            % grid on; xlabel('Tempo (ms)'); ylabel('Amp.');
            % title('2) PAM após RRC (TX)');
            %
            % subplot(4,1,3);
            % plot(t1, real(s_tx(1:Nview)), 'b', 'LineWidth', 1); hold on;
            % plot(t1, imag(s_tx(1:Nview)), 'r', 'LineWidth', 1);
            % grid on; xlabel('Tempo (ms)'); ylabel('Amp.');
            % title('3) FM complexa em banda base (Re e Im)');
            %
            % subplot(4,1,4);
            % plot(t1, pam_rx_mf(1:Nview), 'LineWidth', 1);
            % grid on; xlabel('Tempo (ms)'); ylabel('Amp.');
            % title('4) Saída do filtro casado RRC (pré-decisão)');
            %
            % sgtitle(sprintf('2-FSK por FM: Rs=%g, sps=%d, fs=%g kHz, \\Deltaf=%g kHz', ...
            %     Rs, sps, fs/1e3, fd/1e3));

        end


        function [EVM_rms, EVM_dB] = calc_evm_real(rx_sym, tx_sym)
            %CALC_EVM_REAL Calcula EVM para símbolos reais (+1/-1, PAM binário)
            %
            %   [EVM_rms, EVM_dB] = calc_evm_real(rx_sym, tx_sym)
            %
            % Entradas:
            %   rx_sym -> símbolos recebidos (reais)
            %   tx_sym -> símbolos ideais (+1 / -1)
            %
            % Saídas:
            %   EVM_rms -> valor linear (RMS)
            %   EVM_dB  -> valor em dB

            rx_sym = rx_sym(:);
            tx_sym = tx_sym(:);

            % Ajusta tamanhos
            N = min(length(rx_sym), length(tx_sym));
            rx_sym = rx_sym(1:N);
            tx_sym = tx_sym(1:N);

            % Ajuste de ganho por mínimos quadrados
            alpha = (tx_sym' * rx_sym) / (tx_sym' * tx_sym + eps);
            rx_aligned = rx_sym / alpha;
            %rx_aligned = rx_sym;
            % Erro
            err = rx_aligned - tx_sym;

            % EVM
            EVM_rms = sqrt(mean(err.^2) / mean(tx_sym.^2));
            EVM_dB  = 20*log10(EVM_rms);
        end



        function [theta_deg, B_dB] = beampattern_db(w, M, d, lambda, theta_grid_deg, Gdiag)
            % B(θ) = | w^H * (G * a(θ)) |   (G opcional: erros de ganho/fase por canal)
            if nargin < 6 || isempty(Gdiag), Gdiag = ones(M,1); end
            B = zeros(size(theta_grid_deg));
            for k = 1:numel(theta_grid_deg)
                a = utilssteering_vec(M, d, lambda, theta_grid_deg(k));  % Mx1
                B(k) = abs(w' * (Gdiag(:) .* a));                   % |w^H G a|
            end
            B = B ./ max(B + eps);               % normaliza pico em 1
            B_dB = 20*log10(B + eps);            % dB (amplitude)
            theta_deg = theta_grid_deg(:).';
        end

        function [phi_deg, B_dB] = beampattern_db_uca(w, M, r, lambda, phi_grid_deg, Gdiag)
            if nargin < 6 || isempty(Gdiag)
                Gdiag = ones(M,1);
            end

            theta_fixed_deg = 90;   % plano horizontal do UCA
            B = zeros(size(phi_grid_deg));

            for k = 1:numel(phi_grid_deg)
                a = utils.steering_vec_uca(M, r, lambda, theta_fixed_deg, phi_grid_deg(k));  % Mx1
                B(k) = abs( w' * (Gdiag(:) .* a) );
            end

            B = B ./ max(B + eps);         % normaliza pico
            B_dB = 20*log10(B + eps);      % dB
            phi_deg = phi_grid_deg(:).';
        end

        function C = coupling_matrix_uca(M, r, lambda, alpha0)
            % M: nº de elementos
            % r: raio do UCA (m)
            % lambda: comprimento de onda (m)
            % alpha0: fator de acoplamento (0<α<1)
            beta = 2*pi*(0:M-1)'/M;
            k = 2*pi/lambda;

            C = zeros(M);
            for i = 1:M
                for j = 1:M
                    d = 2*r*sin(abs(beta(i)-beta(j))/2);  % distância entre i e j
                    % decaimento exponencial com distância
                    C(i,j) = alpha0^(d/(lambda/2));  % meia-onda = unidade de decaimento
                end
            end
        end

    end
end
