% =========================================================================
% replot_agregados.m
%
% Script de pos-processamento que regenera os plots agregados a partir das
% variaveis ja existentes na workspace do MATLAB:
%   agg_phi_err_KW, agg_phi_err_DAS, agg_phi_err_MVDR, agg_phi_err_MUSIC
%   agg_eps_F
%   nCoupling, nInits, nSNR, nISR, nPhi, range_*
%   outDir
%
% Corrige os tres bugs:
%   1. RMSE Frobenius == 0 (cenario Oracle) some em escala log
%   2. Cenarios offline (Coupling 2..5) nao aparecem no plot init=kw_offline
%      por nao terem sido populados em agg_eps_F para iInit=2
%   3. Cores SC-CAPON e SC-KW eram muito parecidas (roxo e violeta)
%
% Uso: rode esse script DEPOIS da simulacao principal.
% =========================================================================

% --- (1) Paleta de cores corrigida (alta distincao visual) ---
cmap_cen = [
    0.20 0.20 0.20;   % 0 No-MC      preto
    0.85 0.20 0.20;   % 1 No-Comp    vermelho
    0.10 0.65 0.10;   % 2 Ideal      verde
    0.60 0.45 0.10;   % 3 Pert       marrom-mostarda
    0.10 0.30 0.85;   % 4 KW-1D      azul
    0.20 0.75 0.80;   % 5 KW-MD      ciano
    0.95 0.50 0.10;   % 6 SC-DAS     laranja
    0.55 0.10 0.65;   % 7 SC-CAPON   roxo
    1.00 0.85 0.15;   % 8 SC-MUSIC   amarelo claro
    0.30 0.00 0.50;   % 9 SC-KW      indigo escuro (bem distinto do roxo)
];
styles_cen = {':','--','-.','-.','-','-','-','-','-','-'};
markers_cen = {'none','none','s','d','d','s','o','o','o','o'};
short_labels_cen = {'No-MC', 'No-Comp', 'Ideal', 'Pert', ...
                    'KW-1D', 'KW-MD', 'SC-DAS', 'SC-CAPON', 'SC-MUSIC', 'SC-KW'};
doa_methods = {'KW', 'DAS', 'MVDR', 'MUSIC'};

% --- (2) Replica agg_eps_F para iInit=2 (cenarios Coupling 2..5) ---
% A simulacao soh popula esses cenarios para iInit=1 (init e' irrelevante p/
% calibracao offline). Para os plots aparecerem em init=kw_offline tambem,
% replicamos o conteudo.
for ic = 1:nCoupling
    coup_code = range_coupling(ic);
    if coup_code >= 2 && coup_code <= 5
        for iI = 2:nInits
            % so copia se o destino estiver vazio (NaN)
            slab_src = agg_eps_F(ic, 1, :, :, :, :);
            agg_eps_F(ic, iI, :, :, :, :) = slab_src;
            % idem para erros de DoA
            agg_phi_err_KW(ic, iI, :, :, :, :)    = agg_phi_err_KW(ic, 1, :, :, :, :);
            agg_phi_err_DAS(ic, iI, :, :, :, :)   = agg_phi_err_DAS(ic, 1, :, :, :, :);
            agg_phi_err_MVDR(ic, iI, :, :, :, :)  = agg_phi_err_MVDR(ic, 1, :, :, :, :);
            agg_phi_err_MUSIC(ic, iI, :, :, :, :) = agg_phi_err_MUSIC(ic, 1, :, :, :, :);
        end
    end
end

% --- (3) Floor numerico para evitar -Inf no log (caso Oracle: eps_F == 0) ---
% Plot em escala log de zero some sem aviso. Usamos um floor pequeno e
% adicionamos uma anotacao para o leitor entender que o valor real e' menor.
EPS_FLOOR = 1e-3;   % piso visual para o RMSE Frobenius

% =========================================================================
% Mapeamento "cenario -> metodo de DoA representativo"
% =========================================================================
method_idx = zeros(nCoupling, 1);
for ic_m = 1:nCoupling
    cc = range_coupling(ic_m);
    switch cc
        case 6, method_idx(ic_m) = 2;   % SC-DAS    -> DAS
        case 7, method_idx(ic_m) = 3;   % SC-CAPON  -> MVDR
        case 8, method_idx(ic_m) = 4;   % SC-MUSIC  -> MUSIC
        case 9, method_idx(ic_m) = 1;   % SC-KW     -> KW
        otherwise, method_idx(ic_m) = 4;
    end
end

% =========================================================================
% Loop de plots por init
% =========================================================================
for iInitPlot_agg = 1:nInits
    init_tag = range_selfcal_init{iInitPlot_agg};

    % --- Pre-calcula matrizes (SNR x ISR) por cenario, agregando sobre angulos ---
    rmse_doa_per_cen = cell(nCoupling, 1);
    rmse_eps_per_cen = cell(nCoupling, 1);
    for ic_m = 1:nCoupling
        switch method_idx(ic_m)
            case 1, A_doa = agg_phi_err_KW;
            case 2, A_doa = agg_phi_err_DAS;
            case 3, A_doa = agg_phi_err_MVDR;
            case 4, A_doa = agg_phi_err_MUSIC;
        end
        slab_doa = squeeze(A_doa(ic_m, iInitPlot_agg, :, :, :, :));
        rmse_doa_per_cen{ic_m} = sqrt(squeeze(mean(slab_doa.^2, [3 4], 'omitnan')));

        slab_eps = squeeze(agg_eps_F(ic_m, iInitPlot_agg, :, :, :, :));
        rmse_eps_per_cen{ic_m} = sqrt(squeeze(mean(slab_eps.^2, [3 4], 'omitnan')));
    end

    % =====================================================================
    % RMSE de DoA: vs SNR (por ISR), vs ISR (por SNR)
    % =====================================================================
    for iISR_fix = 1:nISR
        fig_h = figure('Name', sprintf('RMSE DoA vs SNR (ISR=%+d dB, init=%s)', ...
                       range_ISR_dB(iISR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            curve = rmse_doa_per_cen{ic}(:, iISR_fix);
            method_label = doa_methods{method_idx(ic)};
            plot(range_SNR_dB, curve, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', sprintf('%s (%s)', short_labels_cen{coup_code+1}, method_label));
        end
        set(gca, 'YScale', 'log');
        xlabel('SNR (dB)'); ylabel('RMSE de \phi (graus)');
        title(sprintf('RMSE de DoA vs SNR  |  ISR = %+d dB  |  init self-cal: %s', ...
              range_ISR_dB(iISR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_doa_vs_SNR_ISR_%+d_init_%s.png', range_ISR_dB(iISR_fix), init_tag)), ...
            'Resolution', 200);
    end

    for iSNR_fix = 1:nSNR
        fig_h = figure('Name', sprintf('RMSE DoA vs ISR (SNR=%+d dB, init=%s)', ...
                       range_SNR_dB(iSNR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            curve = rmse_doa_per_cen{ic}(iSNR_fix, :);
            method_label = doa_methods{method_idx(ic)};
            plot(range_ISR_dB, curve, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', sprintf('%s (%s)', short_labels_cen{coup_code+1}, method_label));
        end
        set(gca, 'YScale', 'log');
        xlabel('ISR (dB)'); ylabel('RMSE de \phi (graus)');
        title(sprintf('RMSE de DoA vs ISR  |  SNR = %+d dB  |  init self-cal: %s', ...
              range_SNR_dB(iSNR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_doa_vs_ISR_SNR_%+d_init_%s.png', range_SNR_dB(iSNR_fix), init_tag)), ...
            'Resolution', 200);
    end

    % =====================================================================
    % RMSE Frobenius: vs SNR (por ISR), vs ISR (por SNR)
    % Aplica EPS_FLOOR para nao perder a curva Oracle (eps_F = 0) no log.
    % =====================================================================
    for iISR_fix = 1:nISR
        fig_h = figure('Name', sprintf('RMSE Frob. vs SNR (ISR=%+d dB, init=%s)', ...
                       range_ISR_dB(iISR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            if coup_code < 2, continue; end
            curve = rmse_eps_per_cen{ic}(:, iISR_fix);
            if all(isnan(curve)), continue; end
            % Aplica piso visual: substitui zero/NaN por EPS_FLOOR para
            % que a curva apareca no log.
            curve_plot = max(curve, EPS_FLOOR);
            plot(range_SNR_dB, curve_plot, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', short_labels_cen{coup_code+1});
        end
        set(gca, 'YScale', 'log');
        % Marca o piso visual para deixar claro
        yline(EPS_FLOOR, ':', sprintf('piso visual = %.0e', EPS_FLOOR), ...
              'Color',[0.5 0.5 0.5], 'LabelHorizontalAlignment','left');
        xlabel('SNR (dB)'); ylabel('RMSE de ||C_{hat} - C_{true}||_F / ||C_{true}||_F');
        title(sprintf('RMSE Frobenius vs SNR  |  ISR = %+d dB  |  init: %s', ...
              range_ISR_dB(iISR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_frobenius_vs_SNR_ISR_%+d_init_%s.png', range_ISR_dB(iISR_fix), init_tag)), ...
            'Resolution', 200);
    end

    for iSNR_fix = 1:nSNR
        fig_h = figure('Name', sprintf('RMSE Frob. vs ISR (SNR=%+d dB, init=%s)', ...
                       range_SNR_dB(iSNR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            if coup_code < 2, continue; end
            curve = rmse_eps_per_cen{ic}(iSNR_fix, :);
            if all(isnan(curve)), continue; end
            curve_plot = max(curve, EPS_FLOOR);
            plot(range_ISR_dB, curve_plot, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', short_labels_cen{coup_code+1});
        end
        set(gca, 'YScale', 'log');
        yline(EPS_FLOOR, ':', sprintf('piso visual = %.0e', EPS_FLOOR), ...
              'Color',[0.5 0.5 0.5], 'LabelHorizontalAlignment','left');
        xlabel('ISR (dB)'); ylabel('RMSE de ||C_{hat} - C_{true}||_F / ||C_{true}||_F');
        title(sprintf('RMSE Frobenius vs ISR  |  SNR = %+d dB  |  init: %s', ...
              range_SNR_dB(iSNR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_frobenius_vs_ISR_SNR_%+d_init_%s.png', range_SNR_dB(iSNR_fix), init_tag)), ...
            'Resolution', 200);
    end

    % =====================================================================
    % ACURACIA AGREGADA: scatter (um ponto por cenario)
    % =====================================================================
    rmse_phi_global = zeros(nCoupling, 1);
    rmse_eps_global = zeros(nCoupling, 1);
    for ic = 1:nCoupling
        rmse_phi_global(ic) = sqrt(mean(rmse_doa_per_cen{ic}(:).^2, 'omitnan'));
        rmse_eps_global(ic) = sqrt(mean(rmse_eps_per_cen{ic}(:).^2, 'omitnan'));
    end

    fig_h = figure('Name', sprintf('Acuracia agregada (init=%s)', init_tag), ...
                   'NumberTitle','off', 'Position',[100 100 900 700]);
    hold on; grid on;
    for ic = 1:nCoupling
        coup_code = range_coupling(ic);
        if coup_code < 2 || isnan(rmse_eps_global(ic)), continue; end
        method_label = doa_methods{method_idx(ic)};
        scatter(rmse_phi_global(ic), max(rmse_eps_global(ic), EPS_FLOOR), ...
                180, markers_cen{coup_code+1}, 'filled', ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                'DisplayName', sprintf('%s (%s)', short_labels_cen{coup_code+1}, method_label));
    end
    xline(1, 'k:', '1° de erro de DoA');
    yline(0.1, 'k:', '\epsilon_F = 0.1');
    yline(1, 'r--', '\epsilon_F = 1 (sem comp.)');
    set(gca, 'YScale', 'log');
    xlabel('RMSE de \phi (graus)  [agregado]');
    ylabel('RMSE de ||C_{hat} - C_{true}||_F / ||C_{true}||_F  [agregado]');
    title(sprintf('Acuracia agregada: erro de DoA vs erro de C (init=%s)', init_tag));
    legend('Location','best');

    exportgraphics(fig_h, fullfile(outDir, ...
        sprintf('accuracy_agregada_init_%s.png', init_tag)), ...
        'Resolution', 200);
end

fprintf('\n=== Plots agregados regenerados em %s ===\n', outDir);
