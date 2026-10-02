function plot_irf_compare(irfs1, irfs2, varnames, shocknames)
% PLOT_IRF_COMPARE - Plot IRFs for two models on the same axes
%
% Inputs:
%   irfs1      : IRF draws of first model
%   irfs2      : IRF draws of second model
%   varnames   : cell array of variable names (length n)
%   shocknames : cell array of shock names (length n)

    [~, n, ~, horp1] = size(irfs1);
    hor = horp1 - 1;

    quants = [.05, .16, .5, .84, .95];
    horizon = (0:hor)';
    xgraph = [horizon', fliplr(horizon')];

    % ---- Quantiles ----
    IRF_plot1 = squeeze(quantile(irfs1, quants, 1));
    IRF_plot2 = squeeze(quantile(irfs2, quants, 1));

    figure;
    for i = 1:n
        for j = 1:n
            idx = (i-1)*n + j;
            subplot(n, n, idx);
            hold on; grid on;

            % --- First model ---
            plot(horizon, squeeze(IRF_plot1(3,i,j,:)), 'k', 'LineWidth', 1.3);
            Y90_1 = [squeeze(IRF_plot1(1,i,j,:))', fliplr(squeeze(IRF_plot1(5,i,j,:))')];
            fill(xgraph, Y90_1, 'b', 'LineStyle', 'none', 'FaceAlpha',0.2);
            Y68_1 = [squeeze(IRF_plot1(2,i,j,:))', fliplr(squeeze(IRF_plot1(4,i,j,:))')];
            fill(xgraph, Y68_1, 'b', 'LineStyle', 'none', 'FaceAlpha',0.5);

            % --- Second model ---
            plot(horizon, squeeze(IRF_plot2(3,i,j,:)), 'r', 'LineWidth', 1.3);
            Y90_2 = [squeeze(IRF_plot2(1,i,j,:))', fliplr(squeeze(IRF_plot2(5,i,j,:))')];
            fill(xgraph, Y90_2, [1 .5 0], 'LineStyle','none','FaceAlpha',0.2);
            Y68_2 = [squeeze(IRF_plot2(2,i,j,:))', fliplr(squeeze(IRF_plot2(4,i,j,:))')];
            fill(xgraph, Y68_2, [1 .5 0], 'LineStyle','none','FaceAlpha',0.5);

            xticks(0:12:hor);
            xlim([horizon(1), horizon(end)]);
            yline(0, '--', 'Color', [0.4 0.4 0.4]);

            % ---- Add row titles (variables) ----
            if j == 1
                ylabel(varnames{i}, 'Interpreter','latex');
            end

            % ---- Add column titles (shocks) ----
            if i == 1
                title(shocknames{j}, 'Interpreter','latex');
            end

            hold off;
        end
    end
end
