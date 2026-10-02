function plot_irf(irfs)
% PLOT_IRF_SINGLE - Plot IRFs for a single model
%
% Inputs:
%   irfs : [draws x n x n x hor]   Posterior draws of IRFs

    [~, n, ~, horp1] = size(irfs);

    hor = horp1 - 1;

    quants = [.05, .16, .5, .84, .95];
    horizon = (0:hor)';
    xgraph = [horizon', fliplr(horizon')];

    % Quantiles and mean
    IRF_plot = squeeze(quantile(irfs, quants, 1));

    % Plot
    figure;
    for i = 1:n
        for j = 1:n
            idx = (i-1)*n + j;
            subplot(n, n, idx);
            hold on; grid on;

            % Median IRF (black line)
            plot(horizon, squeeze(IRF_plot(3, i,j, :)), 'k', 'LineWidth', 1.3);

            % 90% and 68% credible intervals (blue shaded)
            Y90 = [squeeze(IRF_plot(1,i,j,:))', fliplr(squeeze(IRF_plot(5,i,j,:))')];
            h90 = fill(xgraph, Y90, 'b', 'LineStyle', 'none');
            set(h90, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);

            Y68 = [squeeze(IRF_plot(2,i,j,:))', fliplr(squeeze(IRF_plot(4,i,j,:))')];
            h68 = fill(xgraph, Y68, 'b', 'LineStyle', 'none');
            set(h68, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);

            % Axes and title
            xticks(0:12:hor);
            xlim([horizon(1), horizon(end)]);
            yline(0, '--', 'Color', [0.4 0.4 0.4]);
            title(['IRF(', num2str(i), ',', num2str(j), ')'], 'Interpreter','latex');
            hold off;
        end
    end

end
