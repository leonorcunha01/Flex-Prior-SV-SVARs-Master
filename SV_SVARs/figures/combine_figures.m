function combined_fig = combine_figures(fig_paths, nrows, ncols)
%COMBINE_FIGURES Combine multiple .fig files into one figure with subplots.

    nfigs = numel(fig_paths);
    combined_fig = figure('Color','w');
    for k = 1:nfigs
        fig = openfig(fig_paths{k}, 'invisible');

        old_axes = findall(fig, 'Type', 'axes', '-not', 'Tag', 'legend');
        if isempty(old_axes)
            warning('No axes found in %s', fig_paths{k});
            close(fig);
            continue;
        end

        ax_new = subplot(nrows, ncols, k, 'Parent', combined_fig);
        copyobj(allchild(old_axes(1)), ax_new);

        title(ax_new, get(get(old_axes(1),'Title'),'String'), ...
              'Interpreter', get(get(old_axes(1),'Title'),'Interpreter'));
        xlabel(ax_new, get(get(old_axes(1),'XLabel'),'String'));
        ylabel(ax_new, get(get(old_axes(1),'YLabel'),'String'));

        % Get relevant objects
        relevant_objects = findall(ax_new, ...
            'Type','line','-or','Type','patch','-or','Type','bar');

        % Flip the order so it matches legend entry order
        relevant_objects = flipud(relevant_objects);

        % Get original legend
        old_legend = findall(fig, 'Type', 'legend');
        legend_strings = old_legend.String;

        % Create legend using corrected ordering
        legend(ax_new, relevant_objects(1:numel(legend_strings)), ...
            legend_strings, 'Location','best');

        close(fig);
    end
end

