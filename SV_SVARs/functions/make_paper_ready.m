function make_paper_ready(fig, filename)
%MAKE_PAPER_READY  Format and save a MATLAB figure for in a pretty way (PDF).
%
%   Usage:
%       make_paper_ready()                 
%       make_paper_ready('name')           
%       make_paper_ready(figHandle)        
%       make_paper_ready(figHandle, 'name')
%
%   Automatically crops the PDF to the figure size for LaTeX embedding.
%


    % ---- Input handling ----
    if nargin == 0
        fig = gcf; filename = '';
    elseif nargin == 1
        if ischar(fig) || isstring(fig)
            filename = char(fig); fig = gcf;
        else
            filename = '';
        end
    end

    % ---- Determine filename from title if empty ----
    if isempty(filename)
        ax = findall(fig, 'Type', 'axes', '-not', 'Tag', 'legend');
        if ~isempty(ax)
            t = get(get(ax(1), 'Title'), 'String');
            if iscell(t), t = strjoin(t, '_'); end
            if ~isempty(t)
                filename = matlab.lang.makeValidName(erase(string(t), ' '));
            else
                filename = "figure";
            end
        else
            filename = "figure";
        end
    end
    filename = char(filename);

    % ---- Ensure folder exists ----
    if ~exist('figures', 'dir')
        mkdir('figures');
    end

    % ---- Apply nice formatting ----
    set(fig, 'Color', 'w', 'Units', 'inches');
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        % Basic style
        set(ax(k), 'Box', 'off', 'TickDir', 'out', ...
                   'FontName', 'Helvetica', 'FontSize', 10, 'LineWidth', 1);
        % Remove grids
        grid(ax(k), 'off');
        % Adjust title (normal weight)
        t = get(ax(k), 'Title');
        set(t, 'FontWeight', 'normal', 'FontName', 'Helvetica', 'FontSize', 11);
    end

    % ---- Keep legend boxes ----
    % (no legend box modifications — previous code removed them)

    % ---- Crop PDF to figure size ----
    fig_pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPositionMode', 'auto', ...
             'PaperSize', fig_pos(3:4));  % crop to figure size

    % ---- Save as vector PDF ----
    outpath = fullfile('figures_pdf', [filename '.pdf']);
    print(fig, outpath, '-dpdf', '-painters');
    fprintf('Saved clean, cropped figure to %s\n', outpath);
end
