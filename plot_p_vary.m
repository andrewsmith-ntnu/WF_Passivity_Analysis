function [] = plot_p_vary(f, Z1_array, Z2_array, n_o, n_i, opts)
% Plot passivity and impedance variation across datasets
%
% plot_p_vary(f, Z1_array, Z2_array, n_o, n_i, opts)
%
% Inputs
%   f         - figure handle or figure index
%   Z1_array  - cell array of FRD objects used to evaluate passivity
%   Z2_array  - cell array of FRD objects to plot impedance magnitude/imag
%   n_o, n_i  - output/input indices (for selecting matrix element)
%   opts      - optional struct: fields 'plot_pass' (bool), 'plot_imp' (bool), 'plot_type' ('imag')
%
% This utility plots passivity (min eigenvalue of Hermitian part) and
% either |Z| or Im(Z) for several parameter sweeps contained in Z2_array.
% Single impedances can be passed to Z1 or Z2 without a sweep.

figure(f);
f_size = 16;
persistent colorIdx colorFig


if iscell(Z2_array)
    n_p = length(Z2_array);
    freq_array=Z2_array{1}.Frequency;
    color_plot = colormap(jet(max(1,n_p)));
else
    n_p=1;
    freq_array=Z2_array.Frequency;
    Z2_array={Z2_array};
    Z1_array={Z1_array};
    
    % Reset colorIdx if figure is new
    if isempty(colorFig) || ~isvalid(colorFig) || colorFig ~= gcf
        colorIdx = 1;
        colorFig = gcf;
    end
    ax = gca;
    colorOrder = ax.ColorOrder;
    color_plot = colorOrder(mod(colorIdx-1, size(colorOrder,1)) + 1, :);
    colorIdx = colorIdx + 1;
end



% Decide layout: two rows when plotting passivity, otherwise single plot
if exist('opts', 'var') && isfield(opts, 'plot_pass') && opts.plot_pass == 1 && ((isfield(opts, 'plot_imp') && opts.plot_imp == 1) || ~isfield(opts, 'plot_imp'))
    t = findall(gcf, 'Type', 'tiledlayout');
    if isempty(t)
        t = tiledlayout(2,1);
    end
    t.TileSpacing='compact';
    if length(Z1_array) == 1

        Z1 = Z1_array{1};
        pass = mimo_passivity(Z1);

        % Find zero-crossings in passivity
        signal = pass(:);
        signChanges = diff(sign(signal));
        zeroCrossings = find(signChanges ~= 0);

        % Plot passivity of Z1
        nexttile(t,1);
        p = semilogx(Z1.Frequency, pass); % Use default color
        hold on
        ylabel('Passivity');
        grid on;
        p.LineWidth = 1.5;
        yline(0);
        set(gca,'XTickLabel',[]);

        for i = 1:length(zeroCrossings)
            xline(Z1.Frequency(zeroCrossings(i)), 'k--');
        end

        nexttile(t,2);
        for i = 1:n_p
            color = color_plot(i,:);
            Z2 = Z2_array{i};
            mag = abs(Z2.ResponseData);
            im_part = imag(Z2.ResponseData);

            % Plot magnitude or imaginary part depending on opts
            if exist('opts', 'var') && isfield(opts, 'plot_type') && strcmp(opts.plot_type, 'imag')
                p = semilogx(freq_array, squeeze(im_part(n_o,n_i,:)), 'Color', color);
                p.LineWidth = 1.5;
                if isempty(Z2.Name)
                    ylabel(['$Im[Z_{', Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)]$'], 'Interpreter', 'latex');
                else
                    ylabel(['$Im[',Z2.Name(1:end-1), ',',Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)]$'], 'Interpreter', 'latex');
                end
                hold on
                yline(0,'k');
            else
                p = loglog(freq_array, squeeze(mag(n_o,n_i,:)), 'Color', color);
                p.LineWidth = 1.5;
                if isempty(Z2.Name)
                    ylabel(['$|Z_{', Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)|$'], 'Interpreter', 'latex');
                else
                    ylabel(['$|',Z2.Name(1:end-1), ',',Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)|$'], 'Interpreter', 'latex');
                end
                hold on;
            end
        end
        for i = 1:length(zeroCrossings)
            xline(Z1_array{1}.Frequency(zeroCrossings(i)), 'k--');
        end
    else
        % Multiple Z1 sets: plot each with matching color
        for i = 1:n_p
            color = color_plot(i,:);
            Z1 = Z1_array{i};
            Z2 = Z2_array{i};
            mag = abs(Z2.ResponseData);
            im_part = imag(Z2.ResponseData);
            pass = mimo_passivity(Z1);

            % Plot passivity of Z1
            nexttile(t,1);
            p = semilogx(Z1.Frequency, pass, 'Color', color);
            hold on
            ylabel('Passivity');
            grid on;
            p.LineWidth = 1.5;
            yline(0);
            set(gca,'XTickLabel',[]);

            % Plot magnitude/imag part of Z2
            nexttile(t,2);
            if exist('opts', 'var') && isfield(opts, 'plot_type') && strcmp(opts.plot_type, 'imag')
                p = semilogx(freq_array, squeeze(im_part(n_o,n_i,:)), 'Color', color);
                grid on
                p.LineWidth = 1.5;
                if isempty(Z2.Name)
                    ylabel(['$Im[Z_{', Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)]$'], 'Interpreter', 'latex');
                else
                    ylabel(['$Im[',Z2.Name(1:end-1), ',',Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)]$'], 'Interpreter', 'latex');
                end
                hold on;
                yline(0,'k');
            else
                p = loglog(freq_array, squeeze(mag(n_o,n_i,:)), 'Color', color);
                p.LineWidth = 1.5;
                grid on
                if isempty(Z2.Name)
                    ylabel(['$|Z_{', Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)|$'], 'Interpreter', 'latex');
                else
                    ylabel(['$|',Z2.Name(1:end-1), ',',Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)|$'], 'Interpreter', 'latex');
                end
                hold on;
            end
        end
    end
elseif exist('opts', 'var') && isfield(opts, 'plot_pass') && opts.plot_pass == 1 && isfield(opts, 'plot_imp') && opts.plot_imp == 0 
    t = findall(gcf, 'Type', 'tiledlayout');
    if isempty(t)
        t = tiledlayout(1,1);
    end
    t.TileSpacing='compact';
    nexttile(t,1);
    if length(Z1_array) == 1
        
        Z1 = Z1_array{1};
        pass = mimo_passivity(Z1);

        % Find zero-crossings in passivity
        signal = pass(:);
        signChanges = diff(sign(signal));
        zeroCrossings = find(signChanges ~= 0);

        % Plot passivity of Z1
        p = semilogx(Z1.Frequency, pass); % Use default color
        hold on
        ylabel('Passivity');
        grid on;
        p.LineWidth = 1.5;
        yline(0);

        for i = 1:length(zeroCrossings)
            xline(Z1.Frequency(zeroCrossings(i)), 'k--');
        end

    else
        % Multiple Z1 sets: plot each with matching color
        for i = 1:n_p
            color = color_plot(i,:);
            Z1 = Z1_array{i};
            pass = mimo_passivity(Z1);

            % Plot passivity of Z1
            p = semilogx(Z1.Frequency, pass, 'Color', color);
            hold on
            ylabel('Passivity');
            grid on;
            p.LineWidth = 1.5;
            yline(0);
        end
    end

else
    % Only impedance plots
    t = findall(gcf, 'Type', 'tiledlayout');
    if isempty(t)
        t = tiledlayout(1,1);
    end
    t.TileSpacing='compact';
    nexttile(t,1);

    for i = 1:n_p
        color = color_plot(i,:);
        Z2 = Z2_array{i};
        mag = abs(Z2.ResponseData);
        im_part = imag(Z2.ResponseData);

        if exist('opts', 'var') && isfield(opts, 'plot_type') && strcmp(opts.plot_type, 'imag')
            p = semilogx(freq_array, squeeze(im_part(n_o,n_i,:)), 'Color', color);
            p.LineWidth = 1.5;
            if isempty(Z2.Name)
                ylabel(['$Im[Z_{', Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)]$'], 'Interpreter', 'latex');
            else
                ylabel(['$Im[',Z2.Name(1:end-1), ',',Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)]$'], 'Interpreter', 'latex');
            end
            hold on
            yline(0,'k');
        else
            p = loglog(freq_array, squeeze(mag(n_o,n_i,:)), 'Color', color);
            p.LineWidth = 1.5;
            if isempty(Z2.Name)
                ylabel(['$|Z_{', Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)|$'], 'Interpreter', 'latex');
            else
                ylabel(['$|',Z2.Name(1:end-1), ',',Z2.InputName{n_i}(end), Z2.OutputName{n_o}(end), '}(j\omega)|$'], 'Interpreter', 'latex');
            end
            hold on;
        end
    end
    if length(Z1_array) == 1
        Z1 = Z1_array{1};
        pass = mimo_passivity(Z1);

        %Find zero-crossings in passivity
        signal = pass(:);
        signChanges = diff(sign(signal));
        zeroCrossings = find(signChanges ~= 0);

        for i = 1:length(zeroCrossings)
            xline(Z1_array{1}.Frequency(zeroCrossings(i)), 'k--');
        end
    end

end

grid on
xlabel('Frequency (Hz)', 'interpreter', 'latex');
ax = gca;
ax.FontSize = f_size;


nexttile(t,1);
ax = gca;
ax.FontSize = f_size;

end

