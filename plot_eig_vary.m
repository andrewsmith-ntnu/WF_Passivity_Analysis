function [] = plot_eig_vary(f,eig_array)
% Plot passivity and impedance variation across datasets
%
% plot_p_vary(f, Z1_array, Z2_array, n_o, n_i, opts)
%
% Inputs
%   f         - figure handle or figure index
%   Z1_array  - cell array of FRD objects used to evaluate passivity
%   Z2_array  - cell array of FRD objects to plot impedance magnitude/imag
%   n_o, n_i  - output/input indices (for selecting matrix element)
%   opts      - optional struct: fields 'plot_pass' (bool), 'plot_type' ('imag')
%
% This utility plots passivity (min eigenvalue of Hermitian part) and
% either |Z| or Im(Z) for several parameter sweeps contained in Z2_array.

figure(f);
f_size = 16;

n_p = length(eig_array);
color_plot = colormap(jet(max(1,n_p)));


for i = 1:n_p
    color = color_plot(i,:);
    eig_n = eig_array{1,i};

    p = plot(real(eig_n), imag(eig_n),'o', 'Color', color);
    hold on
end


grid on
xlabel('Real', 'interpreter', 'latex');
ylabel('Imaginary');
ax = gca;
ax.FontSize = f_size;


end

