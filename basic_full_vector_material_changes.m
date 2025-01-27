% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2_initial = 3.305; % Initial ridge index
n2_final = 3.44;    % Final ridge index
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw_initial = 0.325; % Initial Ridge half-width
rw_final = 1.0;     % Final Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

% Make the mesh 8 times less dense:
dx_coarse = dx * 8;
dy_coarse = dy * 8;

lambda = 1.55;      % vacuum wavelength
nmodes = 10;         % number of modes to compute

% Ridge index sweep:
n2_values = linspace(n2_initial, n2_final, 10);
neff_results = zeros(10, nmodes); % Storing the effective indices for all steps

figure;
for i = 1:length(n2_values)
    n2 = n2_values(i); % Current ridge index

    % Generate waveguide mesh for the current ridge index
    [x, y, xc, yc, nx, ny, eps, edges] = waveguidemesh([n1, n2, n3], ...
                                                        [h1, h2, h3], ...
                                                        rh, rw_initial, side, ...
                                                        dx_coarse, dy_coarse);

    % Compute TE modes:
    [Hx, Hy, neff] = wgmodes(lambda, n2, nmodes, dx_coarse, dy_coarse, eps, '000A');
    neff_results(i, :) = neff; % Store neff for this ridge index

    % Plot modes for the first mode (example visualization)
    subplot(2, 5, i);
    contourmode(x, y, Hx(:, :, 1)); % Plot Hx of the first TE mode
    title(['n2 = ' num2str(n2)]);
    xlabel('x'); ylabel('y');
    for v = edges, line(v{:}); end
end

% Plot effective indices (neff) vs ridge refractive index
figure;
plot(n2_values, neff_results(:, 1), '-o');
title('Effective Index (neff) vs Ridge Refractive Index (n2)');
xlabel('Ridge Index (n2)');
ylabel('Effective Index (neff)');
grid on;



