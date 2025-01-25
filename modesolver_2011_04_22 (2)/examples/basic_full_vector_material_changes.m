% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
% n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;         % number of modes to compute

% Ridge index sweep from 3.305 to 3.44 in 10 steps 
ridge_indices = linspace(0.325, 3.44, 10);
neff_values = zeros(size(ridge_indices)); % Stores the neff values

for i = 1:length(ridge_indices)
    n2 = ridge_indices(i); % Current ridge index
   
    % Generate the waveguide mesh
     [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
                                             rh,rw,side,dx,dy); 

    
    % First consider the fundamental TE mode
    [Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000A');
    neff_values(i) = neff; % Stores the effective index
    
    % Plot the mode profiles for the current ridge index 
    figure(i);
    subplot(121);
    contourmode(x,y,Hx(:,:,1));
    title(sprintf('Hx(TE mode, n2=%.3f)', n2))
    xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
   
    subplot(122);
    contourmode(x,y,Hy(:,:,1));
    title(sprintf('Hy(TE mode, n2=%.3f)', n2)); 
    xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
end 

% Plot neff as a function of ridge index 
figure;
plot(ridge_indices,neff_values, '-o', 'LineWidth', 2);
xlabel('Ridge Index (n2)');
ylabel('Effective Index (neff)');
title('Effective Index vs Ridge Index');
grid on;

% ------ Removed code from the basic_fullvector_singlemode --------
% Mesh 8 times less dense 

% dx_coarse = 8 * dx;
% dy_coarse = 8 * dy;
% 
% fprintf('\nCalculating with 8x coarser mesh:\n');
% rw = 1.0; % Testing at the widest ridge width 
% 
% % Generating waveguide mesh
% 
% [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3], [h1,h2,h3],rh,rw,side,dx_coarse,dy_coarse);
% [Hx_coarse, Hy_coarse, neff_coarse] = wgmodes(lambda, n2, nmodes, dx_coarse, dy_coarse, eps, '000A');
% 
% fprintf('Coarse mesh: neff = %.6f\n', neff_coarse);
% 
% figure;
% subplot(121);
% contourmode(x, y, Hx_coarse(:,:,1));
% title('Hx (TE mode, coarse mesh)');
% xlabel('x'); ylabel('y');
% for v = edges, line(v{:}); end
% 
% subplot(122);
% contourmode(x, y, Hy_coarse(:,:,1));
% title('Hy (TE mode, coarse mesh)');
% xlabel('x'); ylabel('y');
% for v = edges, line(v{:}); end