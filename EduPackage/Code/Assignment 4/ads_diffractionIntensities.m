set(groot,'defaulttextinterpreter','latex');
set(groot,'DefaultTextFontname', 'CMU Serif');
set(groot,'DefaultAxesFontName', 'CMU Serif');
set(groot,'DefaultAxesTickLabelInterpreter','latex'); % Using '0' instead of 'groot' is depreciated
clear;

%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% scattering parameters
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% hexagonal plane lattice constant
a = 2.55;

% Incident energy /meV
E = 25;

% Probe particle mass /amu
m = 4;

% Incident angle /rad
theta_i = 0;

% Substrate corrugation depth /A
hs = 1;

% Adsorbate corrugation depth /A
ha = 1;



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% define reciprocal lattice vectors of interest and calculate corresponding wavevector components
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% define reciprocal lattice vectors for [1 0] azimuth of surface ( actually [110] in 3D lattice but we assume an in−plane basis):

% reciprocal point spacing along [10]
dGx = 2*pi/(a/2);

% set of first three x−components in either direction, plus zeroth
Gx = [-3*dGx:dGx:3*dGx];

% y−component zero along [10] azimuth
Gy = zeros(size(Gx));

% beam propagation
[k_iz, k_Gz] = beamprops_1(E, m, theta_i, Gx);



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% define unit cell vectors and reciprocal vectors for substrate and adsorbate layers
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% substrate vectors
[a1s, a2s, b1s, b2s] = auxCls.unitCellVecs(1,a,0);

% adsorbate layers > 3x thelattice parameter
[a1a, a2a, b1a, b2a] = auxCls.unitCellVecs(1,3*a,30);

% define third 'redundant' reciprocal lattice vectors for both substrate and adsorbate
b3s = b2s - b1s;
b3a = b2a - b1a;



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% grid supercell of the ensemble with origin at the centre
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% scaling factor for the number of cells included in the supercell
nCells = 7*sqrt(3);

% Range along x and y axes −− note, norm(a1+2*a2) is the perpendicular height of unit cell
Lx = nCells*[-a/2 a/2]; Ly = nCells*[-norm(a1s+2*a2s)/2 norm(a1s+2*a2s)/2];

% Linspaces across x and y axes
lx = linspace(Lx(1),Lx(2),30*nCells); lx(end) = [];
ly = linspace(Ly(1),Ly(2),50*nCells); ly(end) = [];

% meshgrid across range with desired resolution
[x,y] = meshgrid(lx,ly);
R = [x(:) y(:)];



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% evaluate corrugations for both substrate and adsorbate layers then resize to fit meshgrid formulation
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% individual corrugations
zeta_s = reshape(sample_zeta_2(R, hs, b1s', b2s', b3s'), size(x) );
zeta_a = reshape(sample_zeta_2(R, ha, b1a', b2a', b3a'), size(x));

% total corrugation
zeta = zeta_s + zeta_a;



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% diffraction intensity calculation
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% area of unit cell
S = nCells^2 * norm(cross([a1s 0], [a2s 0]));

for n=1:length(Gx)
    % Amplitude from single scattering centre
    A_j = exp(1i*(Gx(n)*x+Gy(n)*y+(k_Gz(n)-k_iz)*zeta));

    % Total amplitudes
    A = 1/S * trapz(y(:,1),trapz(x(1,:),A_j,2));

    % Intensities
    P(n) = abs(k_Gz(n)/k_iz) * abs(A)^2;
end



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% check validity of approximation then plot the intensity distribution
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% Sum of all intensities
U = sum(P, 'all');

% Bar graph plot
fig3 = figure; hold on
bar(Gx,P)
%axis ([−20 20 0 0.1]) % Useful if you wish to inspect specific regions of the graph