%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% define lattice parameters and compute lattice vectors
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% Lattice parameter /A (Angstroms)
a = 4.52;

% Corrugation depth /A (Angstroms)
h = 1;

% define (hexagonal) unit cell vectors for surface
[a1, a2, b1, b2] = auxCls.unitCellVecs(1,a,0);
b3 = b2 - b1;

%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% grid 2x2 unit cells and evaluate the potential across it
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% a scaling factor for the number of cells included in the supercell
nCells = 2;

% Range along x and y axes −− note, norm(a1+2*a2) is the perpendicular height of unit cell
Lx = nCells*[-a/2 a/2]; Ly = nCells*[-norm(a1+2*a2)/2 norm(a1+2* a2)/2];

% Linspaces across x and y axes
lx = linspace(Lx(1),Lx(2),30*nCells); lx(end) = [];
ly = linspace(Ly(1),Ly(2),50*nCells); ly(end) = [];

[x,y] = meshgrid(lx,ly);
R = [x(:) y(:)];

% evaluate corrugation at at each point on the real space grid
zeta = reshape(sample_zeta_2(R, h, b1', b2', b3'), size(x));

%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% plot the scattering planes on top of the contour plot of the corrugation
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% x−coords of first 3 planes from origin.
planePos = [-a/2 0 a/2]';

% PLOT:
fig2 = figure; hold on;

% Plot corrugation from tutorial 2.1.
pcolor(x,y,zeta); shading flat; axis equal

% Plot first 3 planes as projections on x−y plane (of any length

p = plot([planePos planePos]', repmat([-5 5]', 1, length(planePos)), 'color', '#800000', 'LineWidth', 1.5);

%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% define the reciprocal lattice vectors of interest (along the [1 0] direction)
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% separation of first 3 reciprocal vectors in k−space.
dGx = 2*pi/(a/2);

% x− and y−components of first three [1 0] reciprocal lattice vectors in both direcs
Gx = [-3*dGx:dGx:3*dGx];
Gy = zeros(size(Gx));

% plot the [1 0] azimuth over the above range
g = plot(Gx, Gy, '--', 'color', '#483493', 'LineWidth', 1.5);
legend([p(1) g], 'scattering planes', '<1 0> azimuth');