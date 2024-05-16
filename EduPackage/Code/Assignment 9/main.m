% clears workspace variables
clear;

% prevents console output for fsolve()
o = optimset('Display','off');

%% --------------------------------------------------------------
% initialise scattering parameters
%---------------------------------------------------------------

% reduced Planck's constant units of J, kg and u
hbar = 2.05;

% Incident energy /meV
E = 35;

% Probe particle mass /amu
m = 4;

% Incident angle /rad
theta_i = 0;



%% --------------------------------------------------------------
% define the momentum transfers we're interested in
% ---------------------------------------------------------------

% number of moduli of transfers we wish to consider
nDK = 500;

% set of moduli in given range
mods_DK = linspace(-20,20,nDK);

% set of considered directions
args_DK = 0;

% final set of considered momentum transfers as vectors
DKs = [cos(args_DK); sin(args_DK)] * mods_DK;



%% --------------------------------------------------------------
% beam propagation for each considered momentum transfer
% ---------------------------------------------------------------

% z-components of wavevectors /A^{-1}
[k_iz,k_Gz] = beamprops_3(E,m,theta_i,DKs);



%% --------------------------------------------------------------
% lattice simulation parameters
% ---------------------------------------------------------------

% number of unit cells per side (x and y) of the considered surface
len = 10;

% maximum allowed vertical height of stacked adsorbates
hgt = 1;

% substrate lattice parameter /A
a = 3.912;

% adsorbate adlayer depth /A
c = 1;

% temperature /K
T = 300;

% number of adsorbates
na = 40;

% number of timesteps between each adsorbate
nt = 5;


%% --------------------------------------------------------------
% run simulation
% ---------------------------------------------------------------

% intialise grid
grid = Grid(a, c, T, len, hgt);

% iteratively add adsorbates, giving time (~nt ps) for them to diffuse between each one
for ia =1:na
    
    % add a particle to the surface @ a random location
	grid.addRandParticle();
	
	% allow nt timesteps for the particles to approach equilibrium
	for it = 1:nt
        
        % perform one timestep
        grid.timestep();
	end
end

% Optionally plot positions of adsorbates
% grid.plot();

%% --------------------------------------------------------------
% extract grid of surface
% ---------------------------------------------------------------

% grid of the orthogonal (x,y,z) positions of possible adsorbate sites
R = grid.Rxyz;

% all x-coords
X = R(:, 1);

% all y-coorda
Y = R(:, 2);

% all z-coords
Z = R(:, 3);

% extract index (u,v,w) positions of adsorbates and reshape ( reshaping isn't strictly necessary due to something in Matlab called 'linear indexing', but is done here for clarity)
rs = reshape(grid.rs, [len^2*hgt 1]);



%% --------------------------------------------------------------
% helium/platinum potential
%---------------------------------------------------------------

% parameters

% multiplicative factor, parameter /a.u.
D_e = 2.89*10^(-4);

% exponential prefactor, parameter /a.u.
alpha = 0.52;

% point where Vs(z_m)=0, parameter /a.u.
z_m = 11.46;


% potential for He -> Pt surface
Vs = @(z) D_e*(exp(-2*alpha*(z-z_m)) - exp(-alpha*(z-z_m)));



%% --------------------------------------------------------------
% helium/single silver adatom potential
% ---------------------------------------------------------------

% parameters

% multiplicative factor, parameter /a.u.
C_12 = 6.33 * 10^6;

% multiplicative factor, parameter /a.u.
C_6 = 29.0;

% potential function for He -> single Ag adsorbate
V_He_Ag = @(x0,y0,z0,x,y,z) C_12./((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^6 + C_6./((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^3;



%% --------------------------------------------------------------
% helium/all silver adatoms potential
% ---------------------------------------------------------------

% Sum up all adsorbate potentials:

% initialise Va function
Va = @(x,y,z) 0;

% for each index position in the grid...
for uvw = 1:len^2*hgt

% if there's a particle at (u,v,w)
	if rs(uvw)==1

        % add a contribution to the adsorbate potential with its origin at this position
        Va = @(x,y,z) Va(x,y,z) + V_He_Ag(X(uvw),Y(uvw),Z(uvw),x,y,z);
    end
end



%% ---------------------------------------------------------
% calculate final potential
% ---------------------------------------------------------------

% Redefine meshgrid of surface to ignore the z-dimension:

% exclude z-dimension since z is now a variable
R = R(1:len^2,1:2);

% redefine X,Y to ignore permutations over possible z-coords
X = R(:,1); Y = R(:,2);

% total configuration potential for all adsorbates + surface
Vf = @(z) Vs(z) + Va(X,Y,z);



%% ---------------------------------------------------------
% Solve energy-conservation for corrugation
% ---------------------------------------------------------------

% energy conservation equation
eqn = @(z) Vf(z) - (hbar * k_iz)^2/(2*m);

% initial guess (value far from adsorbates on preliminary tests)
z0 = 0.25*ones([len^2 1]);

% solve numerically for corrugation
zeta = fsolve(eqn, z0);

% grid area /A^2
S = sqrt(3)/2 * (len*a)^2;



%% ---------------------------------------------------------
% calculate intensity for each momentum transfer
% ---------------------------------------------------------------

% Eikonal + Sudden approximation formula for integrand
phi = reshape(exp(1i*(-zeta*(k_Gz - k_iz) + R*DKs)), [len len nDK]);

% Integrate for intensity and apply kinematic correction
I = abs(k_Gz'/k_iz) .* reshape(abs(trapz(Y(1:len), trapz(X(1:len), phi, 1), 2)).^2 ./S^2, [nDK 1]);

% plot graph
figure();
plot(mods_DK,I);