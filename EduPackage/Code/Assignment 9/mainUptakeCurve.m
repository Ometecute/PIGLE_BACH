% clears workspace variables
clear;

% prevents console output for fsolve()
o = optimset('Display','off');

%% --------------------------------------------------------
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



%% --------------------------------------------------------
% specular beam => only consider DK=0 (here we allow for study of other beams)
% --------------------------------------------------------------

% number of momentum transfers of interest => only one here b/c we only consider specular scatter
nDK = 1;

% specular scatter has |DK|=0
mods_DK = 0;

% argument is undefined for specular scatter
args_DK = 0;

% only one vector => DKs = [0; 0], but useful to keep the option to generalise
DKs = [cos(args_DK); sin(args_DK)] * mods_DK;



%% --------------------------------------------------------
% beam propagation for specular beam only
% --------------------------------------------------------------

% Incident wavevector, z-component /A^{-1}
[k_iz,k_Gz] = beamprops_3(E,m,theta_i,DKs);



%% --------------------------------------------------------
% lattice simulation parameters
% --------------------------------------------------------------

% number of unit cells per side (x and y) of the considered surface
len = 5;

% maximum allowed vertical height of stacked adsorbates
hgt = 7;

% lattice parameter /A
a = 3.912;

% adsorbate adlayer depth /A
c = 1;

% temperature /K
T = 300;

% number of adsorbates
na = 11;
% number of timesteps between each adsorbate
nt = 5;



%% -------------------------------------------------------------
% define parameters for looping over configurations
% --------------------------------------------------------------

% number of confiugurations
nbeta = 10;

% matrix of 'intensity' columns for each configuration
I_beta = zeros([na+1 nbeta]);

% matrix of 'no. adsorbate' columns for each configuration
nas_beta = zeros([na+1 nbeta]);



%% -------------------------------------------------------------
% loop over nbeta configurations
% --------------------------------------------------------------

% iterate over nbeta configurations of adding na adparticles randomly to the surface
for ibeta = 1:nbeta

    % initialise clean grid
	grid = Grid(a, c, T, len, hgt);
	
	% specular intensities
	I = zeros([na+1 1]);
	
	% corresponding number of adsorbates
    nas = zeros([na+1 1]);
    
    %-----------------------------------------------------------
    % slowly add adsorbates to the surface, allow them to approach equilibrium with a few timesteps, then calculate spectal intensity
	%-----------------------------------------------------------
    
    % iteratively add adsorbates, giving time (~nt ps) for them to diffuse between each one. For each new adsorbate added, recalculate the spectral intensity
	for ia =1:na+1
        
        % extract grid of surface
        
        % grid of the orthogonal (x,y,z) positions of possible adsorbate sites
        R = grid.Rxyz;
        
        % all x-coords
        X = R(:, 1);
    	
        % all y-coorda
        Y = R(:, 2);
        
        % all z-coords
        Z = R(:, 3);
        
        % extract index (u,v,w) positions of adsorbates and reshape (reshaping isn't strictly necessary due to something in Matlab called 'linear indexing', but is done here for clarity)
        rs = reshape(grid.rs, [len^2*hgt 1]);
        
        % calculate the number of adsorbates at this iteration (not necessarily equal to ia since adsorbate placement may fail)
        nas(ia) = sum(rs);
        
        %-------------------------------------------------------
        % helium/platinum potential
        %-------------------------------------------------------
        
        % multiplicative factor, parameter /a.u.
        D_e = 2.89*10^(-4);
        
        % exponential prefactor, parameter /a.u.
        alpha = 0.52;
        
        % point where Vs(z_m)=0, parameter /a.u.
        z_m = 11.46;
        
        
        % potential for He -> Pt surface
        Vs = @(z) D_e*(exp(-2*alpha*(z-z_m)) - exp(-alpha*(z-z_m)));
        
        
        
        %-------------------------------------------------------
        % helium/single silver adatom potential
        % ------------------------------------------------------
        
        % multiplicative factor, parameter /a.u.
        C_12 = 6.33 * 10^6;
        
        % multiplicative factor, parameter /a.u.
        C_6 = 29.0;
        
        
        % potential function for He -> single Ag adsorbate
        V_He_Ag = @(x0,y0,z0,x,y,z) C_12./((x-x0).^2 + (y-y0).^2+ (z-z0).^2).^6 + C_6./((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^3;
        
        
        
        %-------------------------------------------------------
        % helium/all silver adatoms potential
        %-------------------------------------------------------
        
        % Sum up all adsorbate potentials:
        
        % initialise Va function
        Va = @(x,y,z) 0;
        
        % for each index position in the grid...
        for uvw = 1:len^2*hgt
            % if there's a particle at (u,v,w)
            if rs(uvw)>0
                % add a contribution to the adsorbate potential with its origin at this position
                Va = @(x,y,z) Va(x,y,z) + V_He_Ag(X(uvw),Y(uvw),Z(uvw),x,y,z);
            end
        end
        
        
        
        %-------------------------------------------------------
        % calculate final potential
        %-------------------------------------------------------
        
        % Redefine meshgrid of surface to ignore the z-dimension:
        
        % exclude z-dimension since z is now a variable
        R = R(1:len^2,1:2);
        
        % redefine X,Y to ignore permutations over possible z-coords
        X = R(:,1); Y = R(:,2);
        
        % total configuration potential for all adsorbates + surface
    	Vf = @(z) Vs(z) + Va(X,Y,z);
    	
        
        
        % ------------------------------------------------------
        % Solve energy-conservation for corrugation
        % ------------------------------------------------------
        
        % energy conservation equation
        eqn = @(z) Vf(z) - (hbar * k_iz)^2/(2*m);
        
        % initial estimate of turning point (value with no adsorbates in preliminary tests)
        z0 = 0.25*ones([len^2 1]);
        
        % solve numerically for corrugation
        zeta = fsolve(eqn, z0);
        
        % grid area /A^2
        S = sqrt(3)/2 * (len*a)^2;
        
        
        
        % ------------------------------------------------------
        % calculate intensity using Eikonal+Sudden approximation
        % ------------------------------------------------------
        
        % eikonal approximation for integrand
        phi = reshape(exp(1i*(-zeta*(k_Gz - k_iz) + R*DKs)), [len len nDK]);
        
        % integrate to get intensity (kinematic factor = 1 for specular beam)
        I(ia) = reshape(abs(trapz(Y(1:len), trapz(X(1:len), phi,1), 2)).^2 ./S^2, [nDK 1]);
        
        % add a particle to the surface @ a random location
        grid.addRandParticle();

        % allow nt timesteps for the particles to approach equilibrium
        for it = 1:nt

            % perform one timestep
            grid.timestep();

        end
    end
    % append list of intensities at each iteration from this configuration to a matrix (I_beta)
    I_beta(:,ibeta) = I;

    % append list of adsorbate numbers from this configuration to a matrix (nas_beta)
    nas_beta(:,ibeta) = nas;
end



%% -------------------------------------------------------------
% normalisation and plotting
% --------------------------------------------------------------

% average nas_beta and I_beta over all configurations <what do these represent?>
nas_avg = mean(nas_beta,2);
I_avg = mean(I_beta,2);

% convert nas to theta by dividing by the number of unit cells considered
theta = nas_avg./len^2;

% convert intensity to reflectivity by dividing by intensity on clean surface
Ref = I_avg/I_avg(1);

% plot graph
figure();
plot(theta, Ref);



%% -------------------------------------------------------------
% fitting analytical models for low coverages
% --------------------------------------------------------------

hold on;

% low-coverage limit

% fitting length
fitlen = floor(length(theta)/2);

% define equation of fit
f = fittype('exp(-a*x)');

% approximate parameter a in fit
fit_f = fit(theta, Ref, f);

% plot fitted equation for low-coverage
plot(fit_f);


% lattice gas model

% define equation of fit
g = fittype('(1-x)^a');

% approximate parameter a in fit
fit_g = fit(theta, Ref, g);

% plot fitted equation for lattice gas model
plot(fit_g);
% label axes
xlabel('Coverage $\Theta$');
ylabel('Specular Reflectivity $I/I_0$');


hold off;