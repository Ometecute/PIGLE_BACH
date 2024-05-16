%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% define scattering parameters and compute the wavevector components for each scattering azimuth
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% Incident energy /meV
E = 25;

% Probe particle mass /amu
m = 4;

% Incident angle /rad
theta_i = 0;

% BEAM PROPAGATION:
[k_iz, k_Gz] = beamprops_1(E, m, theta_i, Gx); % We pass in Gx here b/c all Gy(i)=0

% area of unit cells considered
S = nCells^2 * norm(cross([a1 0], [a2 0]));



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% diffraction intensity calculation
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

for n=1:length(Gx)
    % amplitude from single scattering center, j, for a given G:
    A_j = exp(1i*(Gx(n)*x+Gy(n)*y+(k_Gz(n)-k_iz)*zeta));

    % Integrate over unit cell to obtain total amplitude for the given G:
    A = 1/S * trapz(y(:,1),trapz(x(1,:),A_j,2));

    % Total intensity for given G:
    P(n) = abs(k_Gz(n)/k_iz) * abs(A)^2;
end



%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% check validity of model and plot intensity distribution
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% Sum of all intensities
U = sum(P, 'all')

% Bar graph plot
fig3 = figure; hold on
bar(Gx,P)