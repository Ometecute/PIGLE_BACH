% avoids a cluttered workspace
clear;



%% --------------------------------------------------------
% define simulation parameters
% ---------------------------------------------------------

% hop probability
p = 0.01;

% lattice parameter /A
a = 2.5;

% number of timesteps in trajectory simulation
n_t = 100000;

% length of timesteps /micro s
dt = 1;



%% --------------------------------------------------------
% define set of momentum transfers we are interested in
% ---------------------------------------------------------

% number of different moduli DK we wish to find the ISF for
n_DK = 100;

% set of those moduli in given range /A^{-1}
mods_DK = linspace(0, 6, n_DK);

% direction of DK /rad
arg_DK = 0;

% unit vector in given direction
unitvec = [cos(arg_DK) sin(arg_DK)];

% scaled set of vectors with direction arg_DK and moduli mods_DK
DKs = mods_DK' * unitvec;



%% --------------------------------------------------------
% iterate over each momentum transfer to calculate a fitteddephasing rate alpha(DK) for each one
% ---------------------------------------------------------

for j=1:n_DK
    
    DK = DKs(j,:);

    % ---------------------------------------------------------
    % run single MC simulation and calculate scattering amplitude
    % ---------------------------------------------------------

    % simulate
    mcs = MCSimu_sqr(p, a, dt, n_t);
    mcs.nStates_loop(n_t);

    % extract trajectory
    x = mcs.x; y = mcs.y; t = mcs.t;

    % calculate scattering amplitude
    A = exp(-1i*(DK(1)*x + DK(2)*y));
    
    
    
    % ---------------------------------------------------------
    % ISF calculation
    % ---------------------------------------------------------

    % evaluate ISF formula
    step_1 = fft(A);
    step_2 = abs(step_1).^2;
    I = ifft(step_2);
    
    % half timebase
    t = t(1:floor(n_t/2)); I = I(1:floor(n_t/2));

    % take normalised real part
    ReI = real(I)/real(I(1));
    
    
    
    % ---------------------------------------------------------
    % perform exponential fit of ISF
    % ---------------------------------------------------------
    
    % attempt fit
    g = fittype('exp(-a*x)');
    fitg = fit(t, ReI, g, 'startpoint', [0.001]);
    
    % extract dephasing rate alpha for this momentum transfer DK
    alpha(j) = coeffvalues(fitg);
    
end

%% --------------------------------------------------------
% plot dephasing rate curve
% ---------------------------------------------------------

figure();
plot(mods_DK, alpha); hold on;
xlabel('parallel momentum transfer'); ylabel('dephasing rate');



%% --------------------------------------------------------
% fit curve to dingle-jump diffusion model
% ---------------------------------------------------------

% fit single-jump diffusion model
h = fittype('c*(1-cos(b*x))');
fith = fit(mods_DK', alpha', h, 'startpoint', [2*pi/a 0.01]);

% plot fitted curve
plot(fith); hold off;