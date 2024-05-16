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
% define single momentum transfer we are interested in
% ---------------------------------------------------------

% modulus of parallel momentum transfer /A^{-1}
mod_DK = 2;

% direction of DK /rad
arg_DK = 0;

% set of DK Vectors /A^{-1}
DK = mod_DK .* [cos(arg_DK) sin(arg_DK)]';



%% --------------------------------------------------------
% run MC simulation to generate a single trajectory
% --------------------------------------------------------------

% instantiate square simulation object
mcs = MCSimu_sqr(p, a, dt, n_t);

% run simulation
mcs.nStates_loop(n_t);

% extract trajectory
x = mcs.x; y = mcs.y; t = mcs.t;



%% --------------------------------------------------------
% calculate and plot scattering amplitude
% ---------------------------------------------------------

% calculate amplitude
A = exp(-1i*(DK(1)*x + DK(2)*y));

% plot real part of amplitude vs time
figure();
plot(t, real(A));
xlabel('t /microseconds'); ylabel('Re(A) /arb. units');



%% --------------------------------------------------------
% ISF calculation
% ---------------------------------------------------------

% evaluate formula
step_1 = fft(A);
step_2 = abs(step_1).^2;
I = ifft(step_2);

% redefine timebases
t = t(1:floor(n_t/2));
I = I(1:floor(n_t/2));

% take normalised real part
ReI = real(I)/real(I(1));



%% --------------------------------------------------------
% ISF plot
% ---------------------------------------------------------

figure();
plot(t, ReI); hold on;



%% --------------------------------------------------------
% fit ISF to exponential decay and plot
% ---------------------------------------------------------

% fit exponential decay
g = fittype('exp(-a*x)');
fit_g = fit(t, ReI, g, 'startpoint', [0.001]);

% dephasing rate extracted from fit
alpha = coeffvalues(fit_g);

% plot fitted curve
plot(t, exp(-alpha.*t)); % alternatively, write plot(fit_g)
xlabel('t /microseconds'); ylabel('Re(I) /arb. units'); hold off;
