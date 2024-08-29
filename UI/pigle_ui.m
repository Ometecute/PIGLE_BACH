

%% parameters for config_model.m

z_enabled = 0;          % Should z motion be enabled
dKz_include_in_isf = 0; % Calculated the dK perpendicular ISF
theta_enabled = 0;      % Should rotational motion of the adsorbates be included
zero_p_init = 0;        % set initial momentum be set to zero? (if set to 0, p_init will correspond to thermal distribution)
interactions_active = 1;% Should rotational motion of the adsorbates be included
N_runs = 1;             % How many runs of the simultion to perform 
run_parallel = 1;       % Use parralell computing?

% Specify dK as a 2D vector, 3rd dim is azimuths.
dK = [0.2:0.1:1 1.2:0.2:2 2.5:0.25:5];
azim_1 = [1 0];
azim_2 = [1 1];

% specify beam parameters and geometrical parameters for scatering calculations
theta_tot = 44.4; % Degrees
beam_ki = 3.3977; % Angstrom ^{-1} 

% Specify simulation time parameters
% (those will be adjusted by the program, see below if interested)
sample_time = 5e-4;       % time step for simulation, ps
sample_time_clist = 1e-3; % time step for results
isf_sample_time = 5e-2;   % ISF time interval, ps
thermalizing_time = 50;   % time at the begining of the simulation to allow thermalization
stop_time = 1024*10;     % total time of the simulation

% N_steps and N_ISF_steps are calculated after PIGLE adjusts the requested time parameters
max_N_steps = 1e8;        % Total number of allowed steps
max_N_ISF_steps = 6e5;    % Total number of allowed steps in the ISF

