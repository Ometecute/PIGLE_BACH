% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

% describe the job
job_name = 'TempSweep';
proj_name = 'Pyrazine'; % Folder under which MD results will be saved

% Enable and config distributed computing (for using a job scheduler such
% as 'slurm')
distributed_computing = 0;
dist_comp_cmd = 'sbatch ~/pigle_sweep_slurm_submit.peta4-skylake';
n_cores_available = 6; %  How many parallel processes can be used on the target machine

% Define the parameter space and values to be probed.
sweepParams = {'T'};
sweepVal = {{ 230   240   250   260   270   280   290   300}};
