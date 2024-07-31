% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

%% PIGLE_RUN_SINGLE_TASK is designed to be executed by pigle_shell.

% load specific variables for this single task. The file
% pigle_shell_params.m is generated by config_job.m for each case
pigle_shell_params

pigle_path = [sub_job_path '/../../'];

cd(pigle_path)
prep_environment
cd(sub_job_path)
run('run_pigle.m')

disp('done pigle')
