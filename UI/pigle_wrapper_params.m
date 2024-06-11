
%% Wrapper parameters

isISF        = 1; % 0 - just run the simulation (single run), 1 - calculate an averaged ISF
isSave       = 1; % save results (on/off)
ISF2save     = [1]; % which ISF to save? 1-Incoherent, 2-Coherent (needs isSave to be turned-on)
toPlot       = 1; %
reduceData   = 2; % 0 - don't reduce. 1 - remove p and r_supercell. 2 - Leave only one trajectory for each population. 3 - Remove all trajectories from all populations
clearParams  = 0; % Clear the structure containing all the model-configuration parameters from orevious run
