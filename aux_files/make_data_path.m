
% Make data path - script to construct the data path for saving the
% simulation results

global proj_name
global data_path
global verbal

if ~(exist('CustomFuncs','file')==2)
    error('make_data_path:ERROR: CustomFuncs.m does not exist')
end
if ~(exist(CustomFuncs.find_data_path(proj_name),'dir')==7)
    error('make_data_path:ERROR: cannot find project data (directory) path')
end
    
data_path = 'results';%CustomFuncs.find_data_path(proj_name);
serial=define_datafile_mfile('md',0);

if verbal, disp(['serial num for file: ' num2str(serial)]); end
