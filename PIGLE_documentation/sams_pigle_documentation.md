---
useMath: true
mainfont: Linux Biolinum O
fontsize: 10pt
geometry: margin=3.7cm
link-citations: true
colorlinks: true
title: PIGLE Starter Guide
author: SM Lambrick
date: August 2022
---

This guide goes through the steps of installing then getting started running PIGLE (Particles Interacting in Generalized Langevin Equation simulator). It is a practical instruction set on using the software and does not go into the theory or the science of the method, refer to the [associated publication](https://doi.org/10.1016/j.cpc.2019.04.013) for information on the method.

# What is PIGLE

PIGLE (Particles Interacting in Generalized Langevin Equation simulator) is a simulator of adsorbate particle dynamics on surfaces. The adsorbed particles may move in 4 dimensions: *x*, *y*, *z*, and rotations. The simulation of the *z* and rotational motion is optional. In addition adsorbate-adsorbate interactions may be included. The simulator models the *Generalized Langevin Equation*,

$$
m\ddot{\mathbf{x}}_i = -\nabla V[\mathbf{x}_i(t)] - m\int_{-\infty}^t \gamma_i(t - t^\prime)\dot{\mathbf{x}}_i(t^\prime) \mathrm{d}t^\prime + \zeta(t) + \sum_{j\neq i}F_{ij},
$$

which includes a time dependent friction term in addition to the [ordinary Langevin Equation](https://en.wikipedia.org/wiki/Langevin_equation).

## Units

In general the units of PIGLE are:

- ps for time (ps<sup>-1</sup> for inverse time)  
- Å (Angstroms) for spatial distances (Å<sup>-1</sup> for inverse distances)  
- meV for energy

# Downloading/installing

## Prerequisits

PIGLE uses Matlab and Simulik therefore a local install of both is needed. Refer to the [Mathworks website](www.mathworks.com) for instruction on the installation of Matlab and Simulink. PIGLE can be run as a job on a high performance computing unit, however that is not covered in this documentation.

The following toolboxes are needed in addition to the base install of Matlab+Simulink:

- Parallel Computing Toolbox
- Polyspace Bug Finder
- Simulink Coder
- Simulink Compiler
- DSP system toolbox
- (on linux) a supported compiler such as gcc
- (on Windows) the [MinGW-w64 Compiler toolbox](https://uk.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler)

### Operating system

"In theory" PIGLE "should" work on any operating system [officially supported by Matlab/Simulink](https://uk.mathworks.com/support/requirements/matlab-system-requirements.html). In practice I have successfully run PIGLE on

- Windows 10 Pro (I assume Windows 11 will also work) with Matlab R2022a
- Ubuntu 22.04 LTS (not officially supported but does work) Matlab R2022a

Older version of Matlab down to at least R2017b should work but have not been tested.

Officially "supported" Ubuntu 20.04 has issues out of the box when running Simulink due to incompatibilities between compiler versions. Unfortunately I do not have access to a MacOS system to test.

## Getting PIGLE

The permanent record of PIGLE is stored in a [Zenodo repository](http://doi.org/10.5281/zenodo.2025809) while the up-to-date development version is kept on the [Cambrideg Atom Scattering Centre github](https://github.com/Cambridge-Atom-Scattering-Centre/PIGLE). At present it is recommended to use the version on GitHub and to use the branch `sam`. Click on the dropdown menu to see the list of available branches to choose `sam`.

![Changing branch on GitHub. It is recommended to use the 'sam' branch.](change_branch.png)

Once on the correct branch the code may be downloaded by clicking the "Code" dropdown menu and selecting "download zip". Alternativly `git` may be used to download the code by running the command  
`git clone https://github.com/Cambridge-Atom-Scattering-Centre/PIGLE.git`  
followed by  
`git checkout sam`

# Running PIGLE

The basic steps in running PIGLE are:

1. (optional) Generate a potential and save it as a `mat` file.
2. Configure the parameters for your adsorbates and surface in the `.m` files under the subfolder "UI"
3. Run the main script `run_pigle.m`
4. (optional) Use additional scripts for plotting and analysis.

When you first download PIGLE the simulation will be set up for the motion of CO molecules on a Cu(111) surface that are measured by the Cambridge He3 Spine-Echo machine. Therefore you can test run PIGLE before changing the parameters for your system (see setting parameters below).

# Setting parameters

With the exception of the potential energy surface the parameters for a PIGLE simulation are specified in three `.m` files found in the "UI" subfolder. Modify the parameters by changing the value of Matlab variables in these files.

There are many parameters that are 'true/false' or 'on/off' options, in general these use the C convention of `1` being true/on and `0` being false/off.

List of files including parameters:

| File | Purpose |
|------|---------|
| `prep_environment.m` | Defines paths for PIGLE and data saving |
| `pigle_wrapper_params.m`    | Non-scientific general parameters, e.g. plotting. |
| `pigle_ui.m`                | Key scientific simulation parameters. |
| `pigle_ui_surface_params.m` | Parameters for the surface and the adsorbates. |
| `params_for_function_ prepare_potential.m`| Specifying a potential for PIGLE to interpolate. |
| `prepare_params_for_ interactions.m` | Define adsorbate-adsorbate interactions. |
| `params_for_pos_depended_ spatial_friction.m` | Can be used to create a spatial dependence of the friction. |
| `params_for_pos_depended_ spatial_theta_friction.m`| Can be used to create a spatial dependence of the rotational friction. |
| `config_job_params.m` | HPC use parameters. |

## Data path

The current process for setting the data path is convoluted and I haven't fully got my head round the best way to set your own data path, however when `run_pigle.m` is run text appears at the command prompt stating where the data file is being saved.

## Wrapper parameters
`pigle_wrapper_params.m`

The wrapper parameters are related to the running and output of the simultion, they do not change what is actually simulated.

| Parameter  | Possible values | Description |
|------------|-----------------|-------------|
| `isISF`    | 1/0             | should an ISF be calcuated? |
| `isSave`   | 1/0             | should the results be saved? |
| `ISF2save` | list, 1/2       | save the 1-incoherent ISF, 2-coherent |
| `toPlot`   | 1/0             | should plots of the results be made |
|  `reduceData` | 2/1/0        | reduce the amount of data saved (0=keep all) |  
|  `clearParams` | 1/0         | clear previous parameters file |

## Simulation parameters
`pigle_ui.m`

Here the main simulation parameters are defined.

### Overall simulation parameters

The first section of parameters define the overall simulation to be performed, some of these are scientific and some are practical. Key *scientific* options are should *rotations be included*, should *z motion be included*, and the *initial momentum* of the adsorbate (0 or thermal).

| Parameter      | Possible values | Description |
|----------------|-----------------|-------------|
| `z_enabled`    | 1/0 | include z motion of adsorbates? |  
| `dKz_include_in_isf` | 1/0 | calculate a $\Delta K$ perpendicular ISF? |
| `theta_enables` | 1/0 | include adsorbate rotations? |
| `zero_p_init`  | 1/0 | set the initial momentum to be 0 or thermal |
| `N_runs`       | +ve integer | how many runs of the simulation to perform |
| `run_parallel` | 1/0 | should parallel computing be used? |

### Delta K & azimuth

The next set of parameters specify the the $\Delta K$s to be use in the simulation and two azimuths to simulate. The $\Delta K$s are specified in a Matlab array of all desired $\Delta K$ in units of Å<sup>-1</sup> and the azimuths are defined by 2 element arrays that define azimuthal directions using the crystolographic directions.

### Beam

PIGLE is initially set up to simulate the Cambridge He3 Spin-Echo, therefore the beam parameters, total scattering angle and beam incidence wavevector, only need to be changed if you want to simulate a different machine.

### Time and simulation steps

The time and simulation steps may be set via different combinations of variables, either a total time + a time step time, or a time step time + a total number of time steps. The smaller value of the two options will be used by PIGLE.

| Parameter | Description |
|-----------|-------------|
| `sample_time` | Time step for the simulation, ps |
|`sample_time_clist` | ??? |
| `isf_sample_time` | ISF time interval, ps |
|`thermalizing_time` | Time ignored at the beginning of the simulation to allow thermalisation. |
|`stop_time` | Total simulation time |
| `max_N_steps` | Total maximum allowed steps |
| `max_N_ISF_steps` | Total maximum allowed stps in ISF |


## Surface parameters
`pigle_ui_surface_params.m`  
`params_for_function_prepare_potential.m`

*Throughout PIGLE the initialism PES is used to refer to a potential energy surface.*

There are two ways to specify the adsorbate-surface potential in PIGLE: either the potential may be generated externally, saved as a `.mat` file then imported into the correct structure; or heights of specific points on the potential surface may be specified and then PIGLE will interpolate using Fourier components in between them. A full description


### PIGLE interpolated PES

Built into PIGLE is an interpolation function to generate a PES for a close packed hexagonal surface. Four heights and two slopes are needed to specify a close packed surface, they are demonstrated on the figure below.

![The locations of the six potential values needed for PIGLE to interpolate a PES.](potential_specification.png)

These six parameters need to be specified as a column vector in the file
`params_for_function_prepare_potential.m` as `pot_strct(i).V`, where `pot_strct(i)` is the a Matlab structure containing all the information on the potential, `i` means we are setting the potential for the *i*th species on the surface, and `.V` means we are specifying the potential values.

If *z* motion is being considered then the parameters for a Morse potential,  
$$
V(r) = D_e\left[1 - e^(-a(r-r_e))\right]^2,
$$  
need to be specified at each point as well as the potential height. In the specification of the potential all **potential values are in meV** and all **position values are in Angstroms**.

The total code for specifying a potential is shown below,

```
pot_strct(i).ref_De = [top; slope1; slop2; bridge; hcp; fcc];
pot_strct(i).V = [top; slope1; slop2; bridge; hcp; fcc];
pot_strct(i).a = [top; slope1; slop2; bridge; hcp; fcc];
pot_strct(i).r_e = [top; slope1; slop2; bridge; hcp; fcc];

pot_strct(i).is_potval = [1 0 0 1 1 1];
pot_strct(i).f_2D = @hexagonal6interp;
```

`hexagonal6interp` is the function that performs Fourier interpolation between the given potential values to produce a PES and `is_potval` should not need to be edited.

If more than one adsorbate species is being simulated then PES need to be defined for every adsorbate species and added to the array of potentials `pot_strct`.

A very similar specification to that used for the PES is used to define position dependent friction (spatial and theta) in the files:  
`params_for_pos_depended_spatial_friction.m`  
`params_for_pos_depended_spatial_theta_friction.m`    
A struct called `friction_strct(i)` is constructed for each species on the surface with the same set of member variables: `ref_De`, `V`, `a`, `r_e`, `f_2D`, and `is_potval`.

### Specifying your own PES

...

# A PIGLE exercise

This is a very short exercise in running a PIGLE simulation and looking at the data. A more detailed exercise in PIGLE can be found as assignments 13 and 14 part of the Cambridge Atom Scattering Centre [Educational Package](https://atomscattering.phy.cam.ac.uk/resources/educational-package/), alternatively you may wish to get on with simulating your system after this brief introduction.
