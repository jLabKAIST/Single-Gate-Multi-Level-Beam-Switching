# Single Gate Multi Level Beam Switching

## Design code of single-gate multi-level beam switching [metasurface](https://en.wikipedia.org/wiki/Electromagnetic_metasurface) utilizing [Adjoint method](https://en.wikipedia.org/wiki/Adjoint_state_method)


This document provides a design code for single-gate multi-level beam switching based on a metasurface, employing the adjoint method for optimization. This code is associated with the research paper [High-Efficiency Multi-Level Beam Switching with Single-Gate Tunable Metasurfaces Based on Graphene](https://doi.org/10.1002/adom.202500236).


## Figure of merit
<!-- ![plot](./images/env.png) -->
The design aims to maximize both directivity and target diffraction efficiency at specific wavelengths and gate voltages, which are customizable by the user.


## Algorithm
The core of the design utilizes the adjoint method, an efficient gradient calculation technique. The adjoint method offers a significant advantage over [finite-difference](https://en.wikipedia.org/wiki/Finite_difference) approaches by providing gradient values on design variables irrespective of their quantity, thereby accelerating the optimization process. Detailed information regarding the application of the adjoint method in nanophotonic systems can be found in our [review paper](https://www.degruyter.com/document/doi/10.1515/nanoph-2021-0713/html?lang=en&srsltid=AfmBOoopnvQKaBim4U1x62GNLuUmxDV_sV_W0sMarN0pgrkE5Q17UaBR).


## Simulation software
This work employs the open-source MATLAB-based Rigorous Coupled-Wave Analysis [RCWA](https://en.wikipedia.org/wiki/Rigorous_coupled-wave_analysis) software, [RETICOLO](https://zenodo.org/record/3610175#.YBkECS2UGX0). In this repository we already included RETICOLO V9 version for convenience.


## Installation and running the optimization
1. Clone the repository using the following command:
To install, clone this repository in your desired folder:
```
git clone --https://github.com/jLabKAIST/Single-Gate-Multi-Level-Beam-Switching.git
```

2. Execute the MATLAB script 'Optimization_multi_level_beam_switch.m'. Optimization parameters can be adjusted within the script. Example parameters include:
```
%% Optimization parameters
folderindex = 8000; % save folder index
OptParm.wavelength = 8.0e-6; % wavelength used in geometrical parameter calculation
OptParm.sim_wavelength = 8.0e-6; % wavelength of input wave (useful for broadband spectrum analysis)
OptParm.N = 20; % # of design variables
% For two-level beam switching, set OptParm.EFs = [0.2, 0.6]
% For three-level beam switching, set OptParm.EFs = [0.2, 0.6, 1.0]
% For four-level beam switching, set OptParm.EFs = [0.05, 0.35, 0.65, 0.95]
OptParm.EFs = [0.2, 0.6, 1.0]; % Target Fermi levels
% For two- and three-level beam switching, 
% set OptParm.diffraction_channels = [1 0 -1];
% For four-level beam switching, 
% set OptParm.diffraction_channels = [2 1 0 -1 -2];
OptParm.diffraction_channels = [1 0 -1]; % Propagating diffraction channels
OptParm.gradient_type = 'shape-derivative'; % Either 'grayscale' or 'shape-derivative'
OptParm.angle = 80; % Target diffraction angle
% heights = linspace(1.5*pi,3.5*pi,9); % height sweep
heights = 1.5*pi; 
% thicknesses = linspace(1.5*pi,3.5*pi,9); % thickness sweep
thicknesses = 1.5*pi;
repeats = 1:30;
OptParm.aspect_ratio = 10; % maximum aspect ratio
OptParm.spacer_thickness = 30e-9; % HfO2 spacer thickness
OptParm.t_HfO2 = 50e-9; % HfO2 etch-stop layer thickness
OptParm.b_coefficient = 0.3; % trade-off coefficient
%%
```

## Optimized structures
In 'Optimized_Structure' folder, we have included optimized solutions for two-, three-level metasurfaces at 7.0 7.5, 8.0 (main result), 8.5, and 9.0 ums. For the four-level, optimization result for 8.0 um is included.

## Citation
If you utilize this code, please cite our paper: https://doi.org/10.1002/adom.202500236.
