# Single Gate Multi Level Beam Switching

## Design code of single-gate multi-level beam switching [metasurface](https://en.wikipedia.org/wiki/Electromagnetic_metasurface) utilizing [Adjoint method](https://en.wikipedia.org/wiki/Adjoint_state_method)


This document provides a design code for single-gate multi-level beam switching based on a metasurface, employing the adjoint method for optimization. This code is associated with the research paper [High-directivity multi-level beam switching with single-gate tunable metasurfaces based on graphene](https://arxiv.org/abs/2410.00806)


## Figure of merit
<!-- ![plot](./images/env.png) -->
The design aims to maximize both directivity and target diffraction efficiency at specific wavelengths and gate voltages, which are customizable by the user.


## Algorithm
The core of the design utilizes the adjoint method, an efficient gradient calculation technique. The adjoint method offers a significant advantage over [finite-difference](https://en.wikipedia.org/wiki/Finite_difference) approaches by providing gradient values on design variables irrespective of their quantity, thereby accelerating the optimization process. Detailed information regarding the application of the adjoint method in nanophotonic systems can be found in our [review paper](https://www.degruyter.com/document/doi/10.1515/nanoph-2021-0713/html?lang=en&srsltid=AfmBOoopnvQKaBim4U1x62GNLuUmxDV_sV_W0sMarN0pgrkE5Q17UaBR).


## Simulation software
This work employs the open-source MATLAB-based Rigorous Coupled-Wave Analysis [RCWA](https://en.wikipedia.org/wiki/Rigorous_coupled-wave_analysis) software, [RETICOLO] (https://zenodo.org/record/3610175#.YBkECS2UGX0). In this repository we already included RETICOLO V9 version for convenience.


## Installation and running the optimization
1. Clone the repository using the following command:
To install, clone this repository in your desired folder:
```
git clone --https://github.com/jhpark94/Single_Gate_Multi_Level_Beam_Switching.git
```

2. Execute the MATLAB script 'Optimization_multi_level_beam_switch.m'. Optimization parameters can be adjusted within the script. Example parameters include:
```
%% Optimization parameters
folderindex = 8000; % save folder index
OptParm.wavelength = 8.0e-6; % wavelength used in geometrical parameter calculation
OptParm.sim_wavelength = 8.0e-6; % wavelength of input wave (useful for broadband spectrum analysis)
OptParm.N = 20; % # of design variables
OptParm.EFs = [0.2, 0.6, 1.0]; % Target Fermi levels
OptParm.diffraction_channels = [1 0 -1]; % Propagating diffraction channels
OptParm.gradient_type = 'shape-derivative'; % Either 'grayscale' or 'shape-derivative'
OptParm.angle = 80; % Target diffraction angle
heights = 1.5*pi; 
thicknesses = 1.5*pi;
repeats = 1:30;
OptParm.aspect_ratio = 10; % maximum aspect ratio
OptParm.spacer_thickness = 30e-9; % SiN spacer thickness
OptParm.t_HfO2 = 50e-9; % HfO2 thickness
OptParm.b_coefficient = 0.3; % trade-off coefficient
```


## citation
If you utilize this code, please cite our paper: https://arxiv.org/abs/2410.00806. The citation will be updated upon journal acceptance.