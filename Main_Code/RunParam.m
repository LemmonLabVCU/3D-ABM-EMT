function [Param] = RunParam( )
% RunParam.m Initialize constant parameter values being used in the model. 

% Variables for meshgrid and spheroid
Param.n           = 32;               % grid size (units)

% Initializing Cell Characteristics
Param.Csize       = 15;               % size of the cell (um)
Param.Crado       = 2.5;              % outer radius fo the spheroid (units)
Param.Cradi       = 2.0;              % inner radius fo the spheroid (units)

Param.Morphology        = 0; 
Param.CellTypes         = 'Epithelial Only';
Param.FibroblastCount   = 0;

% Setting up Diffusion Coefficients/ ECM representation
Param.Dcell         = 38;                % TGFB Diffusion coefficient for epithelial cells (pixel2/hr)
Param.DThresh       = [ 0.0024, 0.004, 0.0084, 0.01 ]; % EMT-dependent changes in ECM remodeling

Param.DtgfbInhibitionMethod = 1; % 1 = Local; 2 = Global
Param.EMTThresh     = 1; % 0 - constant; 1 - cells can remodel

% Variables defining cell activity
Param.NcadMax          = 3.1515;
Param.ProlifThresh     = 24;   % Time (hr) until proliferation
% Param.ApopThresh     = 5;
Param.PThresh          = [ 0.002, 0.004, 0.008, 0.01 ]; % Probability thresholds governing proliferation, migration, apoptosis


% Defining Intracellular TGFB ODE Signaling Model
% Adapted from Tian et al. (2013). DOI: 10.1016/j.bpj.2013.07.011

Param.conc0       = [0.01, 0.01, 0.38, 0.03, 0.01, 0.25, 3.20, 0];
                 % [snail, SNAIL, miRNA34, zeb1, Zeb1, miRNA200, ECAD, NCAD]

Param.Mconc0      = [0.324, 2.779, 0.035, 0.290, 2.79, 0.015, 0.071, Param.NcadMax];                 
                 % [snail, SNAIL, miRNA34, zeb1, Zeb1, miRNA200, ECAD, NCAD]

Param.HillCoeff_J = [0.06, 1.6, 0.08, 0.15, 0.36, 3.5, 0.06, 5, 0.2, 0.2, 0.5, 0.2, 0.5]; 
    % (uM) [JT, Js, JS, J13, J23, Jz, JZ, J12, J22, J1e, J2e, J1n, J2n] 

Param.HillCoeff_n = [2, 2, 2, 2, 2]; % [nnt, nns, nnz, nnr2, nnr3]

Param.Krates = [0.0659, 1.425, 0.086, 0.0062, 0.03, 0.09, 17, 1.66, ...
    0.0012, 0.012, 0.035, 0.003, 0.06, 0.09, 17, 1.66,...
    0.0002, 0.012, 0.035, 1, 0.6, 0.5, 1, 0.6, 0.5]; 

% Reaction Rates: [k0T, kT, kdT, k0s, ks, kds, kS, kdS, 
%       k03, k3, kd3, k0z, kz, kdz, kZ, kdZ, 
%       k02, k2, kd2, ke1, ke2, kde, kn1, kn2, kdn]


end

