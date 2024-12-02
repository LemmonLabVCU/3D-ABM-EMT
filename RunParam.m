function [Param] = RunParam( )
% RunParam.m Initialize constant parameter values being used in the model. 

% Variables for meshgrid and spheroid
Param.n           = 32;               % grid size (units)

% Initializing Cell Characteristics
Param.Csize       = 15;               % size of the cell (um)
Param.Crado       = 4.5;              % outer radius fo the spheroid (units)
Param.Cradi       = 4.0;              % inner radius fo the spheroid (units)

% Setting up Diffusion Coefficients/ ECM representation
Param.Dcell   = 38;                % TGFB Diffusion coefficient for epithelial cells (pixel2/hr)

% Variables defining cell activity
Param.NcadMax          = 3.1515;
Param.ProlifThresh     = 24;   % Time (hr) until proliferation
% Param.ApopThresh       = 5;
Param.PThresh          = [ 0.002, 0.004, 0.008 ] ; % Probability thresholds governing proliferation, migration, apoptosis

Param.conc0       = [0.01, 0.01, 0.38, 0.03, 0.01, 0.25, 3.20, 0];
                 % [snail, SNAIL, miRNA34, zeb1, Zeb1, miRNA200, ECAD, NCAD]

Param.Mconc0      = [0.324, 2.779, 0.035, 0.290, 2.79, 0.015, 0.071, Param.NcadMax];                 
                 % [snail, SNAIL, miRNA34, zeb1, Zeb1, miRNA200, ECAD, NCAD]


end

