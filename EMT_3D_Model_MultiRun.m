% Modeling TGF-B1 induced Epithelial-Mesenchymal Transition in a
% multicellular system using a hybrid cell automaton model. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMT_3D_Model_master.m

% This script simulates the progression of Epithelial-Mesenchymal 
% Transition in response to exogeneous TGF-B1 added to the 3D system.
% The ODE model of TGF-B1 signaling is adapted from the model described by
% Tian et al (2013), where we expand and integrate it to include spatial
% and temporal changes in each cell defined in the system.

% Created by Kristin Kim: kimkp@vcu.edu
% Last Modified by: Kristin Kim, May 24, 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean Up Section
clear variables;
clc;

%% Set up run options
tic

% Define the name of the xcel sheet the data will save to
Xfilename = 'Exmpl2.xlsx';
% File Name:[DATE_CELL type started with_TFinal_Ncount.xlsx]

DataOutput = [];

% Defining Parameters to change during the run
DiffCoeff = [0, 4, 9, 19, 38]; % Changing the Dcyt value

TGFBConc = [0, 2, 4]; % Changing TGFB concentration added

for Ndtgfb = 1:length(DiffCoeff)

for Ntgfb = 1:length(TGFBConc)
    
    Ncount = 1;     %%%%%% How many repeated runs do you want to have?  
    
    disp(['Dcyt = ', num2str(DiffCoeff(Ndtgfb)), ' TGFBConc = ', num2str(TGFBConc(Ntgfb))]);
    
for Nruns = 1:Ncount

% Setting up the meshgrid/ 3D plot and epithelial cell arrangement:
Param = RunParam;          
n = Param.n;

L           = n*Param.Csize;    % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
dx          = round(L/(n-1));   % distance intervals of the grid (um)
dy = dx; dz = dx;               % distance intervals of the grid (um)

x           = linspace(0, L, n); % converts pixels to um
y           = linspace(0, L, n); % converts pixels to um 
z           = linspace(0, L, n); % converts pixels to um
[X, Y, Z]   = meshgrid(x, y, z); % Setting up the 3D grid

% Setting up the Time Scale
Tfinal      = 10*24;            % total time of the run (hr)
dt          = 1;                % time step (hr)
tskip       = 24;               % define the time intervals shown in the .gif
Tint        = 0*24;             % Time interval of pulsed signaling (hr)

% Setting up Exogeneous TGFB1 variable
TGFB        = TGFBConc(Ntgfb);  % exogeneous TGFB0 (uM)
Ctgfb       = TGFB*ones(n,n,n); % TGFB concentration gradient
TGFBMax     = TGFB + 2;
TGFBstimuli = 'Constant'; 

% Setting up Diffusion Coefficients/ ECM representation
Dcyt    = DiffCoeff(Ndtgfb); % TGFB diffusion coefficient (1/hr)
Dtgfb   = Dcyt*ones(n,n,n);  % Diffusion values matrix (nxnxn)

CellState   = struct();
CellState   = DefineCells(CellState); 
    % Set up structure (CellState) that holds cell-specific variable
 
MorphType = 0;  % 0 - spheroid, 1 - vertical straight tube, 2 = horizontal straight tube, 3 - curved tube
CultureModel = 'Epithelial Only';
FibroblastCount = 0;

[CellState, Ctgfb, Dtgfb, Cstate] = CellMorphology(CellState, Param, Ctgfb, Dtgfb, MorphType, CultureModel, FibroblastCount);

C0tgfb      = Ctgfb;

%% Setting up the time loop

for t = 2:dt:Tfinal

   switch TGFBstimuli

        case 'Removed'
           if t > TimeInt
               for i = 1:n
                    for j = 1:n
                        for k = 1:n
                            if Cstate(i,j,k) == 0
                                Ctgfb(i,j,k) = 0;
                            end
                        end
                    end
                end
            end

        case 'Pulsed'
           if rem(t, TimeInt*2) == 0
               for i = 1:n
                    for j = 1:n
                        for k = 1:n
                            if Cstate(i,j,k) == 0
                                Ctgfb(i,j,k) = Ctgfb(i,j,k) + TGFB;
                            end
                        end
                    end
               end
           elseif rem(t,TimeInt) == 0 && rem(t,TimeInt*2)~= 0
               for i = 1:n
                    for j = 1:n
                        for k = 1:n
                            if Cstate(i,j,k) == 0
                                Ctgfb(i,j,k) = 0;
                            end
                        end
                    end
               end
           end

        case 'Increasing Interval'
           if rem(t, TimeInt) == 0
               for i = 1:n
                    for j = 1:n
                        for k = 1:n
                            if Cstate(i,j,k) == 0
                                Ctgfb(i,j,k) = Ctgfb(i,j,k) + TGFB;
                            end
                        end
                    end
               end
           end

        case 'Constant'
    
   end

    Cmoved = zeros(n,n,n); % monitor movements of cells + update in system
    [CellState, Cstate, Ctgfb, Dtgfb] = EMTStateChange(CellState, Param, Cstate, Ctgfb, Dtgfb, t, Tfinal, dt);

    % Update Cell Behaviors
    flag = 0;
    v = 1; numCells = length(CellState.state);
    
    while flag == 0

        i = CellState.Position(v,4); 
        j = CellState.Position(v,5); 
        k = CellState.Position(v,6);
        
        if Cmoved(i,j,k) == 0 && CellState.state(v) > 0

            [CellState, Cstate, Cmoved, Ctgfb, Dtgfb] = CellBehavior(CellState,...
                Param, Cstate, i, j, k, v, Ctgfb, t, Cmoved, Dtgfb);

        end
        
        v = v + 1; numCells = length(CellState.state);
        
        if v >= numCells
            flag = 1;
        end
    end
    
    % Updating TGFB Diffusion 
    Ctgfb = TGFBDiffusion(Param, Ctgfb, TGFB, Dtgfb);
    
    % Updating population cell concentrations/ values
    eInd = find(CellState.state < 4); % Only use cells that made up the initial spheroid

    CellState.AvgEcad(t)  = mean(CellState.conc(eInd, 7));
    CellState.AvgNcad(t)  = mean(CellState.Ncad(eInd, 2));
    CellState.AvgSnail(t) = mean(CellState.conc(eInd, 2));
    CellState.AvgZeb(t)   = mean(CellState.conc(eInd, 5));
    CellState.AvgR200(t)  = mean(CellState.conc(eInd, 6)); 
    CellState.AvgR34(t)   = mean(CellState.conc(eInd, 3)); 
    
    CellState.AvgTGFB(t)  = mean(CellState.Ctgfb(eInd));
    CellState.AvgDtgfb(t)  = mean(CellState.Dcell(eInd));


    fInd = find(CellState.state == 4);
    Cstate2 = Cstate;

    for nn = 1:length(fInd)
        [posx] =  CellState.Position(fInd(nn), 4:6);
        Cstate2(posx(1), posx(2), posx(3)) = 0;
    end

    % Finding CS area and diameters of spheroid
    Imgdim = sum(Cstate); Imgdim = squeeze(Imgdim); Imgdim = imfill(Imgdim);
    Imgdim2 = permute(Imgdim, [1 3 2]); 
    Imgdim2 = reshape(Imgdim,[],size(Imgdim,2),1);
    Imgdim2 = imbinarize(Imgdim2);
    tempStat = regionprops(Imgdim2, 'Area', 'MajorAxisLength','MinorAxisLength');
    tempStat = struct2table(tempStat);
    [~, maxInd] = max(tempStat.Area);
    
    CellState.Carea(t) = tempStat.Area(maxInd);
    CellState.MajorD(t) = tempStat.MajorAxisLength(maxInd);
    CellState.MinorD(t) = tempStat.MinorAxisLength(maxInd); 
   
end

%% Collecting Data from Run

incr = Nruns + Ncount*(Ntgfb - 1); % for each run, increase index value by # runs per condition to prevent data overwriting
DataOutput = SaveExcelResults(CellState, Param, Cstate, Ctgfb, Dtgfb, incr, DataOutput);

end % Finish Nruns loop

end % Finish Ntgfb loop

% Save outputs from DataOutput in an excel sheet
sheetData = ['Dtgfb = ', num2str(DiffCoeff(Ndtgfb))]; % save values in an xcel sheet
writetable(struct2table(DataOutput), Xfilename, 'Sheet', sheetData, 'Range', 'A1');

end % Finish Ndtgfb loop

totRunTime = toc;

%%
sheetParam = 'Run Parameters';

runParam.ProcessingTime = totRunTime;
runParam.Ncount = Ncount; 
runParam.Tfinal = Tfinal;

runParam.CellSize = Param.Csize;
runParam.InnerRadii = Param.Cradi;
runParam.OuterRadii = Param.Crado;
runParam.NumBlocks = Param.n;

runParam.Dcell = DiffCoeff;
runParam.TGFB  = TGFBConc; 
runParam.Morph = MorphType;
runParam.CultureModel = CultureModel;
runParam.TGFBDynamics = TGFBstimuli;

runParam.ProlifThreshold = Param.ProlifThresh;
runParam.ProbThresholds = Param.PThresh;
runParam.Conc0 = Param.conc0;

ParamT = struct2cell(runParam);

ParamNames = {'Processing Speed', 'Repeat Runs', 'Tfinal', 'Cell Size',...
    'Inner Radii', 'Outer Radii', '# voxels/edge', 'Dcell', 'TGFB', ...
    'Morphology', 'Culture Model', 'TGFB Signaling Dynamics', ...
    'Proliferation Threshold', 'Prob Thresholds', 'ODE Conc(0)'};

writecell([ParamNames', ParamT], Xfilename, 'Sheet', sheetParam, 'Range', 'A1');
