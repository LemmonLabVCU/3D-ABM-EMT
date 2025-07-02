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
% Setting up the meshgrid/ 3D plot and epithelial cell arrangement:
Param = RunParam;          
n = Param.n;

L           = n*Param.Csize;    % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
dx          = round(L/(n-1));   % distance intervals of the grid (um)
dy = dx; dz = dx;               % distance intervals of the grid (um)

% Setting up Plotting Options/ Parameters
xslice      = L/2;              % x plane that will be shown
yslice      = [];               % y plane that will be shown
zslice      = L/2;              % z plane that will be shown

x           = linspace(0, L, n); % converts pixels to um
y           = linspace(0, L, n); % converts pixels to um 
z           = linspace(0, L, n); % converts pixels to um
[X, Y, Z]   = meshgrid(x, y, z); % Setting up the 3D grid

% Do you want to save a .gif image of the run? 
getVideo = 0;   % Record Video = 1   ; No Video = 0
videoName = 'Insert_Video_Name_Here.gif'; % define name of output video

% Setting up the Time Scale
Tfinal      = 10*24;            % total time of the run (hr)
dt          = 1;                % time step (hr)
tskip       = 12;               % define the time intervals shown in the .gif
t           = 1;                % initialize time counter
Tint        = 0*24;             % Time interval of pulsed signaling (hr)

% Setting up Exogeneous TGFB1 variable
TGFB        = 4;                % exogeneous TGFB0 (uM)
Ctgfb       = TGFB*ones(n,n,n); % TGFB concentration gradient
TGFBMax     = TGFB + 2;
TGFBstimuli = 'Constant'; 

% Setting up Diffusion Coefficients/ ECM representation
Dcyt    = 38;     % TGFB diffusion coefficient (1/hr)
Dtgfb   = Dcyt*ones(n,n,n);  % Diffusion values matrix (nxnxn)

CellState   = struct();
CellState   = DefineCells(CellState); 
    % Set up structure (CellState) that holds cell-specific variable
 
MorphType = 0;  % 0 - spheroid, 1 - vertical straight tube, 2 = horizontal straight tube, 3 - curved tube
CultureModel = 'Epithelial Only';
FibroblastCount = 0;

[CellState, Ctgfb, Dtgfb, Cstate] = CellMorphology(CellState, Param, Ctgfb, Dtgfb, MorphType, CultureModel, FibroblastCount);

% Plotting Initial Figure Output
figure(1); set(gcf, 'Position', [293,116,1317,842]);
tcl = tiledlayout(2,3);
title(tcl, replace(videoName, '_', ' '), 'FontSize', 20);

plotFigs_SingleRun(CellState, Param, Cstate, Ctgfb, Dtgfb, Tfinal, 1, TGFBMax, Dcyt);

if getVideo == 1 % To record a .gif of the run
    F1(1) = getframe(gcf); drawnow;
    im = frame2im(F1(1));
    [imind{1}, cm{1}] = rgb2ind(im, 256); % use color
    imwrite(imind, cm, videoName, 'gif', 'Loopcount', inf); % save image to the loop
end

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
                Param, Cstate, i, j, k, v, Ctgfb, t, Cmoved, Dtgfb) ;

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
    CellState.AvgR34(t)   = mean(CellState.conc(eInd,3));

    CellState.AvgTGFB(t)  = mean(CellState.Ctgfb(eInd));
    CellState.AvgDtgfb(t)  = mean(CellState.Dcell(eInd));

    % updating 3D positions of cells
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
   

    % Plotting Data:
    if rem(t,tskip) == 0
        figure(1); tiledlayout(2,3); 

        plotFigs_SingleRun(CellState, Param, Cstate, Ctgfb, Dtgfb, Tfinal, t, TGFBMax, Dcyt);

        if getVideo == 1
            F1(t/tskip) = getframe(gcf);
        end
    end  


end

%% Collecting Data from Run
    
% Setup Video Collection
if getVideo == 1
    for nn = 1:length(F1)
        frame = F1(nn);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, videoName, Filename, 'gif', 'WriteMode', 'append');

    end
end

% Changes in cell movement
for i = 1:length(CellState.Position)
    P1 = CellState.Position(i,1:3); 
    P2 = CellState.Position(i,4:6);   
    DiffMove = sqrt((P2(1) - P1(1))^2 + (P2(2) - P1(2))^2 + (P2(3) - P1(3))^2);
    CellState.movement(i) = DiffMove; 
end

% Collecting Population Dynamics
Cell_total = sum(CellState.Pop, 2); 
CellRatio = CellState.Pop ./ Cell_total;

Outputs = [];
TimeInt = [1:tskip:Tfinal, Tfinal];
dd = 1;

for i = 1:length(TimeInt)
    
    cc = TimeInt(i);

    Outputs(:, i) = [Cell_total(cc)'; CellState.Pop(cc,1:3)'; CellRatio(cc,1:3)'; ...
        CellState.AvgEcad(cc); CellState.AvgSnail(cc); CellState.AvgNcad(cc); ...
        CellState.AvgZeb(cc); CellState.AvgR200(cc); CellState.AvgR34(cc); CellState.AvgTGFB(cc); ...
        CellState.AvgDtgfb(cc)];
    
end

    