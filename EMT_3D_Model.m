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
% Last Modified by: Kristin Kim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean Up Section
clear variables;
clc;


%% Set up run options
% Setting up the meshgrid/ 3D plot and epithelial cell arrangement:
Param = RunParam;          
n = Param.n;

L           = n*Param.Csize;          % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
dx          = round(L/(n-1));   % distance intervals of the grid (um)
dy = dx; dz = dx;               % distance intervals of the grid (um)

% Setting up Plotting Options/ Parameters
xslice      = L/2;              % x plane that will be shown
yslice      = [];               % y plane that will be shown
zslice      = [];              % z plane that will be shown

x           = linspace(0, L, n); % converts pixels to um
y           = linspace(0, L, n); % converts pixels to um 
z           = linspace(0, L, n); % converts pixels to um
[X, Y, Z]   = meshgrid(x, y, z); % Setting up the 3D grid

% Do you want to save a .gif image of the run? 
getVideo = 0;   % Record Video = 1   ; No Video = 0
videoName = 'Insert_Name_Here.gif'; % define name of output video

% Setting up the Time Scale
Tfinal      = 15*24;            % total time of the run (hr)
dt          = 1;                % time step (hr)
tskip       = 24;               % define the time intervals shown in the .gif

% Setting up Exogeneous TGFB1 variable
TGFB        = 0;                % exogeneous TGFB0 (uM)
Ctgfb       = TGFB*ones(n,n,n); % TGFB concentration matrix
TGFBMax     = TGFB + 2;         % Max Concentration shown on colorbar when plotting
TGFBstimuli = 0;                % TGFB Signaling Dynamics: 0: Constant; 1: Removed; 2: Pulsed; 3: Added
TimeInt     = 0*24;             % Time interval of pulsed/ removal of TGFB1

% Setting up Diffusion Coefficients/ ECM representation
Dcyt    = 38;                % TGFB diffusion coefficient (1/hr)
Dtgfb   = Dcyt*ones(n,n,n);  % Diffusion values matrix (nxnxn)

CellState   = struct();
CellState   = DefineCells(CellState); 
    % Set up structure (CellState) that holds cell-specific variable
 
MorphType = 0;  % Cell Organization: 0 - spheroid, 1 - vertical straight tube, 2 = horizontal straight tube, 3 - curved tube
CultureModel = 'Epithelial Only';
FibroblastCount = 0; % for co-culture models

[CellState, Ctgfb, Dtgfb, Cstate] = CellMorphology(CellState, Param, Ctgfb, Dtgfb, MorphType, CultureModel, FibroblastCount);

% Plotting Initial Figure Output
figure(1); set(gcf, 'Position', [1142,42,778,954]); % trilinear : set(gcf, 'Position', [739,42,1178,954]);
tiledlayout(2,2);

plotFigs2(CellState, Param, Cstate, Ctgfb, Dtgfb, Tfinal, 1, TGFBMax);

if getVideo == 1 % To record a .gif of the run
    F1(1) = getframe(gcf); drawnow;
    im = frame2im(F1(1));
    [imind{1}, cm{1}] = rgb2ind(im, 256); % use color
    imwrite(imind, cm, videoName, 'gif', 'Loopcount', inf); % save image to the loop
end

%% Setting up the time loop
for t = 2:dt:Tfinal

    switch TGFBstimuli
        case 1 % Removed
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

        case 2 % Pulsed
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

        case 3 % Increased Addition
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

        case 0 % Constant/ Single Bolus

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
            
    % Plotting Data:
    if rem(t,tskip) == 0
        figure(1); tiledlayout(2,2); 

        plotFigs2(CellState, Param, Cstate, Ctgfb, Dtgfb, Tfinal, t, TGFBMax);

        if getVideo == 1
            F1(t/tskip) = getframe(gcf);
        end
    end  
    
    % Updating population cell concentrations/ values
    eInd = find(CellState.state < 4); % Only use cells that made up the initial spheroid
    CellState.AvgEcad(t)  = mean(CellState.conc(eInd, 7));
    CellState.AvgNcad(t)  = mean(CellState.Ncad(eInd, 2));

    CellState.AvgSnail(t) = mean(CellState.conc(eInd, 2));
    CellState.AvgZeb(t)   = mean(CellState.conc(eInd, 5));

    CellState.AvgR200(t)       = mean(CellState.conc(eInd, 6)); 
    CellState.AvgR34(t)        = mean(CellState.conc(eInd, 3));
        
    fInd = find(CellState.state == 4);
    Cstate2 = Cstate;
    for nn = 1:length(fInd)
        [posx] =  CellState.Position(fInd(nn), 4:6);
        Cstate2(posx(1), posx(2), posx(3)) = 0;
    end
    Imgdim = sum(Cstate2); Imgdim = squeeze(Imgdim); Imgdim = imfill(Imgdim);
    CellState.Carea(t) = length(find(Imgdim > 0));   
   
end

%% Collecting Data from Run

% Changes in cell movement
for i = 1:length(CellState.Position)
    P1 = CellState.Position(i,1:3); 
    P2 = CellState.Position(i,4:6);   
    DiffMove = sqrt((P2(1) - P1(1))^2 + (P2(2) - P1(2))^2 + (P2(3) - P1(3))^2);
    CellState.movement(i) = DiffMove; 
end
    
% Setup Video Collection
if getVideo == 1
    for nn = 1:length(F1)
        frame = F1(nn);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, videoName, Filename, 'gif', 'WriteMode', 'append');

    end
end

%% Setting up other figure outputs
figure(2); 
set(gcf, 'Position', [1100,110,469,823]); % [636,462,546,484]); % [903,42,386,954]); %set(gcf, 'Position', [833,42,485,954]);
Tplots = tiledlayout(3, 1);
Tplots.Padding = 'tight';
Tplots.TileSpacing = 'compact';

% Displays cell movement + EMT
ax1 = nexttile; 
Cellpnts = find(Cstate(:,:,:) == 1);
[C1, C2, C3] = ind2sub([n n n], Cellpnts); 
kind = find(C1 >= h-4 & C1 <= h+4);
k = scatter3(x(C1(kind)), y(C2(kind)), z(C3(kind)), 'MarkerFaceColor',[0 1 0]); view(90, 0)
k.SizeData = 100;
hold on;

Cellpnts = find(Cstate(:,:,:) == 2);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
kind = find(C1 >= h-2 & C1 <= h+2);
k = scatter3(x(C1(kind)), y(C2(kind)), z(C3(kind)), 'MarkerFaceColor',[1 1 0]);
k.SizeData = 100;
hold on;

Cellpnts = find(Cstate(:,:,:) == 3);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
kind = find(C1 >= h-2 & C1 <= h+2);
k = scatter3(x(C1(kind)), y(C2(kind)), z(C3(kind)), 'MarkerFaceColor',[1 0 0]);
k.SizeData = 100;
hold on;        
Cellpnts = find(Cstate(:,:,:) == 4);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
k.SizeData = 50;

% legend('Epithelial', 'Partial', 'Mesenchymal', 'Fibroblast');
subtitle('Cell Morphology', 'FontSize', 18);
xlim([0 L]); ylim([0 L]); zlim([0 L]); grid off;
set(ax1, 'XTickLabel', {}, 'YTickLabel', {}, 'ZTickLabel', {});
ax1.FontSize = 16; 

% Displays TGFB1 gradient
ax2 = nexttile; 
slice(X, Y, Z, Ctgfb, xslice, yslice, zslice); view(90,0); 
        xlim([0 L]); ylim([0 L]); zlim([0 L]);  
subtitle('TGFB Gradient', 'FontSize', 18); 
colormap (ax2, jet); colorbar; caxis([0 TGFBMax]); 
set(ax2, 'XTickLabel', {}, 'YTickLabel', {}, 'ZTickLabel', {});
ax2.FontSize = 16;  

% Displays TGFB Diffusion Coefficient values in the system
ax3 = nexttile;
slice(X, Y, Z, Dtgfb, xslice, yslice, zslice); view(90,0)
colormap(ax3, jet); colorbar; caxis([0 Dcyt]); 
xlim([0 L]); ylim([0 L]); zlim([0 L]);      
subtitle('Diffusion Coefficient', 'FontSize', 18);
set(ax3, 'XTickLabel', {}, 'YTickLabel', {}, 'ZTickLabel', {});
ax3.FontSize = 16;

% figure(3); % set(gcf, 'Position',  [1109,347,731,631]); 
% t = (1:Tfinal)./24;
% plot(t, CellState.AvgEcad, 'g-', 'LineWidth',2);
% title('Changes in EMT Markers', 'FontSize', 18);
% xlabel('Time (hr)', 'FontSize', 16); ylabel('Concentration', 'FontSize', 16);
% hold on; 
% plot(t, CellState.AvgNcad, 'r-','LineWidth',2);
% plot(t, CellState.AvgSnail, 'm-', 'LineWidth',2);
% plot(t, CellState.AvgZeb, 'b-','LineWidth',2); 
% plot(t, CellState.AvgR34, 'k-','LineWidth',2); 
% plot(t, CellState.AvgR200, 'k--','LineWidth',2);
% 
% legend('Ecad', 'Ncad', 'Snail', 'Zeb', 'miR34', 'miR200');
% hold off;

Cell_total = CellState.Pop(:,1) + CellState.Pop(:,2) + CellState.Pop(:,3);
CellRatio = CellState.Pop;
CellRatio = CellState.Pop ./ Cell_total;

%% 
Outputs = [];
TimeInt = [1:tskip:Tfinal, Tfinal];

for i = 1:length(TimeInt)
    
    cc = TimeInt(i);
    Outputs(:, i) = [CellState.Pop(cc,1:3)'; Cell_total(cc); CellRatio(cc,1:3)'; ...
        CellState.AvgEcad(cc); CellState.AvgSnail(cc); CellState.AvgNcad(cc); ...
        CellState.AvgZeb(cc); CellState.AvgR200(cc); CellState.AvgR34(cc); CellState.Carea(cc)];
end
