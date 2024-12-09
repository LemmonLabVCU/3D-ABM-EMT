function [] = plotFigs_SingleRun(CellState, Param, Cstate, Ctgfb, Dtgfb, Tfinal, t, TGFBMax, DMax)

% Setting up the meshgrid and epithelial cell arrangement:
Csize       = Param.Csize;           % size of the cell (um)
Crado       = Param.Crado;           % outer radius fo the spheroid (pixels)
Cradi       = Param.Cradi;           % inner radius fo the spheroid (pixels)
n           = Param.n;               % grid size (units)

L           = n*Csize;          % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
dx          = round(L/(n-1));   % distance intervals of the grid (um)
dy = dx; dz = dx;               % distance intervals of the grid (um

% Setting up Plotting Options/ Parameters
xslice      = L/2;              % x plane that will be shown
yslice      = [];               % y plane that will be shown
zslice      = L/2;               % z plane that will be shown

x           = linspace(0, L, n); % converts pixels to um
y           = linspace(0, L, n); % converts pixels to um 
z           = linspace(0, L, n); % converts pixels to um
[X, Y, Z]   = meshgrid(x, y, z); % Setting up the 3D grid

CellCnt     = length(CellState.state);

% _______________________________________________________________________

% Displays cell movement + EMT
ax1 = nexttile; 
Cellpnts = find(Cstate(:,:,:) == 1);
[C2, C1, C3] = ind2sub([n n n], Cellpnts); 
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[0 1 0]); 
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 2);
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 1 0]);
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 3);
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 0 0]);
k.SizeData = 100;
hold on;        
Cellpnts = find(Cstate(:,:,:) == 4);
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
k.SizeData = 50;        
hold on;
subtitle(sprintf('Cell State: T = %.0f hours', t), 'FontSize', 18);
xlim([0 L]); ylim([0 L]); zlim([0 L]); 
xlabel('x'); ylabel('y'); zlabel('z');
hold off; ax1.FontSize = 16; 


% 2D sideview
ax2 = nexttile; 
Cellpnts = find(Cstate(:,:,:) == 1); % Epithelial Cells
[C2, C1, C3] = ind2sub([n n n], Cellpnts); 
kind = find(C2 >= h-2 & C2 <= h+2);
k = scatter3(z(C1(kind)), y(C2(kind)), x(C3(kind)), 'MarkerFaceColor',[0 1 0]); view(-90, 0)
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 2); % Partial Cells
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
kind = find(C1 >= h & C1 <= h);
k = scatter3(x(C1(kind)), y(C2(kind)), z(C3(kind)), 'MarkerFaceColor',[1 1 0]);
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 3); % Mesenchymal Cells
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
kind = find(C2 >= h-2 & C2 <= h+2);
k = scatter3(x(C1(kind)), y(C2(kind)), z(C3(kind)), 'MarkerFaceColor',[1 0 0]);
k.SizeData = 100;
hold on;               
Cellpnts = find(Cstate(:,:,:) == 4); % Fibroblasts
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
kind = find(C2 >= h-2 & C2 <= h+2);
k = scatter3(x(C1(kind)), y(C2(kind)), z(C3(kind)), 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
k.SizeData = 50;  

subtitle('Cell State Side View', 'FontSize', 18); 
xlim([0 L]); ylim([0 L]); zlim([0 L]); grid off;
xlabel('x'); ylabel('y'); zlabel('z');
legend('Epithelial', 'Partial', 'Mesenchymal', 'FontSize', 16, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
ax2.FontSize = 16; 


% Change in Total Popuation Over time
ax3 = nexttile;
plot(1:t, CellState.Pop(1:t, 1), 'g*-'); hold on; % Epithelial
plot(1:t, CellState.Pop(1:t, 2), '*-', 'MarkerEdgeColor',[0.90 0.60 0.10]); % Partial
plot(1:t, CellState.Pop(1:t, 3), 'r*-'); % Mesenchymal
subtitle('EMT State Ratio', 'FontSize', 18);
xlim([0, Tfinal]); ylim([0 CellCnt]);
xlabel('Time (hr)'); ylabel('EMT State Ratio');
ax3.FontSize = 16;

% % Displays TGFB Gradient
ax4 = nexttile;
slice(X, Y, Z, Ctgfb, xslice, yslice, zslice);  % view(90,0); 
subtitle('TGFB Gradient', 'FontSize', 18);
colormap (ax4, jet); colorbar; caxis([0 TGFBMax]);
xlim([0 L]); ylim([0 L]); zlim([0 L]);
ax4.FontSize = 16;

% Displays TGFB Diffusion Coefficient values in the system
ax5 = nexttile;
slice(X, Y, Z, Dtgfb, xslice, yslice, zslice);  % view(90,0); 
subtitle('Local Dtgfb Values', 'FontSize', 18);
colormap(ax5, jet); colorbar; caxis([0 DMax]);
xlim([0 L]); ylim([0 L]); zlim([0 L]);
ax5.FontSize = 16;

% Displays changes in EMT Markers over time
ax6 = nexttile;
plot(1:t, CellState.AvgEcad, 'g-', 'LineWidth',2); hold on; 
plot(1:t, CellState.AvgNcad, 'r-','LineWidth',2);
plot(1:t, CellState.AvgSnail, 'm-', 'LineWidth',2);
plot(1:t, CellState.AvgZeb, 'b-','LineWidth',2); 
plot(1:t, CellState.AvgR34, 'k-','LineWidth',2); 
plot(1:t, CellState.AvgR200, 'k--','LineWidth',2);
subtitle('EMT Markers', 'FontSize', 18);
xlabel('Time (Days)', 'FontSize', 16); ylabel('Concentration', 'FontSize', 16);
legend('Ecad', 'Ncad', 'Snail1', 'Zeb1', 'miR34', 'miR200', 'Location', 'Eastoutside');
xlim([0 Tfinal]); ylim([0 4]); xlabel('Time (hr)'); ylabel('Concentration');
hold off; 
ax6.FontSize = 16;

end
