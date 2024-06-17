function [] = plotFigs2(CellState, Param, Cstate, Ctgfb, Dtgfb, Tfinal, t, TGFBMax)

% Setting up the meshgrid and epithelial cell arrangement:
Csize       = Param.Csize;               % size of the cell (um)
Crado       = Param.Crado;                % outer radius fo the spheroid (pixels)
Cradi       = Param.Cradi;              % inner radius fo the spheroid (pixels)
n           = Param.n;               % grid size (units)

L           = n*Csize;          % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
dx          = round(L/(n-1));   % distance intervals of the grid (um)
dy = dx; dz = dx;               % distance intervals of the grid (um

% Setting up Plotting Options/ Parameters
xslice      = L/2;              % x plane that will be shown
yslice      = [];               % y plane that will be shown
zslice      = [];              % z plane that will be shown

x           = linspace(0, L, n); % converts pixels to um
y           = linspace(0, L, n); % converts pixels to um 
z           = linspace(0, L, n); % converts pixels to um
[X, Y, Z]   = meshgrid(x, y, z); % Setting up the 3D grid

% _______________________________________________________________________

% Displays cell movement + EMT
ax1 = nexttile;
Cellpnts = find(Cstate(:,:,:) == 1);
[C1, C2, C3] = ind2sub([n n n], Cellpnts); 
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[0 1 0]); 
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 2);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 1 0]);
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 3);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 0 0]);
k.SizeData = 100;
hold on;        
Cellpnts = find(Cstate(:,:,:) == 4);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
k.SizeData = 50;        
hold on;
subtitle(sprintf('Cell State: T = %.0f hours', t), 'FontSize', 18);
xlim([0 L]); ylim([0 L]); zlim([0 L]);
xlabel('x'); ylabel('y'); zlabel('z');
legend('Epithelial', 'Partial', 'Mesenchymal', 'Fibroblast', 'FontSize', 16);

hold off;

% 2D sideview
ax2 = nexttile; 
Cellpnts = find(Cstate(:,:,:) == 1);
[C1, C2, C3] = ind2sub([n n n], Cellpnts); 
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[0 1 0]); view(90, 0)
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 2);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 1 0]);
k.SizeData = 100;
hold on;
Cellpnts = find(Cstate(:,:,:) == 3);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 0 0]);
k.SizeData = 100;
hold on;        
Cellpnts = find(Cstate(:,:,:) == 4);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(x(C1), y(C2), z(C3), 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
k.SizeData = 50;        
subtitle(sprintf('Cell State Side View: T = %.0f hours', t), 'FontSize', 18); 
xlim([0 L]); ylim([0 L]); zlim([0 L]); grid off;
ax2.FontSize = 16; 

% Displays TGFB Gradient
ax3 = nexttile;
slice(X, Y, Z, Ctgfb, xslice, yslice, zslice); 
subtitle(sprintf('TGFB Gradient: T = %.0f hours', t), 'FontSize', 18);
colormap (ax3, jet); colorbar; caxis([0 TGFBMax]);

% Displays TGFB Diffusion Coefficient values in the system
ax4 = nexttile;
slice(X, Y, Z, Dtgfb, xslice, yslice, zslice); 
% contourslice(X,Y,Z, Dtgfb, [], [], 0:L); view(-11,14)
subtitle(sprintf('Dtgfb Gradient: T = %.0f hr', t), 'FontSize', 18);
colormap(ax4, jet); colorbar; caxis([0 38]);

end