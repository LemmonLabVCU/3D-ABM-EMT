function plotFigs(CellState, Param, Cstate, Ctgfb, Dtgfb, TGFB, Tfinal,...
    t, UIAxes, UIAxes2, UIAxes3, UIAxes4, UIAxes5, UIAxes6)

% Setting up the meshgrid and epithelial cell arrangement:
Csize       = Param.Csize;      % size of the cell (um)
Crado       = Param.Crado;      % outer radius fo the spheroid (pixels)
Cradi       = Param.Cradi;      % inner radius fo the spheroid (pixels)
n           = Param.n;          % grid size (units)

TGFBMax     = TGFB + 2;
L           = n*Csize;          % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
dx          = round(L/(n-1));   % distance intervals of the grid (um)
dy = dx; dz = dx;               % distance intervals of the grid (um

% Setting up Plotting Options/ Parameters
xslice      = L/2;              % x plane that will be shown
yslice      = [];               % y plane that will be shown
zslice      = L/2;              % z plane that will be shown

x           = linspace(0, L, n); % converts pixels to um
y           = linspace(0, L, n); % converts pixels to um 
z           = linspace(0, L, n); % converts pixels to um
[X, Y, Z]   = meshgrid(x, y, z); % Setting up the 3D grid

CellCnt     = length(CellState.state);

% _______________________________________________________________________
    
% Display cell movement + EMT
Cellpnts = find(Cstate(:,:,:) == 1);
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), 'MarkerFaceColor',[0 1 0]);
k.SizeData = 100; hold(UIAxes, 'on');

Cellpnts = find(Cstate(:,:,:) == 2);
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 1 0]);
k.SizeData = 100; hold(UIAxes, 'on');

Cellpnts = find(Cstate(:,:,:) == 3);
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 0 0]);
k.SizeData = 100; hold(UIAxes, 'on');

Cellpnts = find(Cstate(:,:,:) == 4);
[C2, C1, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), ...
    'MarkerFaceColor',[0.4940 0.1840 0.5560]);
k.SizeData = 50; 
xlim(UIAxes, [0 L]); ylim(UIAxes, [0 L]); zlim(UIAxes, [0 L]);

hold(UIAxes, 'off');

% Display changes in EMT cell population
plot(UIAxes2, 1:t, CellState.Pop(1:t, 1), 'g*-'); % Epithelial
hold (UIAxes2, 'on'); 
plot(UIAxes2,  1:t, CellState.Pop(1:t, 2), '*-', ...
    'MarkerEdgeColor',[0.90 0.60 0.10]); % Partial 

plot(UIAxes2,  1:t, CellState.Pop(1:t, 3), 'r*-'); % Mesenchymal
legend(UIAxes2, 'Epithelial', 'Partial', 'Mesenchymal', ...
    'Location', 'southoutside');
xlim(UIAxes2, [0 Tfinal]); ylim(UIAxes2, [0 CellCnt]);
hold (UIAxes2, 'off'); 

% Display changes in Average EMT concentrations over time
plot(UIAxes3, 1:t, CellState.AvgEcad, 'g-', 1:t, CellState.AvgNcad, ...
    'r-', 1:t, CellState.AvgSnail, 'b-', 1:t, CellState.AvgR34, ...
    'k-', 1:t, CellState.AvgR200, 'k--');

xlabel(UIAxes3, 'Time (hr)'); ylabel(UIAxes3, 'Concentration');
xlim(UIAxes3, [0 Tfinal]); ylim(UIAxes3, [0 4]);
legend(UIAxes3, 'Ecad', 'Ncad', 'Snail1', 'Zeb1', 'miR34', 'miR200', ...
    'Location', 'southoutside', 'orientation', 'horizontal','NumColumns',2);

% Display Changes in Ecadherin v. Ncadherin expression
Cellpnts = find(CellState.state(:,:,:) == 1);
plot(UIAxes4, CellState.conc(Cellpnts, 7), CellState.Ncad(Cellpnts, 2), 'g*');
hold (UIAxes4, 'on'); 

Cellpnts = find(CellState.state(:,:,:) == 2);
plot(UIAxes4, CellState.conc(Cellpnts, 7), CellState.Ncad(Cellpnts, 2), ...
    '*', 'MarkerEdgeColor',[0.90 0.60 0.10]);

Cellpnts = find(CellState.state(:,:,:) == 3);
plot(UIAxes4, CellState.conc(Cellpnts, 7), CellState.Ncad(Cellpnts, 2), 'r+');      

% Define Thresholds of EMT States
yline(UIAxes4, 1.5); 
yline(UIAxes4, 2.4); 
xlabel(UIAxes4, 'ECAD'); ylabel(UIAxes4, 'NCAD');
xlim(UIAxes4, [0 3.5]); ylim(UIAxes4, [0 3.5]);
hold (UIAxes4, 'off'); 

% Display TGFB gradient
slice(UIAxes5, X, Y, Z, Ctgfb, xslice, yslice, zslice);
% view(UIAxes5, 90, 0);
colormap(UIAxes5, jet); colorbar(UIAxes5); clim(UIAxes5, [0 TGFBMax]);

% Display TGFB Diffusion Coefficient Values (Dtgfb)
slice(UIAxes6, X, Y, Z, Dtgfb, xslice, yslice, zslice); 
%view(UIAxes6, 90, 0);
colormap(UIAxes6, jet); colorbar(UIAxes6); clim(UIAxes6, [0 38]);

end
