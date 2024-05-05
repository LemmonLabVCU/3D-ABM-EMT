function plotFigs(CellState, Param, Cstate, Ctgfb, Dtgfb, TGFB, Tfinal, t, UIAxes, UIAxes2, UIAxes3, UIAxes4, UIAxes5, UIAxes6)

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
zslice      = [];              % z plane that will be shown

x           = linspace(0, L, n); % converts pixels to um
y           = linspace(0, L, n); % converts pixels to um 
z           = linspace(0, L, n); % converts pixels to um
[X, Y, Z]   = meshgrid(x, y, z); % Setting up the 3D grid

% _______________________________________________________________________
    
% Display cell movement + EMT
Cellpnts = find(Cstate(:,:,:) == 1);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), 'MarkerFaceColor',[0 1 0]);
k.SizeData = 100;
hold(UIAxes, 'on');

Cellpnts = find(Cstate(:,:,:) == 2);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 1 0]);
k.SizeData = 100;
hold(UIAxes, 'on');

Cellpnts = find(Cstate(:,:,:) == 3);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), 'MarkerFaceColor',[1 0 0]);
k.SizeData = 100;
hold(UIAxes, 'on');

xlim(UIAxes, [0 L]); ylim(UIAxes, [0 L]); zlim(UIAxes, [0 L]);

Cellpnts = find(Cstate(:,:,:) == 4);
[C1, C2, C3] = ind2sub([n n n], Cellpnts);
k = scatter3(UIAxes, x(C1), y(C2), z(C3), 'MarkerFaceColor',[0.4940 0.1840 0.5560]);
k.SizeData = 50; % legend(UIAxes, 'Epithelial', 'Partial', 'Mesenchymal', 'Fibroblast');
hold(UIAxes, 'off');

% Display hanges in EMT cell population
plot(UIAxes2, 1:t, CellState.Pop(1:t, 1), 'g*-'); % Epithelial
hold (UIAxes2, 'on'); 
plot(UIAxes2,  1:t, CellState.Pop(1:t, 2), '*-', 'MarkerEdgeColor',[0.90 0.60 0.10]); % Partial
hold (UIAxes2, 'on'); 
plot(UIAxes2,  1:t, CellState.Pop(1:t, 3), 'r*-'); % Mesenchymal

legend(UIAxes2, 'Epithelial', 'Partial', 'Mesenchymal'); % 0.449435941221106,0.879648578811369,0.130242825607064,0.085271317829457
xlim(UIAxes2, [0 Tfinal]); ylim(UIAxes2, [0 250]);
hold (UIAxes2, 'off'); 

% Display changes in Average EMT concentrations over time
plot(UIAxes3, 1:t, CellState.AvgEcad, 'g-', 1:t, CellState.AvgNcad, 'r-', 1:t, CellState.AvgSnail, 'b-');
xlabel(UIAxes3, 'Time (hr)'); ylabel(UIAxes3, 'Concentration'); xlim(UIAxes3, [0 Tfinal]);
legend(UIAxes3, 'Ecad', 'Ncad', 'Snail'); % 0.88021366468985,0.886108532060332,0.084988962472406,0.085271317829457

% Display Changes in Ecadherin v. Ncadherin expression
Cellpnts = find(CellState.state(:,:,:) == 1);
plot(UIAxes4, CellState.Ecad(Cellpnts), CellState.Ncad(Cellpnts, 2), 'g*');
hold (UIAxes4, 'on'); 

Cellpnts = find(CellState.state(:,:,:) == 2);
plot(UIAxes4, CellState.Ecad(Cellpnts), CellState.Ncad(Cellpnts, 2), '*', 'MarkerEdgeColor',[0.90 0.60 0.10]);

Cellpnts = find(CellState.state(:,:,:) == 3);
plot(UIAxes4, CellState.Ecad(Cellpnts), CellState.Ncad(Cellpnts, 2), 'r+');      

yline(UIAxes4, 1.5); 
yline(UIAxes4, 2.4); 
xlabel(UIAxes4, 'ECAD'); ylabel(UIAxes4, 'NCAD');
xlim(UIAxes4, [0 3.5]); ylim(UIAxes4, [0 3.5]);
hold (UIAxes4, 'off'); 

% Display TGFB gradient
slice(UIAxes5, X, Y, Z, Ctgfb, xslice, yslice, zslice); view(UIAxes5, 90, 0);
colormap(UIAxes5, jet); colorbar(UIAxes5); clim(UIAxes5, [0 TGFBMax]);

% Display TGFB Diffusion Coefficient Values (Dtgfb)
slice(UIAxes6, X, Y, Z, Dtgfb, xslice, yslice, zslice); view(UIAxes6, 90, 0);
colormap(UIAxes6, jet); colorbar(UIAxes6); clim(UIAxes6, [0 38]);

end