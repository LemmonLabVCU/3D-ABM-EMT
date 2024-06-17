

function [CellState, Ctgfb, Dtgfb, Cstate] = CellMorphology(CellState, Param, Ctgfb, Dtgfb, MorphType, CultureModel, FibroblastCount)
%%
% CellMorphology.m: Modifies the starting cell structure.
n = Param.n;
h = round(n/2); 
cnt = 1;

EConc = Param.conc0;
MConc = Param.Mconc0;

switch MorphType
    case 0
        % Initializing epithelial cells as a sphere
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    if ((i-h)^2 + (j-h)^2 + (k-h)^2 <= Param.Crado^2) && ((i-h)^2 + (j-h)^2 + (k-h)^2 >= Param.Cradi^2)
                       CellState.Position(cnt, :) = [i j k];   % add cells
                       CellState.Divide(cnt)     = 0;
                       CellState.conc(cnt, :)    = EConc; % initial ODE conditions
                       CellState.state(cnt)      = 1;

                       CellState.TimesMoved(cnt) = 0;

                       CellState.p(cnt)      = CellState.conc(cnt, 8) / Param.NcadMax; % probability threshold dependent on Ncad
                       Dtgfb(i,j,k)          = Param.Dcell*(1-CellState.p(cnt)); % defining diffusion coefficient for each cell
                                             
                       CellState.Ctgfb(cnt)         = Ctgfb(i,j,k); % TGFB concentration gradient
                       
                       CellState.Pop(1) =  CellState.Pop(1) + 1; % counts total cells in the system
                       cnt       = cnt + 1;       % increase cell count               
                    end
                end
            end
        end
        
        
        
    case 1
        % Initializing epithelial cells as a straight hollow cylinder -
        % VERTICAL
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    if ((i-h)^2 + (j-h)^2 <= Param.Crado^2) && ((i-h)^2 + (j-h)^2 >= Param.Cradi^2) % ((i-h)^2 + (j-h)^2 + (k-h)^2 <= Crado^2) && ((i-h)^2 + (j-h)^2 + (k-h)^2 >= Cradi^2)
                       CellState.Position(cnt, :) = [i j k];   % add cells
                       CellState.Divide(cnt)     = 0;
                       CellState.conc(cnt, :)    = EConc; % initial ODE conditions
                                % [snail, SNAIL, miRNA34, zeb1, Zeb1, miRNA200, ECAD, NCAD]
                       CellState.state(cnt)      = 1;

                       CellState.TimesMoved(cnt) = 0;

                       CellState.p(cnt)      = CellState.conc(cnt, 8) / Param.NcadMax; % probability threshold dependent on Ncad
                       Dtgfb(i,j,k)          = Param.Dcell*(1-CellState.p(cnt)); % defining diffusion coefficient for each cell

                       CellState.Pop(1) =  CellState.Pop(1) + 1; % counts total cells in the system
                       cnt       = cnt + 1;       % increase cell count               
                    end
                end
            end
        end
        
        
     case 2
        % Initializing epithelial cells as a straight hollow cylinder -
        % HORIZONTAL
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    if ((j-h)^2 + (k-h)^2 <= Param.Crado^2) && ((j-h)^2 + (k-h)^2 >= Param.Cradi^2) % ((i-h)^2 + (j-h)^2 + (k-h)^2 <= Crado^2) && ((i-h)^2 + (j-h)^2 + (k-h)^2 >= Cradi^2)
                       CellState.Position(cnt, :) = [i j k];   % add cells
                       CellState.Divide(cnt)     = 0;
                       CellState.conc(cnt, :)    = EConc; % initial ODE conditions
                                % [snail, SNAIL, miRNA34, zeb1, Zeb1, miRNA200, ECAD, NCAD]
                       CellState.state(cnt)      = 1;
  
                       CellState.TimesMoved(cnt) = 0;

                       CellState.p(cnt)      = CellState.conc(cnt, 8) / Param.NcadMax; % probability threshold dependent on Ncad
                       Dtgfb(i,j,k)          = Param.Dcell*(1-CellState.p(cnt)); % defining diffusion coefficient for each cell

                       CellState.Pop(1) =  CellState.Pop(1) + 1; % counts total cells in the system
                       cnt       = cnt + 1;       % increase cell count               
                    end
                end
            end
        end 
    
        
    case 3
        % Initializing Cell Structure: Curved Tubule
        xp = 1:n-3; 
        yp = [11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,13,14,15,16,18,19,22,25,28];
        zp = (n/2)*ones(1, length(xp));

        curve = [xp; yp; zp];   np = Param.Crado;  npoints = 1;
        for k = 2:(size(curve,2)-1)
            if norm(curve(:,k) - curve(:,npoints)) > Param.Crado/2
              npoints = npoints + 1;
              curve(:,npoints) = curve(:,k);
            end
        end

        % Always include endpoint
        if norm(curve(:,end) - curve(:,npoints)) > 0
            npoints = npoints+1;
            curve(:,npoints) = curve(:,end);
        end

        % deltavecs: average for internal points. first strecth for endpoitns.
        dv = curve(:,[2:end,end]) - curve(:,[1,1:end-1]);

        % make nvec not parallel to dv(:,1)
        nvec = zeros(3,1); [buf,idx] = min(abs(dv(:,1))); 
        nvec(idx) = 1; xyz = repmat([0],[3,n + 1,npoints + 2]);

        % precalculate cos and sing factors:
        cfact = repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
        sfact = repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
        
        % Main loop: propagate the normal (nvec) along the tube
        for k = 1:npoints
            convec = cross(nvec,dv(:,k));   convec = convec./norm(convec);
            nvec = cross(dv(:,k),convec);   nvec = nvec./norm(nvec);
            %update xyz:
            xyz(:,:,k+1) = repmat(curve(:,k),[1,n+1])+ cfact.*repmat(Param.Crado*nvec,[1,n+1])...
                + sfact.*repmat(Param.Crado*convec,[1,n+1]);
        end

        %finally, cap the ends:
        xyz(:,:,1) = repmat(curve(:,1),[1,n+1]);
        xyz(:,:,end) = repmat(curve(:,end),[1,n+1]);
        xyz = round(xyz);

        % extract results:
        xt = squeeze(xyz(1,:,:)); yt = squeeze(xyz(2,:,:)); zt = squeeze(xyz(3,:,:));
        xt = reshape(xt, [495, 1]); yt = reshape(yt, [495, 1]); zt = reshape(zt, [495, 1]);
        curvet = unique([xt, yt, zt], 'rows');

        for i = 1:length(curvet)
           CellState.Position(cnt, :) = curvet(i, :);   % add cells
           CellState.Divide(cnt)     = 0;
           CellState.conc(cnt, :)    = EConc; % initial ODE conditions
           
           CellState.state(cnt)      = 1;
           CellState.TimesMoved(cnt) = 0;

           CellState.p(cnt)      = CellState.conc(cnt, 8) / Param.NcadMax; % probability threshold dependent on Ncad
           Dtgfb(curvet(i, 1), curvet(i, 2), curvet(i, 3)) = Param.Dcell*(1-CellState.p(cnt)); % defining diffusion coefficient for each cell

           CellState.Pop(1)      = CellState.Pop(1) + 1; % counts total cells in the system
           cnt                   = cnt + 1;       % increase cell count  
        end

end


% If they want to set up a co-culture model
if contains(CultureModel, 'Epithelial Only')
    
else
    for v = 1:FibroblastCount  
        if contains(CultureModel, 'Embedded')
            pos = randsample(n, 3);
        elseif contains(CultureModel, 'Overlay')
            pos = randsample(n, 2); pos(3) = n - 2;
        end
        
        if ((pos(1)-h)^2 + (pos(2)-h)^2 + (pos(3)-h)^2 > Param.Crado^2)
           CellState.Position(cnt, :) = [pos(1) pos(2) pos(3)];   % add cells
           CellState.Divide(cnt)     = 0;
           CellState.conc(cnt, :)    = MConc; % initial ODE conditions
           CellState.state(cnt)      = 4;
    
           CellState.TimesMoved(cnt) = 0;     
           
           CellState.p(cnt) = CellState.conc(cnt, 8) / Param.NcadMax; 
           Dtgfb(pos(1), pos(2), pos(3)) = Param.Dcell*(1 - CellState.p(cnt));
           
           CellState.Pop(4) =  CellState.Pop(4) + 1; % counts total cells in the system
           cnt = cnt + 1; % increase cell count 
           
        end
    end
end 


%% Create a nxnxn matrix that tracks cell movement + EMT changes
Cstate = zeros(n,n,n); 

for i = 1:length(CellState.state)
   [posx] =  CellState.Position(i, :);
   Cstate(posx(1), posx(2), posx(3)) = CellState.state(i);
end

% Update CellState Structure based on initialized morphology
CellState.Position = [CellState.Position, CellState.Position]; 

CellState.Ncad = [CellState.conc(:, 8), CellState.conc(:, 8)];
CellState.Ecad         = CellState.conc(:, 7);

CellState.Snail        = CellState.conc(:, 2);
CellState.Zeb1         = CellState.conc(:, 5); 

CellState.R200        = CellState.conc(:, 6);
CellState.R34         = CellState.conc(:, 3);   

CellState.AvgEcad(1)       = mean(CellState.conc(:, 7)); 
CellState.AvgNcad(1)       = mean(CellState.conc(:, 8));

CellState.AvgSnail(1)      = mean(CellState.conc(:, 2)); 
CellState.AvgZeb1(1)       = mean(CellState.conc(:, 5));

CellState.AvgR200(1)       = mean(CellState.conc(:, 6)); 
CellState.AvgR34(1)        = mean(CellState.conc(:, 3));

Imgdim = sum(Cstate); Imgdim = squeeze(Imgdim); Imgdim = imfill(Imgdim);
Imedge = edge(Imgdim); 

CellState.Carea(1) = length(find(Imgdim > 0));
CellState.Cperim(1) = length(find(Imedge > 0));
    
end

