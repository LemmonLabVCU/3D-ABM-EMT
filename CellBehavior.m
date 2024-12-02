% FUNCTION: CellBehavior.m
% Goal: Defines a 3D moore neighborhood around the cell in question at
% position (i,j,k). Inputs TGFB gradients and neighboring cells within the
% defined 3x3x3 window that is used to determine cell behavior. 

% [CellState, Cstate, Cmoved, Ctgfb, Dtgfb] = CellBehavior(
% CellState, Cstate, i, j, k, v, ProlifThreshold, Ctgfb, t, Cmoved, Dtgfb)

% INPUTS: 
% _______________________________________________________________________
%   CellState:  Structure that holds the characteristics of all cells
%   defined in the system.

%   Cstate: nxnxn matrix that defines cell position + EMT state

%   [i,j,k]: Center position in the 3x3x3 neighborhood. 

%   v:  index value of the cell being modified in the CellState structure

%   ProlifThreshold: Expected time for cell to proliferate 
%                   (based on cell type being modeled)       

%   Ctgfb:  TGFB concentration (nxnxn matrix)

%   Cmoved: value that states if a cell has moved or not to prevent
%   repeated cells.

%   Dtgfb: TGFB diffusion coefficient matrix (nxnxn)

% OUTPUTS: 
% _________________________________________________________________________

% __________________________________________________________________________

% NOTES:

% _________________________________________________________________________

function [CellState, Cstate, Cmoved, Ctgfb, Dtgfb] = CellBehavior(CellState,...
    Param, Cstate, i, j, k, v, Ctgfb, t, Cmoved, Dtgfb) 
       
% Defining the Moore Neighborhood: 
%       (3x3x3 centered around the cell at position (i,j,k))

n = Param.n; % boundaries of the system
localcells = 4 + zeros(1,27);
Tgrad = zeros(1, 27);

localcells(14) = 0; % we can assume it can stay in one spot? 

if i ~= 1 
    localcells(17) = Cstate(i-1,j,k);
    Tgrad(17)      = Ctgfb(i-1,j,k);
    
    if j ~= 1
        localcells(15) = Cstate(i,j-1,k);
        Tgrad(15)      = Ctgfb(i-1,j,k);
        localcells(18) = Cstate(i-1,j-1,k);
        Tgrad(18)      = Ctgfb(i-1,j-1,k);
        
        if k ~= 1
            localcells(23) = Cstate(i,j,k-1);
            Tgrad(23)      = Ctgfb(i,j,k-1);
            localcells(24) = Cstate(i,j-1,k-1);
            Tgrad(24)      = Ctgfb(i, j-1, k-1);
            localcells(26) = Cstate(i-1,j,k-1);
            Tgrad(26)      = Ctgfb(i-1,j,k-1);
            localcells(27) = Cstate(i-1, j-1,k-1);
            Tgrad(27)      = Ctgfb(i-1, j-1,k-1);
        end
        
        if k ~= n
            localcells(5) = Cstate(i,j,k+1);
            Tgrad(5)      = Ctgfb(i,j,k+1);
            localcells(6) = Cstate(i,j-1,k+1);
            Tgrad(6)      = Ctgfb(i,j-1,k+1);
            localcells(9) = Cstate(i-1,j-1,k+1);
            Tgrad(9)      = Ctgfb(i-1,j-1,k+1);
        end
    end

    if j ~= n
        localcells(13)  = Cstate(i, j+1, k);
        Tgrad(13)       = Ctgfb(i, j+1, k);
        localcells(16)  = Cstate(i-1, j+1, k);
        Tgrad(16)       = Ctgfb(i-1, j+1, k);
        
        if k ~= 1
            localcells(23)  = Cstate(i,j,k-1);
            Tgrad(23)       = Ctgfb(i,j,k-1);
            localcells(22)  = Cstate(i,j+1,k-1);
            Tgrad(22)       = Ctgfb(i, j+1, k-1);
            localcells(26)  = Cstate(i-1,j,k-1);
            Tgrad(26)       = Ctgfb(i-1,j,k-1);
            localcells(25)  = Cstate(i-1,j+1,k-1);
            Tgrad(25)       = Ctgfb(i-1, j+1, k-1);
        end
        
        if k~= n
            localcells(5) = Cstate(i,j,k+1);
            Tgrad(5)      = Ctgfb(i,j, k+1);
            localcells(4) = Cstate(i, j+1, k+1);
            Tgrad(4)      = Ctgfb(i, j+1, k+1);
            localcells(8) = Cstate(i-1,j,k+1);
            Tgrad(8)      = Ctgfb(i-1,j,k+1);
            localcells(7) = Cstate(i-1,j+1,k+1);
            Tgrad(7)      = Ctgfb(i-1,j+1, k+1);
        end
    end
end

if i ~= n
    localcells(11) = Cstate(i+1,j,k);
    Tgrad(11)      = Ctgfb(i+1,j,k);

    if j ~= 1
        localcells(15) = Cstate(i,j-1,k);
        Tgrad(15)      = Ctgfb(i,j-1,k);
        localcells(12) = Cstate(i+1,j-1,k);
        Tgrad(12)      = Ctgfb(i+1,j-1,k);
        
        if k ~= 1
            localcells(23) = Cstate(i,j,k-1);
            Tgrad(23)      = Ctgfb(i,j,k-1);
            localcells(24) = Cstate(i,j-1,k-1);
            Tgrad(24)      = Ctgfb(i,j-1,k-1);
            localcells(20) = Cstate(i+1,j,k-1);
            Tgrad(20)      = Ctgfb(i+1,j,k-1);
            localcells(21) = Cstate(i+1, j-1,k-1);
            Tgrad(21)      = Ctgfb(i+1,j-1,k-1);
        end
        
        if k ~= n
            localcells(5) = Cstate(i,j,k+1);
            Tgrad(5)      = Ctgfb(i,j,k+1);
            localcells(6) = Cstate(i,j-1,k+1);
            Tgrad(6)      = Ctgfb(i,j-1,k+1);
            localcells(2) = Cstate(i+1,j,k+1);
            Tgrad(2)      = Ctgfb(i+1,j,k+1);
            localcells(3) = Cstate(i+1,j-1,k+1);
            Tgrad(3)      = Ctgfb(i+1,j-1,k+1);
        end
    end

    if j ~= n
        localcells(13) = Cstate(i, j+1, k);
        Tgrad(13)      = Ctgfb(i,j+1,k);
        localcells(10) = Cstate(i+1, j+1, k);
        Tgrad(10)      = Ctgfb(i+1,j+1,k);
        
        if k ~= 1
            localcells(23) = Cstate(i,j,k-1);
            Tgrad(23)      = Ctgfb(i,j,k-1);
            localcells(22) = Cstate(i,j+1,k-1);
            Tgrad(22)      = Ctgfb(i,j+1,k-1);
            localcells(20) = Cstate(i+1,j,k-1);
            Tgrad(20)      = Ctgfb(i+1,j,k-1);
            localcells(19) = Cstate(i+1,j+1,k-1);
            Tgrad(19)      = Ctgfb(i+1,j+1,k-1);
        end
        
        if k ~= n
            localcells(5) = Cstate(i,j,k+1);
            Tgrad(5)      = Ctgfb(i,j,k+1);
            localcells(4) = Cstate(i, j+1, k+1);
            Tgrad(4)      = Ctgfb(i,j+1,k+1);
            localcells(2) = Cstate(i+1,j,k+1);
            Tgrad(2)      = Ctgfb(i+1,j,k+1);
            localcells(1) = Cstate(i+1,j+1,k+1);
            Tgrad(1)       = Ctgfb(i+1,j+1,k+1);
        end
    end
end

%% Find empty sites in the grid based on Cstate values
emptyLoc = [];
emptyTGFB = [];

movement = [1 1 1; 1 0 1; 1 -1 1; 0 1 1; 0 0 1; 0 -1 1; -1 1 1; -1 0 1; ...
    -1 -1 1; 1 1 0; 1 0 0; 1 -1 0; 0 1 0; 0 0 0; 0 -1 0; -1 1 0; -1 0 0; ...
    -1 -1 0; 1 1 -1; 1 0 -1; 1 -1 -1; 0 1 -1; 0 0 -1; 0 -1 -1; -1 1 -1; ...
    -1 0 -1; -1 -1 -1]; 

for nn = 1:length(localcells)
    
    if localcells(nn) == 0
        emptyLoc  = [emptyLoc; movement(nn,:)];
        emptyTGFB = [emptyTGFB; Tgrad(nn)];
    end
end

numEmpty = size(emptyLoc,1); % number of free spots with no cells

%% Changes in Cell Behavior

% Cell Migration: based on probability thresholds

    if CellState.state(v) > 0 && rem(t, Param.ProlifThresh) ~= 0 % make sure cell is present
        
        if CellState.state(v) == 1 % epithelial cells don't move
            threshold =  Param.PThresh(1);
        elseif CellState.state(v) == 2 % assume pEMT cells don't move 
            threshold =  Param.PThresh(2); 
        else
            threshold =  Param.PThresh(3); % mesenchymal cells move
        end 
        
        if (rand(1) < threshold)
            
            maxTGFB        = max(emptyTGFB);
            CellState.TimesMoved(v) = CellState.TimesMoved(v) + 1;
            ind            = find(emptyTGFB == maxTGFB);
            
            if length(ind) > 1
               ind = ind(randi([1, length(ind)], 1));
            end

            newi = i + emptyLoc(ind, 1);
            newj = j + emptyLoc(ind, 2);
            newk = k + emptyLoc(ind, 3); 

           if i == newi && k == newk && k == newk
              Cmoved(i,j,k) = 1; 
              
           else
                CellState.Position(v, 4:6) = [newi, newj, newk]; 
                Cstate(newi, newj, newk) = Cstate(i,j,k); 
                Cstate(i,j,k)            = 0;
                
                Dtgfb(newi, newj, newk)  = Dtgfb(i,j,k); % moving ECM along, not breaking previous
                CellState.Ctgfb(v)       = Ctgfb(newi, newj, newk);
                CellState.Dcell(v)       = Dtgfb(newi, newj, newk);
                
                Cmoved(newi, newj, newk) = 1;
           end
        end
        
    % Cell Proliferation: Expected to happen at certain time points/ based on probability
    
    elseif rem(t,  Param.ProlifThresh) == 0
    % 1. Set Threshold Probabilities
         if CellState.state(v) == 1 % epithelial will proliferate
            threshold =  Param.PThresh(3);
            cellType = 1;
        elseif CellState.state(v) == 2 % pEMT still proliferates 
            threshold =  Param.PThresh(2);
            cellType = 2;
         elseif CellState.state(v) == 3
            threshold =  Param.PThresh(1); % mesenchymal cells don't proliferate
            cellType = 3;
         else
             threshold =  Param.PThresh(1);
             cellType = 4;
         end 
    
    if (rand(1) < threshold) % It will proliferate if threshold is reached
        CellState.Divide(v) = CellState.Divide(v) + 1; % increase counter of proliferation
        
        % Define new cell as an exact copy of the original
        siteRand = rand*numEmpty;   
        
        if numEmpty > 0
            siteInd = max(ceil(siteRand),1);
            newi    = i + emptyLoc(siteInd,1);
            newj    = j + emptyLoc(siteInd, 2);
            newk    = k + emptyLoc(siteInd,3);

            Cstate(newi, newj, newk)        = Cstate(i,j,k);
            CellState.Position(end + 1,:)    = [newi, newj, newk, newi, newj, newk];
            
            CellState.state(end + 1)        = CellState.state(v);
            CellState.conc(end + 1, :)      = CellState.conc(v, :);
            
            CellState.Ncad(end + 1, :)      = [CellState.Ncad(v, 2), CellState.Ncad(v, 2)];
            
            CellState.p(end + 1)            = CellState.p(v);
            CellState.Divide(end + 1)       = 0;
            CellState.TimesMoved(end + 1)   = 0; 

            Dtgfb(newi,newj,newk)           = Dtgfb(i,j,k); % reuse same ECM type
            
            CellState.Dcell(end + 1)        = Dtgfb(newi, newj, newk);
            CellState.Ctgfb(end + 1)        = Ctgfb(newi, newj, newk);
    
            CellState.Pop(t, cellType) = CellState.Pop(t, cellType) + 1; 
                % assume new cell a replica of the original 

            Cmoved(i,j,k) = 1; % cell proliferated / not expected to move
            Cmoved(newi, newj, newk) = 1; % cell technically moved
        end
        
    end
    
    end

end
