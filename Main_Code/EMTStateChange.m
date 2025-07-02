function [CellState, Cstate, Ctgfb, Dtgfb] = EMTStateChange(CellState, Param, Cstate, Ctgfb, Dtgfb, t, Tfinal, dt)
% EMTStateChange: Function that runs EMT characterization based on changes in concentration of N-Cadherin. ODE is solved
%using Euler Method.

numCells = length(CellState.state);
CellState.Pop(t, :) = [0, 0, 0, 0]; % CellState.Pop(t - 1, :);
    
for v = 1:numCells % for the number of cells in the system
        
    if CellState.state(v) > 0
        i = CellState.Position(v,4); 
        j = CellState.Position(v,5); 
        k = CellState.Position(v,6);
    
        
        [Ctgfb(i,j,k), CellState.conc(v,:)] = TGFB_ODE(Param, t, dt:dt:Tfinal, dt, CellState.conc(v,:), Ctgfb(i,j,k));
        
        CellState.Ncad(v, 2) = CellState.conc(v,8);
        CellState.p(v)       = CellState.Ncad(v,2) / Param.NcadMax;

        if CellState.p(v) > 1
            CellState.p(v) = 1;
        end

        Dtgfb(i,j,k)         = Param.Dcell*(1- (Param.EMTThresh*CellState.p(v)));
        CellState.Dcell(v)   = Dtgfb(i,j,k);
        CellState.Ctgfb(v)   = Ctgfb(i,j,k);

        if CellState.state(v) < 4
            % Changing Cell State
            % E: < 1.5;    P: 1.5 < X < 2.4;  M: 2.4 < X < 3.1515
            if CellState.Ncad(v,2) >= 2.4 || CellState.state(v) == 3 % Mesenchymal State
                CellState.state(v) = 3;
                Cstate(i,j,k) = 3;
        
            elseif CellState.Ncad(v,2) >= 1.5 && CellState.Ncad(v,2) < 2.4 % Partial State
                CellState.state(v) = 2;
                Cstate(i,j,k) = 2;
        
            elseif CellState.Ncad(v,2) < 1.5 % Epithelial State
                CellState.state(v) = 1;
                Cstate(i,j,k) = 1;
        
            end

        tempInd = CellState.state(v);
        CellState.Pop(t, tempInd) = CellState.Pop(t, tempInd) + 1;
        
        end   
    end
end

CellState.Ncad(:,1) = CellState.Ncad(:,2); % Update Ncad concentrations
  
end

