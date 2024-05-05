function [CellState, Cstate, Ctgfb, Dtgfb] = EMTStateChange(CellState, Param, Cstate, Ctgfb, Dtgfb, t, Tfinal, dt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

numCells = length(CellState.state);
CellState.Pop(t, :) = CellState.Pop(t - 1, :);
    
for v = 1:numCells % for the number of cells in the system
        
    if CellState.state(v) > 0
        i = CellState.Position(v,4); 
        j = CellState.Position(v,5); 
        k = CellState.Position(v,6);
    
        [Ctgfb(i,j,k), CellState.conc(v,:)] = TGFB_ODE(t, dt:dt:Tfinal, dt, CellState.conc(v,:), Ctgfb(i,j,k));
        
        CellState.Ncad(v, 2) = CellState.conc(v,8);
        CellState.Ecad(v) = CellState.conc(v,7);
        CellState.Snail(v) = CellState.conc(v,2);
        
        CellState.p(v)       = CellState.Ncad(v,2) / Param.NcadMax;
        Dtgfb(i,j,k)         = Param.Dcell*(1-CellState.p(v));
        CellState.Ctgfb(v)   = Ctgfb(i,j,k);

      % Changing Cell State
            % E: < 1.5;    P: 1.5 < X < 2.4;  M: 2.4 < X < 3.1515
            if CellState.Ncad(v,2) >= 2.4 && CellState.state(v) < 4 % Mesenchymal State
                CellState.state(v) = 3;
                Cstate(i,j,k) = 3;

                if CellState.Ncad(v,1) < 2.4 && CellState.Ncad(v,2) >= 2.4 % Partial to Mesenchymal
                    CellState.Pop(t, 2) =  CellState.Pop(t ,2) - 1;
                    CellState.Pop(t, 3) =  CellState.Pop(t, 3) + 1;
                end

            elseif CellState.Ncad(v,2) >= 1.5 && CellState.Ncad(v,2) < 2.4 && CellState.state(v) < 4% Partial State
                CellState.state(v) = 2;
                Cstate(i,j,k) = 2;

                if CellState.Ncad(v,1) < 1.5 % E to P
                    CellState.Pop(t, 1) = CellState.Pop(t, 1) - 1;
                    CellState.Pop(t, 2) = CellState.Pop(t, 2) + 1;
                end

            elseif CellState.Ncad(v,2) < 1.5 && CellState.state(v) < 4 % Epithelial State
                CellState.state(v) = 1;
                Cstate(i,j,k) = 1;
                
                if CellState.Ncad(v,1) >= 1.5 && CellState.Ncad(v,1) < 2.4 % P to E
                    CellState.Pop(t, 1) = CellState.Pop(t, 1) + 1;
                    CellState.Pop(t, 2) = CellState.Pop(t, 2) - 1; 
                end         
            end
     end           
end

CellState.Ncad(:,1) = CellState.Ncad(:,2); % Update Ncad concentrations
  
end

