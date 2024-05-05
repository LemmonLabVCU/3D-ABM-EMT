function DataOutput = SaveExcelResults(CellState, Param, Cstate, Ctgfb, C0tgfb, TGFB, Dtgfb, Tfinal, incr)

n           = Param.n;               % grid size (units)
L           = n*Param.Csize;          % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)
DataOutput  = [];

% Total Run Outputs
Change_TGFB = Ctgfb - C0tgfb; % Change in TGFB concentration

% Changes in cell movement
for i = 1:length(CellState.Position)
    P1 = CellState.Position(i, 1:3); 
    P2 = CellState.Position(i, 4:6);   
    DiffMove = sqrt((P2(1) - P1(1))^2 + (P2(2) - P1(2))^2 + (P2(3) - P1(3))^2);
    CellState.movement(i) = DiffMove; 
end

% Changing Cell Position in Spheroid structure
FPos 	= CellState.Position(:, 4:6); % For easier indexing
CellLoc = [0 0 0]; % InitPos, Inside, Outside

for v = 1:length(FPos)
    if ((FPos(v, 1)-h)^2 + (FPos(v, 2)-h)^2 + (FPos(v, 3)-h)^2 <= Param.Crado^2) && ((FPos(v, 1)-h)^2 + (FPos(v, 2)-h)^2 + (FPos(v, 3)-h)^2 == Param.Cradi^2)
        CellLoc(1) = CellLoc(1) + 1; % On spheroid edge
        
    elseif ((FPos(v, 1)-h)^2 + (FPos(v, 2)-h)^2 + (FPos(v, 3)-h)^2 < Param.Cradi^2) % Inside Lumen
        CellLoc(2) = CellLoc(2) + 1; 
        
    else
        CellLoc(3) = CellLoc(3) + 1;  % outside spheroid area
    end
end

DataOutput.TGFB(incr,1) = TGFB;
DataOutput.Dtgfb(incr,1) = Param.Dcell;
DataOutput.Tfinal(incr, 1) = Tfinal;


eInd = find(CellState.state < 4); % find values for epithelial cells

DataOutput.TotState(incr, :)  = [sum(CellState.Pop(end, 1:3)), CellState.Pop(end, 1:4)];
DataOutput.CellLoc(incr, 1:3) = CellLoc;
  
DataOutput.EcadAvg(incr, 1)  = mean(CellState.Ecad(eInd));
DataOutput.NcadAvg(incr, 1)  = mean(CellState.Ncad(eInd,2));

DataOutput.SnailAvg(incr, 1) = mean(CellState.Snail(eInd));    
DataOutput.Zeb1Avg(incr, 1)  = mean(CellState.Zeb1(eInd));

DataOutput.AvgmiR34(incr, 1) = mean(CellState.conc(eInd, 3));
DataOutput.AvgmiR200(incr, 1) = mean(CellState.conc(eInd, 6));

DataOutput.CellTgfbAvg(incr, 1) = mean(mean(mean(Ctgfb)));
DataOutput.TgfbAvgChange(incr, 1) = mean(mean(mean(Change_TGFB))); 

DataOutput.CellDtgfbAvg(incr, 1) = mean(mean(mean(Dtgfb)));     

DataOutput.DistanceMoved(incr, 1) = mean(CellState.movement(eInd));
DataOutput.AvgTimesMoved(incr, 1) = mean(CellState.TimesMoved(eInd));

DataOutput.SpheroidArea(incr, 1) = CellState.Carea(end);
DataOutput.SpheroidPerim(incr, 1) = CellState.Cperim(end);

for i = 2:5
    ind = find(CellState.state == (i-1));
    ind2 = find(Cstate == i - 1);         
    
    if ind > 0 % defining outputs of individual cell states
        DataOutput.EcadAvg(incr, i) = mean(CellState.Ecad(ind));
        DataOutput.NcadAvg(incr, i) = mean(CellState.Ncad(ind, 2));

        DataOutput.SnailAvg(incr, i)      = mean(CellState.Snail(ind)); 
        DataOutput.Zeb1Avg(incr, i)       = mean(CellState.Zeb1(ind));
        
        DataOutput.AvgmiR34(incr, i)      = mean(CellState.conc(ind, 3));
        DataOutput.AvgmiR200(incr, i)     = mean(CellState.conc(ind, 6));

        DataOutput.CellTgfbAvg(incr, i)   = mean(Ctgfb(ind2));
        DataOutput.TgfbAvgChange(incr, i) = mean(Change_TGFB(ind2)); 
             
        DataOutput.CellDtgfbAvg(incr, i)  = mean(Dtgfb(ind2));
        
        DataOutput.DistanceMoved(incr, i) = mean(CellState.movement(ind));
        DataOutput.AvgTimesMoved(incr, i) = mean(CellState.TimesMoved(ind));

    else % if there are no cells in specified state
        DataOutput.EcadAvg(incr, i)       = 0;
        DataOutput.NcadAvg(incr, i)       = 0;
        
        DataOutput.SnailAvg(incr, i)      = 0;
        DataOutput.Zeb1Avg(incr, i)       = 0;
        
        DataOutput.AvgmiR34(incr, i)      = 0;
        DataOutput.AvgmiR200(incr, i)     = 0;

        DataOutput.CellTgfbAvg(incr, i)   = 0;
        DataOutput.TgfbAvgChange(incr, i) = 0;
             
        DataOutput.CellDtgfbAvg(incr, i)  = 0;
        
        DataOutput.DistanceMoved(incr, i) = 0;
        DataOutput.AvgTimesMoved(incr, i) = 0; 
        
        DataOutput.TgfbAvgChange(incr, i) = 0;
    end             
end

end

