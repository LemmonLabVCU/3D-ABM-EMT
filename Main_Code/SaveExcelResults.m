function DataOutput = SaveExcelResults(CellState, Param, Cstate, Ctgfb, Dtgfb, incr, DataOutput)
% SaveExcelResults.m: Saves outputs to an excel file


% rename commonly used variables:
n           = Param.n;               % grid size (units)
L           = n*Param.Csize;          % length of the mesh (um)
h           = round(n/2);       % center of the spheroid (pixels)

% Quantify changes in cell movement
for i = 1:length(CellState.Position)
    P1 = CellState.Position(i, 1:3); 
    P2 = CellState.Position(i, 4:6);   
    DiffMove = sqrt((P2(1) - P1(1))^2 + (P2(2) - P1(2))^2 + (P2(3) - P1(3))^2);
    CellState.movement(i) = DiffMove; 
end

% Quantify Changes in Cell Position
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


% Quantify Apical v Basal TGFB and Dtgfb Expression
TGFB_Apical = 0; Dtgfb_Apical = 0; ApicalCnt = 1;
TGFB_Basal = 0; Dtgfb_Basal = 0; BasalCnt = 1;

for i = 1:n
        for j = 1:n
            for k = 1:n
                if Cstate(i,j,k) > 0 && Cstate(i,j,k) < 4 % Apical
                    TGFB_Apical = TGFB_Apical + Ctgfb(i,j,k);
                    Dtgfb_Apical = Dtgfb_Apical + Dtgfb(i,j,k);
                    ApicalCnt = ApicalCnt + 1;

                elseif ((i-h)^2 + (j-h)^2 + (k-h)^2 <= Param.Crado^2) && ((i-h)^2 + (j-h)^2 + (k-h)^2 >= Param.Cradi^2) % inside spheroid 
                    TGFB_Apical = TGFB_Apical + Ctgfb(i,j,k);
                    Dtgfb_Apical = Dtgfb_Apical + Dtgfb(i,j,k);
                    ApicalCnt = ApicalCnt + 1;

                elseif ((i-h)^2 + (j-h)^2 + (k-h)^2 >= Param.Crado^2) % Basal
                    TGFB_Basal = TGFB_Basal + Ctgfb(i,j,k);
                    Dtgfb_Basal = Dtgfb_Basal + Dtgfb(i,j,k);
                    BasalCnt = BasalCnt + 1;
                end
            end
        end
end

% Quantify Spheroid Morphology
[CellState] = UpdateSpheroidMorphology(CellState, Cstate);


% Organizing Data into DataOutput variable:

eInd = find(CellState.state < 4); % find values for epithelial cells

DataOutput.TotState(incr, :)  = [sum(CellState.Pop(end, 1:3)), CellState.Pop(end, 1:4)];
DataOutput.CellLoc(incr, 1:3) = CellLoc;
  
DataOutput.EcadAvg(incr, 1)   = mean(CellState.conc(eInd, 7));
DataOutput.NcadAvg(incr, 1)   = mean(CellState.Ncad(eInd, 2));

DataOutput.snailTAvg(incr, 1) = mean(CellState.conc(eInd, 1));
DataOutput.SnailAvg(incr, 1)  = mean(CellState.conc(eInd, 2));

DataOutput.AvgzebT(incr, 1)   = mean(CellState.conc(eInd, 4));
DataOutput.AvgZeb1(incr, 1)    = mean(CellState.conc(eInd, 5));
    
DataOutput.AvgR200(incr, 1)   = mean(CellState.conc(eInd, 6)); 
DataOutput.AvgR34(incr, 1)    = mean(CellState.conc(eInd, 3));

DataOutput.CellTGFBAvg(incr, 1)   = mean(Ctgfb, 'all');

DataOutput.TGFBApicalAvg(incr, 1) = TGFB_Apical/ ApicalCnt;
DataOutput.TGFBBasalAvg(incr, 1) = TGFB_Basal / BasalCnt;

DataOutput.TGFBApicalSum(incr, 1) = TGFB_Apical;
DataOutput.TGFBBasalSum(incr, 1) = TGFB_Basal;

DataOutput.DtgfbAvg(incr, 1)      = mean(Dtgfb, 'all');

DataOutput.DtgfbApicalAvg(incr, 1) = Dtgfb_Apical / ApicalCnt;
DataOutput.DtgfbBasalAvg(incr, 1) = Dtgfb_Basal / BasalCnt;

DataOutput.DtgfbApicalSum(incr, 1) = Dtgfb_Apical;
DataOutput.DtgfbBasalSum(incr, 1) = Dtgfb_Basal;

DataOutput.DistanceMoved(incr, 1) = mean(CellState.movement(eInd));
DataOutput.AvgTimesMoved(incr, 1) = mean(CellState.TimesMoved(eInd));

DataOutput.CArea(incr, 1) = CellState.Carea(end);
DataOutput.LumenArea(incr, 1) = CellState.LumenArea(end);
DataOutput.CSCellCount(incr, 1) = CellState.CSCount(end);

DataOutput.MajDia(incr, 1) = CellState.MajorD(end);
DataOutput.MinDia(incr, 1) = CellState.MinorD(end);

DataOutput.TotalVolume(incr, 1) = CellState.TotalVolume(end);

for i = 2:5
    ind = find(CellState.state == (i-1));
    ind2 = find(Cstate == i - 1);         
    
    if ind > 0 % defining outputs of individual cell states
        DataOutput.EcadAvg(incr, i) = mean(CellState.conc(ind, 7));
        DataOutput.NcadAvg(incr, i) = mean(CellState.Ncad(ind, 2));

        DataOutput.snailTAvg(incr, i) = mean(CellState.conc(ind, 1));
        DataOutput.SnailAvg(incr, i)  = mean(CellState.conc(ind, 2));
        
        DataOutput.AvgzebT(incr, i)   = mean(CellState.conc(ind, 4));
        DataOutput.AvgZeb1(incr, i)    = mean(CellState.conc(ind, 5));

        DataOutput.AvgR200(incr, i)   = mean(CellState.conc(ind, 6)); 
        DataOutput.AvgR34(incr, i)    = mean(CellState.conc(ind,3));
        
        DataOutput.CellTGFBAvg(incr, 1)   = mean(CellState.Ctgfb(ind));
        DataOutput.CellDtgfb(incr, 1)     = mean(CellState.Dcell(ind));
                
        DataOutput.DistanceMoved(incr, i) = mean(CellState.movement(ind));
        DataOutput.AvgTimesMoved(incr, i) = mean(CellState.TimesMoved(ind));

    else % if there are no cells in specified state
        DataOutput.EcadAvg(incr, i)       = 0;
        DataOutput.NcadAvg(incr, i)       = 0;

        DataOutput.snailTAvg(incr, i)     = 0;
        DataOutput.SnailAvg(incr, i)      = 0;
        
        DataOutput.AvgzebT(incr, i)       = 0;
        DataOutput.AvgZeb1(incr, i)       = 0;

        DataOutput.AvgR200(incr, i)       = 0; 
        DataOutput.AvgR34(incr, i)        = 0;

        DataOutput.CellTGFBAvg(incr, 1)   = 0;
        DataOutput.CellDtgfb(incr, 1)     = 0;

        DataOutput.DistanceMoved(incr, i) = 0;
        DataOutput.AvgTimesMoved(incr, i) = 0; 

    end             
end  % Finish formatting data

end
