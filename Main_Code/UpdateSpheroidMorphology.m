function [CellState] = UpdateSpheroidMorphology(CellState, Cstate)

% Quantify Changes in Spheroid Morphology:
% 1. Max Cross-Sectional Area
% 2. Max Lumen Area
% 3. Total Volume
% 4. Cell Count at Max CS Area


% Remove Fibroblasts if present
fInd = find(CellState.state == 4);

Cstate2 = Cstate;
for nn = 1:length(fInd)
    [posx] =  CellState.Position(fInd(nn), 4:6);
    Cstate2(posx(1), posx(2), posx(3)) = 0;
end

% Convert to binary matrix
Cstate2(Cstate2 ~= 0) = 1;
tempCS = imfill(Cstate2);

% Flatten 3D Matrix
Imgdim = sum(Cstate2); Imgdim = squeeze(Imgdim);     

Imgdim2 = permute(Imgdim, [1 3 2]); 
Imgdim2 = reshape(Imgdim,[],size(Imgdim,2),1);
Imgdim2 = imbinarize(Imgdim2);
Imgdim3 = edge(Imgdim2);

tempStat = regionprops(Imgdim2, 'Area', 'MajorAxisLength','MinorAxisLength', 'Perimeter');

tempStat = struct2table(tempStat);
[~, maxInd] = max(tempStat.Area);

% Quantify Morphological Measurement (pixel)
CellState.Carea = tempStat.Area(maxInd);
CellState.MajorD = tempStat.MajorAxisLength(maxInd);
CellState.MinorD = tempStat.MinorAxisLength(maxInd);

CellState.CSCount = sum(Imgdim3, 'all'); 
CellState.LumenArea = sum(Imgdim2 - Imgdim3, 'all'); 
CellState.TotalVolume = (4/3)*pi*CellState.MajorD*CellState.MinorD;

end