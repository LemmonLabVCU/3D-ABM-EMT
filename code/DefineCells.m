% FUNCTION: DefineCells.m
% Goal: Set up the CellState structure that stores cell parameters changing
% over time. 

% [CellState] = DefineCells(CellState)
% _______________________________________________________________________


function [CellState] = DefineCells(CellState)
% DefineCells: Define the parameters defining each cell stated in the
% system

CellState.Divide  = []; % Number of cell divisions

CellState.conc = []; % Initial ODE values of each cell
CellState.Ncad = [];
CellState.Ecad = [];
CellState.Snail = [];
CellState.Zeb1  = [];

CellState.state = []; % defines each cell state (1 = E, 2 = P, 3 = M) 
CellState.Pop = [0, 0, 0, 0]; % tracks total # cells in each state

CellState.Position = []; % position of each cell in the system
CellState.TimesMoved = [];

CellState.AvgEcad  = [];
CellState.AvgNcad  = [];
CellState.AvgSnail = [];
CellState.AvgZeb1  = [];

end

