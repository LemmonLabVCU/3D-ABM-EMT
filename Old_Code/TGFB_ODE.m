% FUNCTION: TGFB_ODE.m
% Goal: Calculates TGFB signaling dynamics using the euler method.

% [TGFB conc, EMT marker conc] = TGFB_ODE(i, t, tau, C0, T)

% INPUTS: 
% _______________________________________________________________________
%   i:  time (t)

%   t:  total time interval (dt:dt:Tfinal)

%   tau: time step (dt)

%   C0: Cell-specific concentrations of EMT markers defined in
%   CellState.conc
        % [snail, SNAIL, miRNA34, zeb1, Zeb1, miRNA200, ECAD, NCAD]
 
%   T: TGFB concentration defined at point [i,j,k].        

% Other Variables used: 

%   k:  TGFB1 signaling kinetics/ reaction rates (defined as 1/hr or um/hr)
%   k = [k0T, kT, kdT, k0s, ks, kds, kS, kdS, 
%       k03, k3, kd3, k0z, kz, kdz, kZ, kdZ, 
%       k02, k2, kd2, ke1, ke2, kde, kn1, kn2,
%        kdn, kbT]

%   J and n: Hill Function parameters

% OUTPUTS: 
% _________________________________________________________________________
%   T: New TGFB concentration at position (i,j,k) based on ODE outputs

%   C: New Cell EMT markers concentrations based on ODE outputs
% ________________________________________________________________________

% NOTES:
% Functions and parameters used here are based on the Tian model: 
%   Source: Tian, X.-J., Zhang, H., & Xing, J. (2013). 
%           Coupled Reversible and Irreversible Bistable Switches 
%           Underlying TGFβ-induced Epithelial to Mesenchymal Transition. 
%           Biophysical Journal, 105(4). doi:10.1016/j.bpj.2013.07.011

% Parameters modified from initial paper: TGFB, k1, k2, k3
% _________________________________________________________________________

function [T, C] = TGFB_ODE(i, t, tau, C0, T)
 
% Hill Function values:
J = [0.06, 1.6, 0.08, 0.15, 0.36, 3.5, 0.06, 5, 0.2, 0.2, 0.5, 0.2, 0.5]; 
    % (uM) [JT, Js, JS, J13, J23, Jz, JZ, J12, J22, J1e, J2e, J1n, J2n] 
n = [2, 2, 2, 2, 2]; % [nnt, nns, nnz, nnr2, nnr3]

% Reaction Rates 
k = [0.0659, 1.425, 0.086, 0.0062, 0.03, 0.09, 17, 1.66, ...
    0.0012, 0.012, 0.035, 0.003, 0.06, 0.09, 17, 1.66,...
    0.0002, 0.012, 0.035, 1, 0.6, 0.5, 1, 0.6, ...
    0.5]; 

%   k = [k0T, kT, kdT, k0s, ks, kds, kS, kdS, 
%       k03, k3, kd3, k0z, kz, kdz, kZ, kdZ, 
%       k02, k2, kd2, ke1, ke2, kde, kn1, kn2,
%        kdn, kbT]

h = tau; % time step (hr)

x = [T, C0]; % initial conditions

%Model equations
sT = @(t,x)  k(1) + (k(2) / (1 + (x(7) / J(1))^n(1)))  - k(3)*(x(1)); %soluble TGFB

s = @(t,x) k(4) + k(5)*(((x(1) / J(2)).^n(2) / (1 + (x(1) / J(2)).^n(2)))) - k(6)*x(2); % Snail Transcription factor (s)
 
S =  @(t,x) ( (k(7)*x(2)) / ( 1 + ( x(4)/ J(3)).^n(5)) ) - k(8)*x(3); % Snail Protein (S)

R3 = @(t,x) k(9) + k(10)*(1 / (1 + (x(3)/J(4)).^n(2) + (x(6)/J(5)).^n(3))) - k(11)*x(4); % miRNA-34 (R3)

z = @(t,x) k(12) + k(13)*(( (x(3)/ J(6)).^n(2) ) / ( 1 + (x(3)/ J(6)).^n(2) )) - k(14)*x(5); % Zeb1 transcription factor (z)

Z = @(t,x) k(15)*x(5)*(1 / ( 1 + (x(7) / J(7)).^n(4) )) - k(16)*x(6); % Zeb1 Protein (Z)

R2 = @(t,x) k(17) + k(18)*(1 /(1 + (x(3)/J(8)).^n(2) + (x(6)/J(9)).^n(3))) - k(19)*x(7); % miRNA-200 (R2)

Ecad = @(t,x) k(20)*(1 / ( (x(3)/J(10)).^n(2) + 1 )) + k(21)*(1 / ( (x(6)/J(11)).^n(3) + 1 )) - k(22)*x(8); % E-cadherin (E)

Ncad = @(t,x) k(23)*( ((x(3) / J(11)).^n(2)) / (((x(3) / J(11)).^n(2)) + 1)) + k(24)*( ((x(6) / J(12)).^n(3)) / (((x(6) / J(12)).^n(3)) + 1)) - k(25)*x(9); % N-cadherin (N)

% Euler Method
    x(1) = x(1) + h*sT(t(i), x(:)); %last y value + slope
    x(2) = x(2) + h*s(t(i), x(:)) ; 
    x(3) = x(3) + h*S(t(i), x(:)) ; 
    x(4) = x(4) + h*R3(t(i), x(:)) ; 
    x(5) = x(5) + h*z(t(i), x(:)) ; 
    x(6) = x(6) + h*Z(t(i), x(:)) ; 
    x(7) = x(7) + h*R2(t(i), x(:)) ; 
    x(8) = x(8) + h*Ecad(t(i), x(:)) ;
    x(9) = x(9) + h*Ncad(t(i), x(:)) ;
    
C = x(2:9);
T = x(1);

end
