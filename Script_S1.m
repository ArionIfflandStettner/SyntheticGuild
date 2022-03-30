%% Script for running FBA with metabolic models

clear
close all
addpath('/Library/gurobi811/mac64/matlab') % Gurobi path
load vhep % Load model (rpom, thala, or vhep)

% Run FBA
bmIndex = find(strcmp('BiomassRxn',m.rxns));
cIndex = find(strcmp('EX_glc',m.rxns)); % Glucose as carbon source
m.lb(cIndex) = -10;
b = zeros(1,size(m.S,1));
c = zeros(1,size(m.S,2));
c(bmIndex) = 1;
contypes = char( '='*ones( 1,size(m.S,1) ) )';
constr_changes = [];
[v, val, v2, val2] = FBA_Gurobi(c,m.S,b,contypes,m.lb,m.ub,constr_changes);
