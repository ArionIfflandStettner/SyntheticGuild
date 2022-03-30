% Function to run FBA using Gurobi
% c - objective vector
% S - stoichiometric matrix
% b - vector of changes in metabolite concentrations (typically zeros)
% contypes - types of constraints for Gurobi
% lb, ub - vectors of lower and upper bounds
% constr_changes - matrix of changes to the mass balance constraints b

function [v, val, v_min, val2] = FBA_Gurobi(c,S,b,contypes,lb,ub,constr_changes)
% Set output options
opts.FeasibilityTol = 1e-6;
opts.OutputFlag = 0;
opts.DisplayInterval = 1;

% Set variable types
vartypes = char( 'C'*ones( 1,size(S,2) ) )';

% Check if there are any changes to metabolite concentrations
if ~isempty(constr_changes)
    % First column is index of mass balance constraint, second column is
    % type of constraint change, third column is value of constraint
    for i = 1:numel(constr_changes.index)
        contypes( constr_changes.index(i) ) = constr_changes.type(i);
        b( constr_changes.index(i) ) = constr_changes.eps(i);
    end
end

% Perform maximization
model = struct('A',sparse(S),'obj',c,'sense',contypes','rhs',b,'lb',lb,'ub',ub,'vtype',vartypes,'modelsense','max');
params = struct('outputflag',opts.OutputFlag,'FeasibilityTol',opts.FeasibilityTol,'DisplayInterval',opts.DisplayInterval);
result = gurobi(model,params);
if ~strcmp(result.status,'OPTIMAL')
    fprintf('Error: Nonoptimal solution in initial optimization!\n')
    v = [];
    val = [];
    v_min = [];
    val2 = [];
else
    v = result.x;
    val = result.objval;
    
    % Minimize the sum of the absolute values of fluxes
    S = [S -S];
    
    % Add an additional constraint for the flux just maximized
    S = [ S; [c -c] ];
    b = [b';c*v];
    contypes = vertcat(contypes,'=');
    
    % Add extra variables
    vartypes = char( 'C'*ones( 1,size(S,2) ) )';
    
    % Lower bound on forward fluxes is max(zero, LB), lower bound on reverse fluxes is max(zero, -UB)
    Lbound = [max(lb, zeros(numel(lb),1)); max(-ub, zeros(numel(lb),1))];
    
    % Upper bound on forward fluxes is max(zero, UB), upper bound on reverse fluxes is max(zero, -LB)
    Ubound = [max(ub, zeros(numel(lb),1)); max(-lb, zeros(numel(lb),1))];
    
    % Revise objective
    new_obj = ones( size( S,2 ), 1 );
    
    % Rerun simulation to minimize sum of fluxes
    model = struct('A',sparse(S),'obj',new_obj,'sense',contypes','rhs',b,'lb',Lbound,'ub',Ubound,'vtype',vartypes,'modelsense','min');
    params = struct('outputflag',opts.OutputFlag,'FeasibilityTol',opts.FeasibilityTol,'DisplayInterval',opts.DisplayInterval);
    result = gurobi(model,params);
    if ~strcmp(result.status,'OPTIMAL')
        % Solution not optimal, print error
        fprintf('Error: Nonoptimal solution when minimizing absolute values!\n')
        v_min = [];
        val2 = [];
    else
        v2 = result.x;
        val2 = result.objval;
        
        % Merge positive and negative values
        v_min = v2(1:end/2) - v2(end/2+1:end);
        
        % Check if objectives match
        if c*(v_min-v) > 1e-5
            fprintf('Error: Minimization of absolute values of fluxes failed.\n')
        end
    end
end