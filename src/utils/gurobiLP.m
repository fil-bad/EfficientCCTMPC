function [x, fval, exitflag, output] = gurobiLP(f, Aineq, bineq, Aeq, beq, modelsense)
    % GUROBI_LP Solves an LP problem using Gurobi
    %   [x, fval, exitflag, output] = gurobi_lp(f, Aineq, bineq, Aeq, beq, modelsense)
    %   solves the LP problem:
    %       minimize    f'*x
    %       subject to  Aineq*x <= bineq
    %                   Aeq*x == beq
    
    if (nargin == 3); Aeq = []; beq = []; end
    if nargin < 6 || isempty(modelsense)
        modelsense = 'min';
    end
    
    % Create a Gurobi model
    model.A = sparse([Aineq; Aeq]); % Combine inequality and equality constraints
    model.obj = f;
    model.modelsense = modelsense; % Use provided modelsense ('min' or 'max')
    model.rhs = [bineq; beq]; % Combine right-hand sides
    model.sense = [repmat('<', size(Aineq, 1), 1); repmat('=', size(Aeq, 1), 1)]; % Inequality and equality constraints
    
    % Define variable types (continuous in this case)
    model.vtype = 'C';
    
    % Set Gurobi parameters to suppress verbose output 
    params.outputflag = 0; 
    
    % Solve the model 
    result = gurobi(model, params);
    
    % Extract the results
    x = result.x;
    fval = result.objval;
    
    % Interpret Gurobi status codes
    if strcmp(result.status, 'OPTIMAL')
        exitflag = 1; % Solution found
    else
        exitflag = 0; % No solution found or other issues
    end
    
    % Output additional information
    output = result;
    
%             % Display the results (for debugging purposes, can be removed in production)
%             disp('Optimal Solution:');
%             disp(x);
%             disp('Objective Value:');
%             disp(fval);
end