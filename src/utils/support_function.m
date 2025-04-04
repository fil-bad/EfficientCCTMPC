function h = support_function(A, b, D)
% SUPPORT_FUNCTION Computes the support function for a polytope in H-representation.
%
%   h = SUPPORT_FUNCTION(A, b, D)
%
%   Computes the support function for the polytopic set
%
%       P = { x in R^n : A * x <= b }
%
%   in the direction(s) specified by D. Each row of D is interpreted as a
%   direction vector. The support function in a direction d is defined as:
%
%       h_P(d) = max { d' * x : A*x <= b }
%
%   Input:
%       A - m x n matrix, where each row defines a half-space.
%       b - m x 1 vector, the right-hand side of the inequality.
%       D - k x n matrix, where each row is a direction vector.
%
%   Output:
%       h - k x 1 vector of support function values. For each direction
%           D(i,:), h(i) is the maximum value of d'*x over x in P.
%
%   Example:
%       % Define a unit box in R^2:
%       A = [ -1  0;
%              0 -1;
%              1  0;
%              0  1 ];
%       b = [0; 0; 1; 1];
%
%       % Define two direction vectors:
%       D = [1 1;
%            1 -1];
%
%       % Compute support function values:
%       h = support_function(A, b, D)
%
%   Notes:
%       - This function uses MATLAB's 'linprog' to solve the linear
%         program. If the LP is found unbounded in a given direction, the
%         support function returns Inf for that direction.
%
%   Author: [Your Name]
%   Date: [Today's Date]

    % Determine the number of variables (n) and number of directions (k)
    n = size(A, 2);
    k = size(D, 1);
    
    % Preallocate output vector h
    h = zeros(k, 1);
    
    % Check if Gurobi is available; if not, fall back to linprog.
    isGurobi = (exist('gurobi', 'file') == 2);
    if ~isGurobi
        options = optimoptions('linprog', 'Display', 'none');
    end
    
    if isGurobi
        if k > 1
            % --- Begin batch optimization with Gurobi for speed ---
            % Construct block diagonal constraint matrix for k LPs
            % (Each block is A)
            rep_A = repmat({A}, 1, k);
            A_batch = blkdiag(rep_A{:});
            % Stack each direction vector (as column) into one objective vector
            obj_batch = reshape(D', [], 1);
            % Replicate b for each LP block and set constraint sense
            rhs_batch = repmat(b, k, 1);
            sense_batch = repmat('<', length(b)*k, 1);
            
            model.A = sparse(A_batch);
            model.obj = obj_batch;
            model.rhs = rhs_batch;
            model.sense = sense_batch;
            model.modelsense = 'max';
            
            try
                result = gurobi(model);
            catch gurobiError
                error('Gurobi Error: %s', gurobiError.message);
            end
            
            if strcmp(result.status, 'OPTIMAL')
                % Partition the solution vector and compute h(i) = d_i' * x_i
                x_sol = result.x;
                for i = 1:k
                    xi = x_sol((i-1)*n+1:i*n);
                    h(i) = D(i, :) * xi;
                end
            else
                error('Batch Gurobi optimization did not converge. Status: %s', result.status);
            end
            % --- End batch optimization ---
        else
            % Single direction: use the original per-iteration Gurobi call.
            d = D(1, :)';
            model.A = sparse(A);
            model.obj = d;                  % maximize d' * x directly
            model.rhs = b;
            model.sense = repmat('<', size(b));
            model.modelsense = 'max';
            
            try
                result = gurobi(model);
            catch gurobiError
                error('Gurobi Error: %s', gurobiError.message);
            end
            
            if strcmp(result.status, 'OPTIMAL')
                h(1) = result.objval;
            elseif strcmp(result.status, 'UNBOUNDED') || strcmp(result.status, 'INF_OR_UNBD')
                h(1) = Inf;
            else
                error('Gurobi did not converge for direction %d. Status: %s', 1, result.status);
            end
        end
    else
        % Use linprog (processing one direction per iteration)
        for i = 1:k
            d = D(i, :)';
            f = -d;
            [~, fval, exitflag] = linprog(f, A, b, [], [], [], [], options);
            if exitflag == 1
                h(i) = -fval;
            elseif exitflag == -3
                h(i) = Inf;
            else
                error('LP did not converge for direction %d. Exit flag: %d', i, exitflag);
            end
        end
    end
end
