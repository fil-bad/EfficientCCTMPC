function normSquared = weighted2NormSquared(x, W)
    if nargin < 2
        % No weighting matrix provided; compute standard squared 2-norm
        normSquared = x' * x;
    else
        % Weighting matrix provided; compute weighted squared 2-norm
        normSquared = x' * W * x;
    end
end
