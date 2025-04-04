% see https://www.jstor.org/stable/27594324 for the MM algorithm.
% see also https://arxiv.org/pdf/1902.00220 and later work 
% https://arxiv.org/pdf/2105.00401 for the Predefined Evenly-Distributed 
% Class Centroids (PEDCC) 
function points = numericalSpherePoints(nx, pts,randSeed, max_iter, tol, plt)
    % n: number of points
    % d: dimension of the hypersphere (d-1)
    % max_iter: maximum number of iterations
    % tol: tolerance for convergence

    % NOTE: we need at least a simplex to guarantee to have a non-empty
    % interior!
    str = sprintf("We need at least a simplex, that is %d vertices in R^%d",nx+1,nx);
    assert(pts >= nx+1, str);

    if nargin < 3
        rng("default");
    else
        rng(randSeed);
    end
    if nargin < 4
        max_iter = 1e6;
    end
    if nargin < 5
        tol = 1e-6;
    end
    if nargin < 6
        plt = false;
    end

    % Step 1: Initialize points randomly on the hypersphere
    points = initialize_points_on_hypersphere(pts, nx);
    
    for i = 1:max_iter
        old_points = points;
        
        % Step 2: Minorization step - construct surrogate function
        surrogate_function = construct_surrogate_function(points);
        
        % Step 3: Maximization step - optimize surrogate function
        points = optimize_surrogate_function(surrogate_function);
        
        % Use squared Frobenius norm to avoid sqrt overhead
        if sum((points - old_points).^2, 'all') < tol^2
            if plt
            disp("Iterations required: " + i)
            disp("Squared Frobenius Norm: " + sum((points - old_points).^2, 'all'))    
            end
            break;
        end
    end
    
    % Plotting the points and edges for 2D and 3D cases
    if not(plt)
        return;
    end
    if nx == 2
        plot_2d(points);
    elseif nx == 3
        plot_3d(points);
    else
        disp('Plotting is only available for 2D and 3D cases.');
    end
end

function points = initialize_points_on_hypersphere(n, d)
    % Generate n random points on a d-dimensional hypersphere
    % points = zeros(1,d); points(1) = 1;
    % points = [points; randn(n-1, d)];
    points = randn(n, d);
    points = points ./ sqrt(sum(points.^2, 2));
end

function surrogate_function = construct_surrogate_function(points)
    % Vectorized computation of the surrogate function
    n = size(points, 1);
    diff = reshape(points, [n, 1, size(points,2)]) - reshape(points, [1, n, size(points,2)]);
    distances = sqrt(sum(diff.^2, 3));
    distances(1:n+1:end) = Inf;  % Avoid division by zero on diagonal
    normDiff = diff ./ distances;
    surrogate_function = squeeze(sum(normDiff, 2));
end

function points = optimize_surrogate_function(surrogate_function)
    % Optimize the surrogate function to update the points
    points = surrogate_function ./ sqrt(sum(surrogate_function.^2, 2));
end

function plot_2d(points)
    % Plot the points on a 2D circle
    theta = linspace(0, 2*pi, 100);
    x = cos(theta);
    y = sin(theta);
    
    figure;
    plot(x, y, 'k--'); % Plot the circle
    hold on;
    plot(points(:, 1), points(:, 2), 'bo');
    axis equal;
    title('Equidistant Points on a 2D Circle');
    xlabel('x');
    ylabel('y');
    
    % Plot external edges between points
    k = convhull(points(:, 1), points(:, 2));
    for i = 1:length(k)-1
        plot([points(k(i), 1), points(k(i+1), 1)], [points(k(i), 2), points(k(i+1), 2)], 'k-');
    end
    hold off;
end

function plot_3d(points)
    % Plot the points on a 3D sphere
    [X, Y, Z] = sphere(100);
    
    figure;
    surf(X, Y, Z, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % Plot the sphere
    hold on;
    scatter3(points(:, 1), points(:, 2), points(:, 3), 'filled');
    axis equal;
    title('Equidistant Points on a 3D Sphere');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    % Plot external edges between points
    k = convhull(points);
    trisurf(k, points(:, 1), points(:, 2), points(:, 3), 'FaceColor', 'none', 'EdgeColor', 'k');
    hold off;
end
