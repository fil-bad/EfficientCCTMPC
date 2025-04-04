% Parameters
n = 3;   % Dimension
m = 10;  % Controls the number of points

% Angle discretization
theta = cell(n-1, 1);
for k = 1:n-2
    theta{k} = linspace(0, pi, m);
end
theta{n-1} = linspace(0, 2*pi, m);

% Generate grid of angles
[Theta{1:n-1}] = ndgrid(theta{:});

% Flatten grids
for k = 1:n-1
    Theta{k} = Theta{k}(:);
end

% Number of points
numPoints = numel(Theta{1});

% Initialize vertices matrix
vertices = zeros(numPoints, n);

% Calculate Cartesian coordinates
for idx = 1:numPoints
    sinProd = 1;
    for dim = 1:n
        if dim == n
            angle = Theta{n-1}(idx);
            vertices(idx, dim) = sinProd * sin(angle);
        elseif dim == n-1
            angle = Theta{n-1}(idx);
            vertices(idx, dim) = sinProd * cos(angle);
        else
            angle = Theta{dim}(idx);
            vertices(idx, dim) = sinProd * cos(angle);
            sinProd = sinProd * sin(angle);
        end
    end
end

% Remove duplicate points caused by angle discretization
vertices = unique(round(vertices * 1e6) / 1e6, 'rows');

% Compute convex hull facets
facets = convhulln(vertices);

% Visualization for n = 3
if n == 3
    figure;
    trisurf(facets, vertices(:,1), vertices(:,2), vertices(:,3), ...
        'FaceColor', [0.8, 0.9, 1.0], 'EdgeColor', 'k');
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Tessellation of the Unit Sphere with m = ', num2str(m)]);
    grid on;
end

