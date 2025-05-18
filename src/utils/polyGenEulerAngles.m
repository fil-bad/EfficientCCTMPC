function [vertices] = polyGenEulerAngles(nx, pts, plt)
    % nx: dimensionality
    % pts: #points of discretization for each angle (suggested: 100 nx=2, 10 nx=3)
    if nargin < 3
        plt = false;
    end

    % Angle discretization
    theta = cell(nx-1, 1);
    for k = 1:nx-2
        theta{k} = linspace(0, pi, pts);
    end
    % Adjust last angle to avoid duplicate between 0 and 2*pi
    temp = linspace(0, 2*pi, pts+1);
    theta{nx-1} = temp(1:end-1);

    % Generate grid of angles
    [Theta{1:nx-1}] = ndgrid(theta{:});
    for k = 1:nx-1
        Theta{k} = Theta{k}(:);
    end

    % Number of points
    numPoints = numel(Theta{1});

    % Vectorized Cartesian conversion
    % Assemble angles into a matrix: each row corresponds to one point, columns to angles.
    ThetaMat = cell2mat(Theta'); 
    ThetaMat = reshape(ThetaMat, numPoints, nx-1);
    
    vertices = zeros(numPoints, nx);
    % Standard spherical coordinate conversion:
    vertices(:,1) = cos(ThetaMat(:,1));
    sinProd = sin(ThetaMat(:,1));
    for dim = 2:nx-1
        vertices(:,dim) = sinProd .* cos(ThetaMat(:,dim));
        sinProd = sinProd .* sin(ThetaMat(:,dim));
    end
    vertices(:,nx) = sinProd;
    
    % Remove duplicate points caused by angle discretization
    vertices = unique(round(vertices * 1e6) / 1e6, 'rows');

    % Visualization (activated only if plt is true)
    if plt
        if nx == 3
            facets = convhulln(vertices);
            figure;
            trisurf(facets, vertices(:,1), vertices(:,2), vertices(:,3), ...
                'FaceColor', [0.8, 0.9, 1.0], 'EdgeColor', 'k');
            axis equal;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title(['Tessellation of the Unit Sphere with m = ', num2str(pts)]);
            grid on;
        elseif nx == 2
            k = convhull(vertices(:,1), vertices(:,2));
            figure;
            plot(vertices(:,1), vertices(:,2), 'bo');
            hold on;
            plot(vertices(k,1), vertices(k,2), 'r-');
            axis equal;
            xlabel('X'); ylabel('Y');
            title(['Tessellation of the Unit Circle with m = ', num2str(pts)]);
            grid on;
            hold off;
        end
    end
end

