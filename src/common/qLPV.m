classdef qLPV < handle
    % The system is expressed as x+ = Ax+Bu+w, where (A,B) \in \Delta,
    % convex hull of n_m different models.

    properties (SetAccess = private)
        A_convh % Convex hull of state matrices
        B_convh % Convex hull of input matrices
        Bw      % Disturbance to state matrix
        W_dist  % Disturbance set
        X       % State constraints 
        U       % Input constraints
        
    end

    properties (Access = private)
        A_curr_han, B_curr_han  % function handle to get current system realization
    end

    properties(Dependent)
        nx, nu, ny  % state input and output dimensions
        nm          % #vertices of convex hull (models)
    end

    methods % GETTER methods
        function nx = get.nx(obj)
            nx = size(obj.A_convh{1},1); 
        end
        function nu = get.nu(obj)
            nu = size(obj.B_convh{1},2);
        end
        function nm = get.nm(obj)
            nm = numel(obj.A_convh);
        end
    end

    methods (Access = public)
        function obj = qLPV(A_convh, B_convh, Bw, X, U, W_dist, A_curr_han, B_curr_han)

            assert(numel(A_convh)==numel(B_convh),"Define each model as (A_i,B_i) pair.")
            obj.A_convh = A_convh;
            obj.B_convh = B_convh;
            obj.Bw = Bw;
            obj.W_dist = W_dist;
            obj.X = X;
            obj.U = U;

            if nargin == 6 % Pick one model among all (works also for LTI)
                obj.A_curr_han = @(A_hull,~,~) A_hull{randi([1,numel(A_hull)])};
                obj.B_curr_han = @(B_hull,~,~) B_hull{randi([1,numel(B_hull)])};
            elseif nargin == 8
                obj.A_curr_han = A_curr_han; %{A_convh, x_curr, ext_params}
                obj.B_curr_han = B_curr_han; %{B_convh, u_curr, ext_params}
            else
                error("Define both A_curr and B_curr function handles.")
            end
        end

        function updateSysMatrices(obj, A_convh_new, B_convh_new)
            % used for online model refinement
            obj.A_convh = A_convh_new;
            obj.B_convh = B_convh_new;
        end

        function x_next = step_noDist(obj, x, u, varargin)
            % varargin handles external parameter vector for an LPV system
            x_next = obj.A_curr(x, varargin{:})*x + obj.B_curr(x,varargin{:})*u;
        end

        function x_next = step(obj, x, u, varargin)
            % varargin handles external parameter vector for an LPV system
            x_next = obj.step_noDist(x,u,varargin{:}) + obj.Bw*obj.W_dist.randomPoint();
        end

        function A_curr = A_curr(obj,varargin)
            % varargin: {x_curr, ext_params}
            A_curr = obj.A_curr_han(obj.A_convh, varargin{:});
        end

        function B_curr = B_curr(obj,varargin)
            % varargin: {u_curr, ext_params}
            B_curr = obj.B_curr_han(obj.B_convh, varargin{:});
        end

    end

    % methods(Abstract) % to be implemented in a concrete subclass
    %     A_curr(obj,varargin) % varargin: {state, ext_param_vector}
    %     B_curr(obj,varargin) % varargin: {state, ext_param_vector}
    % end

end