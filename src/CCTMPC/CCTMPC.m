classdef CCTMPC < handle
    % CCTMPC implementation. Based on OptimalController class. For clarity
    % of presentation, the QP used for computing the closed-loop control
    % law is included in this class instead of solving a single QP. This
    % class acts also as a wrapper for the OCP solver.

    properties (SetAccess = private)
        var_convh % flag, enable variability for system convex hulls

        ocpSol % cell array, solution from the main QP problem
        lambdaSol % solution of convex combination of the vertex control law

        hausDist % Hausdorff distance from state constraint ||v_x_i -z0||

    end

    properties (Access = private)
        OCP % YALMIP optimizer, Optimal Controller
        mpcLaw % YALMIP optimizer, compute the input in state (and not parameter) space
    end

    properties(Dependent)
        Lyapunov_cost % as defined in [(ii),Corollary 1] from paper
        feasRegion % Polytope, (approximated) feasible region of the underlying OCP

        rciSol % cell array, solution from RCI

        y0_sol % compute the first step closed-loop parameter vector
    end

    methods % GETTER methods
        function Lyap_cost = get.Lyapunov_cost(obj)
            Lyap_cost = obj.ocpSol{end};% - obj.rciSol{end};
        end

        function feasReg = get.feasRegion(obj)
            feasReg = obj.OCP.feasRegion;
        end

        function rciSol = get.rciSol(obj)
            rciSol = {obj.OCP.ys, obj.OCP.us, obj.OCP.RCI_setcost};
        end

        function y0Sol = get.y0_sol(obj)
            y0Sol = obj.ocpSol{1}(:,1);
        end
    end


    methods (Access = public)
        function obj = CCTMPC(sys, ccPoly, CFM, N, gamma, sdp_opts, OCP_opts)
            if nargin == 6; OCP_opts = struct("useBasicOCP",false, "var_convh",false); end
            if not(isfield(OCP_opts,"useBasicOCP")) % useBasicOCP not defined
                OCP_opts.useBasicOCP = false;
                if not(isfield(OCP_opts,"var_convh")) % fixed convh of models
                    OCP_opts.var_convh = false;
                end
            end

            if OCP_opts.useBasicOCP
                obj.OCP = CC_OCP_basic(sys, ccPoly, CFM, N, gamma, sdp_opts);
            else
                obj.OCP = CC_OCP(sys, ccPoly, CFM, N, gamma, sdp_opts, OCP_opts.var_convh);
            end
            obj.mpcLaw = obj.initMPCLaw();

            obj.var_convh = OCP_opts.var_convh;

            obj.ocpSol = cell(3,1); % {y_t,u_t,costOCP_t} fields

            obj.computeHausdorffDist(sys);
        end


        function uMPC = solve(obj, x_curr, varargin)
            % solve OCP (in parameter space)
            switch class(obj.OCP)
                case 'CC_OCP'
                    obj.solve_fast(x_curr, varargin{:});
                case 'CC_OCP_basic'
                    obj.solve_basic(x_curr, varargin{:});
                otherwise
                    error("Unsupported OCP class: %s", class(obj.OCP))
            end
            % compute the optimal input (in state space) and lambda multipliers
            mpcSol = obj.mpcLaw(x_curr, obj.ocpSol{1}(:,1), obj.ocpSol{2}(:,:,1));
            [uMPC, obj.lambdaSol] = mpcSol{:};
            if isnan(uMPC)
                feas_x = obj.findFeasibleX0(x_curr,eye(obj.OCP.sys.nx));
                if sum(abs(feas_x{1}-x_curr) > 1e-5)
                    error("The current state does not lie inside the feasible"+...
                        " region of the control scheme. Consider changing initial condition.")
                end
                error("uMPC returned NaN. Consider enabling debug mode.")
            end
        end


        function feasRegion = computeFeasRegion(obj, n_facets, varargin)
            % interface method to the underlying OCP class
            feasRegion = obj.OCP.computeFeasRegion(n_facets, varargin{:});
        end


        function feasX0 = findFeasibleX0(obj,x_curr,varargin)
            % interface method to the underlying OCP class
            switch length(varargin)
                case 0
                    feasX0 = obj.OCP.findFeasibleX0(x_curr,eye(obj.OCP.sys.nx));
                    feasX0 = feasX0{1};
                case 1
                    feasX0 = obj.OCP.findFeasibleX0(x_curr,varargin{:});
                    feasX0 = feasX0{1};
                otherwise
                    error("%d arguments provided (max 1 expected).",length(varargin));
            end
        end
    end

    methods (Access = private)

        %% NOTE: by having removed the RCI computation, we cannot update
        % this set to take advantages when matrices are updated!

        function solve_fast(obj, x_curr, varargin)
            % two possible varargin cases:
            % length(varargin)==0: for fixed systems, RCI computed only if
            %       current reference is different from previous one
            % length(varargin)==2: system dynamics (A_convh,B_convh) are
            %       updated.
            switch length(varargin)
                case 0
                    % solve OCP (solution in parameter space)
                    if obj.var_convh
                        obj.ocpSol = obj.OCP.solve_OCP(x_curr, obj.OCP.sys.A_convh, obj.OCP.sys.B_convh);
                    else
                        obj.ocpSol = obj.OCP.solve_OCP(x_curr);
                    end
                case 2 % (A_convh,B_convh)
                    assert(obj.var_convh == true,...
                        "Dynamic system matrices have been initialized as fixed, cannot be modified online!");
                    obj.OCP.sys.updateSysMatrices(varargin{:});

                    % solve OCP (solution in parameter space)
                    obj.ocpSol = obj.OCP.solve_OCP(x_curr, varargin{:});

                otherwise % nargin-1 as obj is counted as argument
                    error("Wrong number of arguments: %d. Only 2 and 4 supported.", nargin-1);
            end
        end


        function solve_basic(obj, x_curr, varargin)
            % two possible varargin cases:
            % length(varargin)==0: for fixed systems, RCI computed only if
            %       current reference is different from previous one
            % length(varargin)==2: system dynamics (A_convh,B_convh) are
            %       updated; RCI always computed.

            % nargin-1 as obj is counted as argument
            assert(any(length(varargin)==[0,2]),...
                "Wrong number of arguments: %d. Only 2 and 4 supported.", nargin-1)

            % update system matrices
            if length(varargin) == 2
                obj.OCP.sys.updateSysMatrices(varargin{:});
            end

            % solve OCP (solution in parameter space)
            obj.ocpSol = obj.OCP.solve_OCP(x_curr);
        end


        function mpcLaw = initMPCLaw(obj)
            % define the YALMIP optimizer for the control-law computation
            weightSq2Norm = @(vector, mat) vector'*mat*vector;

            f = obj.OCP.ccPoly.f;
            v = obj.OCP.ccPoly.v;
            nu = obj.OCP.sys.nu;
            nx = obj.OCP.sys.nx;
            Vi_s = obj.OCP.ccPoly.Vi_s;

            y0Opt = sdpvar(f, 1,'full');
            u0Opt = sdpvar(nu, v,'full');
            x0 = sdpvar(nx, 1,'full');

            lambda = sdpvar(v, 1,'full');

            constr = [];
            constr = [constr; sum(lambda)==1; lambda >= 0];
            x_eq_sum = 0;
            for j=1:v
                x_eq_sum = x_eq_sum + lambda(j)*Vi_s{j}*y0Opt;
            end
            constr = [constr; x0 == x_eq_sum];

            cost = weightSq2Norm(lambda, eye(size(lambda,1))); % squared 2-norm
            % cost = nnz(lambda); % 0-norm

            mpcLaw = optimizer(constr, cost, ...
                obj.OCP.sdp_opts, ...
                {x0, y0Opt, u0Opt}, ... % input parameters
                {u0Opt*lambda, lambda} ... % output parameters
                ); % (nu x v) * (v x 1) -> (nu x 1)
        end

        function computeHausdorffDist(obj,sys)
            x_verts = sys.X.V';
            d = 0;
            for v = 1:size(x_verts,2)
                feasX0 = obj.OCP.findFeasibleX0(x_verts(:,v), eye(sys.nx));
                d = max(d, sqrt(feasX0{2})); % cost
            end
            obj.hausDist = d;
        end
    end

end