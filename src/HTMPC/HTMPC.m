classdef HTMPC < handle
    % CCTMPC implementation. Based on OptimalController class. For clarity
    % of presentation, the QP used for computing the closed-loop control
    % law is included in this class instead of solving a single QP. This
    % class acts also as a wrapper for the OCP solver.

    properties (SetAccess = private)
        centralInput % to identify chosen input parameterization

        ocpSol % cell array, solution from the main QP problem
        lambdaSol % solution of convex combination of the vertex control law

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
            % y = Fz + \alpha*ys
            y0Sol = obj.OCP.ccPoly.F*obj.ocpSol{1}(:,1) + obj.ocpSol{2}(1)*obj.rciSol{1};
        end
    end


    methods (Access = public)
        function obj = HTMPC(sys, ccPoly, CFM, N, gamma, sdp_opts, centralInput)

            if centralInput
                obj.OCP = H_OCP_cin(sys, ccPoly, CFM, N, gamma, sdp_opts, false);
            else
                obj.OCP = H_OCP(sys, ccPoly, CFM, N, gamma, sdp_opts,false);
            end
            obj.centralInput = centralInput;
            obj.mpcLaw = obj.initMPCLaw();

            obj.ocpSol = cell(4,1); % {z_t, alpha_t, u_t, costOCP_t} fields
        end


        function uMPC = solve(obj, x_curr, varargin)
            % solve OCP (in parameter space)
            obj.ocpSol = obj.OCP.solve_OCP(x_curr);

            % compute the optimal input (in state space) and lambda multipliers
            % mpcSol = obj.mpcLaw(x_curr, obj.ocpSol{1}(:,1), obj.ocpSol{2}(1), obj.ocpSol{3}(:,:,1));
            mpcSol = obj.mpcLaw(x_curr, obj.y0_sol, obj.ocpSol{3}(:,:,1));
            [uMPC, obj.lambdaSol] = mpcSol{:};
            if isnan(uMPC)
                if sum(abs(obj.findFeasibleX0(x_curr,eye(obj.OCP.sys.nx)) - x_curr)) > 1e-5
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
                case 1
                    feasX0 = obj.OCP.findFeasibleX0(x_curr,varargin{:});
                otherwise
                    error("%d arguments provided (max 1 expected).",length(varargin));
            end
        end
    end

    methods (Access = private)

        function mpcLaw = initMPCLaw(obj)
            % define the YALMIP optimizer for the control-law computation

            f = obj.OCP.ccPoly.f;
            v = obj.OCP.ccPoly.v;
            nu = obj.OCP.sys.nu;
            nx = obj.OCP.sys.nx;
            Vi_s = obj.OCP.ccPoly.Vi_s;

            % z0Opt = sdpvar(nx, 1,'full');
            % alpha0Opt = sdpvar(1);
            y0Opt = sdpvar(f, 1,'full');
            u0Opt = sdpvar(nu, v,'full');
            x0 = sdpvar(nx, 1,'full');

            lambda = sdpvar(v, 1,'full');

            constr = [];
            constr = [constr; sum(lambda)==1; lambda >= 0];
            x_eq_sum = 0;
            for j=1:v
                x_eq_sum = x_eq_sum + lambda(j)*Vi_s{j}*y0Opt;
                % x_eq_sum = x_eq_sum + lambda(j)*(z0Opt + alpha0Opt*Vi_s{j}*obj.OCP.ys);
            end
            constr = [constr; x0 == x_eq_sum];

            cost = weighted2NormSquared(lambda, eye(size(lambda,1))); % squared 2-norm
            % cost = nnz(lambda); % 0-norm

            mpcLaw = optimizer(constr, cost, ...
                obj.OCP.sdp_opts, ...
                {x0, y0Opt, u0Opt}, ... % input parameters
                {u0Opt*lambda, lambda} ... % output parameters
                ); % (nu x v) * (v x 1) -> (nu x 1)

            % mpcLaw = optimizer(constr, cost, ...
            %    sdpsettings(obj.OCP.sdp_opts{:}), ...
            %    {x0, z0Opt, alpha0Opt, u0Opt}, ... % input parameters
            %    {u0Opt*lambda, lambda} ... % output parameters
            %    ); % (nu x v) * (v x 1) -> (nu x 1)
        end

    end

end