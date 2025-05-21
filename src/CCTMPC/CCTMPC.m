classdef CCTMPC < handle
    % CCTMPC implementation. For clarity of presentation, the QP used for
    % computing the closed-loop control law is included in this class
    % instead of solving a single QP. This class acts also as a wrapper for
    % the OCP solver.

    properties (SetAccess = private)
        ocpSol % cell array, solution from the main QP problem
        lambdaSol % solution of convex combination of the vertex control law

        hausDist % Hausdorff distance from state constraint ||v_x_i -z0||
        convh_x0s
    end

    properties (Access = private)
        OCP % YALMIP optimizer, Optimal Controller
        mpcLaw % YALMIP optimizer, compute the input in state (and not parameter) space
    end

    properties(Dependent)
        Lyapunov_cost % as defined in [(ii),Corollary 1] from paper

        rciSol % cell array, solution from RCI
        y0_sol % compute the first step closed-loop parameter vector
    end

    methods % GETTER methods
        function Lyap_cost = get.Lyapunov_cost(obj)
            Lyap_cost = obj.ocpSol{end};
        end

        function rciSol = get.rciSol(obj)
            rciSol = {obj.OCP.ys; obj.OCP.us};
        end

        function y0Sol = get.y0_sol(obj)
            y0Sol = obj.ocpSol{1}(:,1);
        end

        function hausDist = get.hausDist(obj)
            if isnan(obj.hausDist)
                disp("Computing Hausdorff distance...")
                obj.getHausdorffDist();
            end
            hausDist = obj.hausDist;
        end
    end


    methods (Access = public)
        function obj = CCTMPC(sys, ccPoly, CFM, N, gamma, sol_opts, compileFunc, y_m,u_m)
            if nargin == 7
                [y_m,u_m] = precompute_mRCI_set(ccPoly,CFM.cost_RCI);
            end
            obj.OCP = CC_OCP_casadi(sys, ccPoly, CFM, N, gamma, sol_opts, compileFunc, y_m,u_m);
            obj.mpcLaw = obj.initMPCLaw(sol_opts,compileFunc);
            obj.ocpSol = cell(3,1); % {y_t,u_0,costOCP_t} OCP fields
            obj.hausDist = nan;
        end


        function uMPC = solve(obj, x_curr)
            % solve OCP (in parameter space)
            [yOpt,u0Opt,ocpCost] = obj.OCP.solve_OCP(x_curr);
            obj.ocpSol = {full(yOpt),full(u0Opt),full(ocpCost)};

            % compute the optimal input (in state space)
            uMPC = full( obj.mpcLaw(x_curr, obj.ocpSol{1}(:,2)) );% x0, y1Opt

            if isnan(uMPC)
                warning("uMPC returned NaN. Consider enabling debug mode.")
                feas_x = full( obj.findFeasibleX0(x_curr,eye(obj.OCP.sys.nx)));
                if sum(abs(feas_x-x_curr) > 1e-5)
                    error("The current state does not lie inside the feasible"+...
                        " region of the control scheme. Consider changing initial condition.")
                end
            end
        end

        function feasRegion = computeFeasRegion(obj, n_facets)
            % interface method to the underlying OCP class
            feasRegion = obj.OCP.computeFeasRegion(n_facets);
        end

        function feasX0 = findFeasibleX0(obj,x_curr, weightMat)
            % interface method to the underlying OCP class
            if nargin < 3
                weightMat = eye(obj.OCP.sys.nx);
            end
            feasX0 = full(obj.OCP.getFeasX0(x_curr,weightMat));
        end

        function hausD = getHausdorffDist(obj,OuterSet)
            % parallel computation of hausdorff distance wrt OuterSet
            if nargin < 2
                OuterSet = obj.OCP.sys.X;
            end
            % normalize directions
            c_dirs = [OuterSet.V; OuterSet.A]; 
            c_dirs = (c_dirs./vecnorm(c_dirs,2,2))';  
            n_dirs = size(c_dirs,2);
            % initialize csdFun.map()
            parHausFun = obj.OCP.implicitSupportFun.map(n_dirs,'thread',feature('numcores'));
            
            [x0_s,supp_MPC] = parHausFun(c_dirs);
            x0_s = full(x0_s); supp_MPC = full(supp_MPC)'; 
            
            obj.convh_x0s = x0_s;

            % being a convex polytope, the max will be at one vertex or an edge
            supp_Outer = [max(c_dirs'*OuterSet.V',[],2)];

            hausD = max(supp_Outer-supp_MPC);
            if nargin < 2
            obj.hausDist = hausD;
            end
        end
    end

    methods (Access = private)
        function mpcLaw = initMPCLaw(obj,~,compileFun)
            % CasADi-based optimizer for computing the control-law.
            % Returns a casadi.Function with inputs {x0, y1Opt} and output u0Opt.

            % Get dimensions and necessary matrices
            f  = obj.OCP.ccPoly.f;
            nu = obj.OCP.sys.nu;
            nx = obj.OCP.sys.nx;
            nm = obj.OCP.sys.nm;

            % Retrieve matrices from the OCP structure
            F_mat = obj.OCP.ccPoly.F;  % Constraint matrix for the polytope (dimensions: ? x f)
            d_vec = obj.OCP.ccPoly.d;  % Vector offset for the polytope
            U_A   = obj.OCP.sys.U.A;    % Input constraint matrix
            U_b   = obj.OCP.sys.U.b;    % Input constraint bounds

            A_convh = obj.OCP.sys.A_convh; % cell array (each: nx x nx)
            B_convh = obj.OCP.sys.B_convh; % cell array (each: nx x nu)

            % Create a CasADi Opti instance
            opti = casadi.Opti('conic');

            % Decision variable: control input u0Opt (nu×1)
            u0Opt = opti.variable(nu,1);

            % Parameters (to be provided at each function call):
            x0    = opti.parameter(nx,1);   % initial state (parameter)
            y1Opt = opti.parameter(f,1);    % desired/target y1 (parameter)

            % For each model, enforce RCT condition:
            for i = 1:nm
                opti.subject_to(F_mat*(A_convh{i}*x0 + B_convh{i}*u0Opt) + d_vec <= y1Opt);
            end
            % Add the input constraint: U_A*u0Opt <= U_b.
            opti.subject_to(U_A*u0Opt <= U_b);

            % Define the cost function: squared 2-norm of u0Opt.
            cost_expr = u0Opt.' * u0Opt;
            opti.minimize(cost_expr);

            % set solver Options
            opti.solver('daqp')

            % Convert the fixed-structure problem into a CasADi Function.
            % This function can be repeatedly called with different x0 and y1Opt values.
            mpcLaw = opti.to_function('mpcLaw', {x0, y1Opt}, {u0Opt}).expand();

            if compileFun

                mfilePath = mfilename('fullpath');
                if contains(mfilePath,'LiveEditorEvaluationHelper')
                    mfilePath = matlab.desktop.editor.getActiveFilename;
                end
                bak_dir = cd(fileparts(mfilePath));

                mpcLaw.generate('mpcLaw.c',struct('mex', true));
                libPath = casadi.GlobalOptions.getCasadiPath();
                incPath = casadi.GlobalOptions.getCasadiIncludePath();
                mex('-v',...        % verbose mode
                    ['-I' incPath],...      % includes
                    ['-L' libPath],...      % used lib path
                    '-ldaqp',...           % osqp lib for linking
                    '-lcasadi',...
                    '-Dcasadi_real_min',...
                    'mpcLaw.c'...
                    );

                mpcLaw = str2func('mpcLaw');
                cd(bak_dir)
            end
        end

        


            

            % x_verts = obj.OCP.sys.X.V'; n_verts = size(x_verts,2);
            % 
            % % initialize csdFun.map()
            % parHausFun = obj.OCP.getFeasX0.map(n_verts,'thread',feature('numcores'));
            % 
            % x0_s_feas = parHausFun(x_verts,kron(ones(1,n_verts),eye(obj.OCP.sys.nx)));
            % obj.hausDist = max(vecnorm(x_verts-full(x0_s_feas)));
            % obj.convh_x0s = full(x0_s_feas);
        
    end

end