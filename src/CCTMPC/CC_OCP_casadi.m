classdef CC_OCP_casadi < handle
    % CC_OCP_casadi implements the OCP (Equation (19) from the paper)
    % using CasADi’s Opti interface. It builds a fixed‐structure problem and
    % then converts it to a function (via to_function) for repeated calls.

    properties (SetAccess = private)
        sys         % dynamical system structure
        N           % prediction horizon length
        ccPoly      % Configuration-Constrained Polytope structure
        cfm         % instance of CostFunctionManagerCasadi
        gamma       % user-specified discount factor

        solve_OCP   % function handle returning {y, u, cost}

        ys          % mRCI parameter solution (computed offline)
        us          % corresponding control law (for mRCI)

        getFeasX0   % define Optimizer for feasible initial condition
        feasRegion % Polytope, (approximated) feasible region for the OCP

        implicitSupportFun % compute farthest x0 point for a given direction

        csd_opts
    end

    methods % GETTER methods
        function feasReg = get.feasRegion(obj)
            if class(obj.feasRegion) == "Polyhedron"
                feasReg = obj.feasRegion;
            else
                % Hard-coded #facets, often enough to get a good approximation
                obj.feasRegion = obj.computeFeasRegion(100);
                feasReg = obj.feasRegion;
            end
        end
    end



    methods
        function obj = CC_OCP_casadi(sys, ccPoly, costFunMan, N, gamma, csd_opts, compileFunc, y_m,u_m)
            obj.sys = sys;
            obj.ccPoly = ccPoly;
            obj.cfm = costFunMan;
            obj.N = N;
            obj.gamma = gamma;

            obj.csd_opts = csd_opts;

            if nargin == 7 % Compute the mRCI set via the CasADi-based routine.
                obj.compute_mRCI_set();
            else % use the given RCI set
                obj.ys = y_m; obj.us = u_m;
            end

            % Initialize Opt.Problems and get function handles
            obj.solve_OCP = obj.init_OCP(compileFunc);
            obj.getFeasX0 = obj.initFeasibleX0(compileFunc);

            obj.implicitSupportFun = obj.initImplicitSupportFun(compileFunc);
            obj.feasRegion = nan; % computation (slow) delayed, as only used in plots
        end


        function feasRegion = computeFeasRegion(obj, n_facets)
            % compute a feasible region given by n_facets directions.
            if isscalar(n_facets)
                % uniformly sampling (nx-1)-unit sphere for directions
                H_feas = numericalSpherePoints(obj.sys.nx, n_facets);
            else % passing already a template matrix
                assert(size(n_facets,2)==obj.sys.nx,"Wrong dimensional template")
                H_feas = n_facets;
            end
            % h_feas = zeros(size(H_feas,1),1);

            parSupportFun = obj.implicitSupportFun.map(size(H_feas,1),...
                'thread',feature('numcores'));
            [~,h_feas] = parSupportFun(H_feas');
            
            % for i=1:size(H_feas,1)
            %     x0_max = full(obj.implicitSupportFun(H_feas(i,:)'));
            %     h_feas(i) = H_feas(i,:)*x0_max;
            % end
            feasRegion = Polyhedron(... % RoA  = X \cap feasReg
                obj.sys.X.intersect(Polyhedron(H_feas, full(h_feas)')).minHRep.V);
        end

    end

    methods (Access = private)
        function ocpFunc = init_OCP(obj, compileFun)
            % Build the OCP using CasADi's Opti.
            %
            % Decision variables:
            %   y ∈ ℝ^(f x N)    -- parameter trajectory
            %   u ∈ ℝ^(nu x (v) x N)  -- vertex control law at each time step
            %
            % Parameter:
            %   x0 ∈ ℝ^(nx x 1)   -- current state
            %
            % The constraints include:
            %   - initial condition: ccPoly.F * x0 <= y(:,1)
            %   - stage constraints: for each k = 1,...,N-1, the constraint set S is
            %       enforced using our helper constrSetS_casadi.
            %   - terminal constraint: y(:,N) must satisfy a terminal condition,
            %       here expressed using the computed mRCI parameter (ys) and a discount.
            %
            % The cost is constructed as
            %   ∑_{t=1}^{N-1} cost_L_k(y(:,t), u(:,:,t)) + cost_L_N(y(:,N), u(:,:,N))

            % Create an Opti instance.
            opti = casadi.Opti('conic');

            % Decision variables.
            y = opti.variable(obj.ccPoly.f, obj.N);
            u = arrayfun(@(k) opti.variable(obj.sys.nu, obj.ccPoly.v),1:obj.N,"UniformOutput",false);

            % Cost function
            cost = 0;
            for t = 1:obj.N-1
                cost = cost + obj.cfm.cost_L_k_CC(y(:,t), u{t});
            end
            cost = cost + obj.cfm.cost_L_N_CC(y(:,obj.N), u{obj.N});
            opti.minimize(cost);

            % Parameter: initial state.
            x0 = opti.parameter(obj.sys.nx, 1);

            % Initial constraint
            opti.subject_to(obj.ccPoly.F * x0 <= y(:,1));
            % Stage + Terminal constraints
            for k = 1:obj.N-1
                opti.subject_to(obj.constrSetS(y(:,k), u{k}, y(:,k+1)) )
            end
            opti.subject_to(obj.constrSetS(y(:,obj.N), u{obj.N}, obj.ys+obj.gamma*(y(:,obj.N)-obj.ys)))

            % Set solver options (here we consider to have daqp)
            opti.solver(obj.csd_opts{:})

            % Convert the fixed structure into a CasADi function.
            % The function takes one input: x0_param, and returns {y,u,cost}.
            ocpFunc = opti.to_function('solveOCP', {x0}, {y, u{1}, cost},...
                {'x0'},{'y','u_t=0','cost'}).expand();

            if compileFun

                mfilePath = mfilename('fullpath');
                if contains(mfilePath,'LiveEditorEvaluationHelper')
                    mfilePath = matlab.desktop.editor.getActiveFilename;
                end
                bak_dir = cd(fileparts(mfilePath));

                ocpFunc.generate('ocpFunc.c',struct('mex', true));
                libPath = casadi.GlobalOptions.getCasadiPath();
                incPath = casadi.GlobalOptions.getCasadiIncludePath();
                mex('-v',...        % verbose mode
                    ['-I' incPath],...      % includes
                    ['-L' libPath],...      % used lib path
                    '-ldaqp',...           % solver lib for linking
                    '-lcasadi',...
                    '-Dcasadi_real_min',...
                    'ocpFunc.c'...
                    );

                ocpFunc = str2func('ocpFunc');
                cd(bak_dir)
            end
        end

        function feasX0 = initFeasibleX0(obj, compileFun)
            % Create a CasADi Opti instance for finding a feasible initial state.
            opti = casadi.Opti('conic');

            % Decision variable: feasible initial state x0 (dimension: nx×1)
            x0_var = opti.variable(obj.sys.nx, 1);
            % Decision variables used to impose the OCP constraints:
            % y: parameter trajectory (dimension: f×N)
            y = opti.variable(obj.ccPoly.f, obj.N);
            % u: vertex control laws; we use a cell array, one per time step
            u = arrayfun(@(k) opti.variable(obj.sys.nu, obj.ccPoly.v), 1:obj.N, "UniformOutput", false);

            % Parameters for the optimizer: desired initial state and cost weighting.
            x0_des = opti.parameter(obj.sys.nx, 1);
            costMat = opti.parameter(obj.sys.nx, obj.sys.nx);

            % Cost: weighted squared norm of (x0_var - x0_des)
            cost = squared2norm(x0_var - x0_des, costMat);
            opti.minimize(cost);

            % Initial constraint
            opti.subject_to(obj.ccPoly.F * x0_var <= y(:,1));
            % Stage + Terminal constraints
            for k = 1:obj.N-1
                opti.subject_to(obj.constrSetS(y(:,k), u{k}, y(:,k+1)) )
            end
            opti.subject_to(obj.constrSetS(y(:,obj.N), u{obj.N}, obj.ys+obj.gamma*(y(:,obj.N)-obj.ys)))

            % Set solver options. NOTE: we need to handle semidefinite cost
            % opti.solver('daqp',struct(),struct('eps_prox',1e-3))
            opti.solver('gurobi',struct(),struct('outputflag',0))



            % Convert the fixed structure into a CasADi function.
            % This function accepts two inputs: the desired x0 (x0_des) and costMatrix,
            % and returns {x0_var, cost_expr}.
            feasX0 = opti.to_function('feasX0', ...
                {x0_des, costMat}, {x0_var},...
                {'x0_des','weightMat'},{'x0_feas'}).expand();

            if compileFun

                mfilePath = mfilename('fullpath');
                if contains(mfilePath,'LiveEditorEvaluationHelper')
                    mfilePath = matlab.desktop.editor.getActiveFilename;
                end
                bak_dir = cd(fileparts(mfilePath));

                feasX0.generate('feasX0.c',struct('mex', true));
                libPath = casadi.GlobalOptions.getCasadiPath();
                incPath = casadi.GlobalOptions.getCasadiIncludePath();
                mex('-v',...        % verbose mode
                    ['-I' incPath],...      % includes
                    ['-L' libPath],...      % used lib path
                    '-ldaqp',...           % solver lib for linking
                    '-lcasadi',...
                    '-Dcasadi_real_min',...
                    'feasX0.c'...
                    );

                feasX0 = str2func('feasX0');
                cd(bak_dir)

            end
        end

        function implSupportFun = initImplicitSupportFun(obj,compileFun)
            % Compute the farthest feasible x0 in a given direction.
            opti = casadi.Opti('conic');

            % Input parameter: the direction vector (c_vec, nx×1)
            c_vec = opti.parameter(obj.sys.nx,1);

            % Decision variable: candidate feasible initial state x0 (nx×1)
            x0 = opti.variable(obj.sys.nx,1);

            % Additional decision variables to enforce the OCP constraints:
            % y: parameter trajectory (f×N)
            y = opti.variable(obj.ccPoly.f, obj.N);
            % u: vertex control laws (nu×v) for each time step (store in cell array)
            u = arrayfun(@(k) opti.variable(obj.sys.nu, obj.ccPoly.v), 1:obj.N, "UniformOutput", false);

            % Objective: maximize c_vec' * x0
            cost = c_vec.' * x0;
            opti.minimize( - cost );

            % Constraints:
            % Initial condition: polytope constraint on x0, enforced via the first stage parameter y(:,1)
            opti.subject_to( obj.ccPoly.F * x0 <= y(:,1) );
            % Stage + Terminal constraints
            for k = 1:obj.N-1
                opti.subject_to(obj.constrSetS(y(:,k), u{k}, y(:,k+1)) )
            end
            opti.subject_to(obj.constrSetS(y(:,obj.N), u{obj.N}, obj.ys+obj.gamma*(y(:,obj.N)-obj.ys)))

            % Set solver options. NOTE: we need to handle semidefinite cost
            % opti.solver('daqp',struct(),struct('eps_prox',1e-3))
            opti.solver('gurobi',struct(),struct('outputflag',0))



            % Convert the problem into a CasADi function.
            % Input: direction vector c_vec; Output: corresponding x0.
            implSupportFun = opti.to_function('implSupportFun', {c_vec}, {x0,cost},...
                {'c_dir'},{'x0','supp_cost'}).expand();
            
            if compileFun
                mfilePath = mfilename('fullpath');
                if contains(mfilePath,'LiveEditorEvaluationHelper')
                    mfilePath = matlab.desktop.editor.getActiveFilename;
                end
                bak_dir = cd(fileparts(mfilePath));

                implSupportFun.generate('implSupportFun.c',struct('mex', true));
                libPath = casadi.GlobalOptions.getCasadiPath();
                incPath = casadi.GlobalOptions.getCasadiIncludePath();
                mex('-v',...        % verbose mode
                    ['-I' incPath],...      % includes
                    ['-L' libPath],...      % used lib path
                    '-ldaqp',...           % solver lib for linking
                    '-lcasadi',...
                    '-Dcasadi_real_min',...
                    'implSupportFun.c'...
                    );

                implSupportFun = str2func('implSupportFun');
                cd(bak_dir)
            end
        end

        function compute_mRCI_set(obj)
            opti = casadi.Opti('conic');

            % Define decision variables with CasADi
            y_m = opti.variable(obj.ccPoly.f, 1);
            u_m = opti.variable(obj.sys.nu, obj.ccPoly.v);

            % Define cost function
            cost = obj.cfm.cost_RCI(y_m, u_m);

            % Build constraints
            opti.subject_to(obj.constrSetS(y_m, u_m, y_m));

            % Define the optimization problem
            opti.minimize(cost);

            % opti.solver(obj.csd_opts{:})
            % it's often convenient use gurobi when #vertices is large
            opti.solver('gurobi',struct(),struct('outputflag',0))

            sol = opti.solve();
            obj.ys = sol.value(y_m); obj.us = sol.value(u_m);
        end

        function constr = constrSetS(obj, y, u, yp)
            % This function replicates the Configuration-Constrained RFIT set constraints in CasADi.
            % It returns a cell array of constraint expressions.
            constr = {};
            for j = 1:obj.ccPoly.v   % loop over vertices
                for i = 1:obj.ccPoly.sys.nm   % loop over models
                    constr{end+1} = obj.ccPoly.F * (obj.ccPoly.sys.A_convh{i} * obj.ccPoly.Vi_s{j} * y + ...
                        obj.ccPoly.sys.B_convh{i} * u(:,j)) + obj.ccPoly.d <= yp;
                end
                constr{end+1} = obj.ccPoly.sys.X.A * obj.ccPoly.Vi_s{j} * y <= obj.ccPoly.sys.X.b;
                constr{end+1} = obj.ccPoly.sys.U.A * u(:,j) <= obj.ccPoly.sys.U.b;
            end
            constr{end+1} = obj.ccPoly.E * y <= 0;
        end


    end
end
