classdef Partial_CC_OCP < handle
    % OCP as defined in Equation (19) from the paper. Problem written using
    % YALMIP. Preferred solver is Gurobi. We also define a flag variable
    % to disable matrix update if not required and speeding up computation
    % (2-3x improvement).

    properties (SetAccess = private)
        sys % dynamical system
        N % prediction horizon length
        ccPoly % Configuration-Constrained Polytope
        sdp_opts % YALMIP optimizer options (solver, verbosity...)
        costFunMan % support class for cost function definition
        solve_OCP  % main QP problem
        solve_DRTO % QP for solving RCI set only. Used for plotting reasons
        gamma % user specified discount factor

        feasRegion % Polytope, (approximated) feasible region for the OCP

        ys, us % parameter vector and vertex control law defining the mRCI set
        RCI_setcost % base cost to be subtracted to have a Lyapunov function

        H_verts, H_params % cell arrays of indices for the vertices (and then their
        % relative parameters) for which we enforce an Homothetic transformation
        E_partial
    end

    properties (Access = private)
        directionalFeasRegion % YALMIP optimizer, compute the farthest feasible
        % initial state with respect to a specific facet.
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


    methods (Access = public)
        function obj = Partial_CC_OCP(sys, ccPoly, H_verts, CFM, N, gamma, sdp_opts, var_convh)
            obj.sys = sys;
            obj.N = N;
            obj.ccPoly = ccPoly;
            obj.sdp_opts = sdp_opts;
            obj.costFunMan = CFM;
            obj.gamma = gamma;

            obj.H_verts = H_verts;

            % pre-compute the mRCI to which we converge
            obj.compute_mRCI_set();

            % find vertex facet relationship
            obj.findIndices_yOpt();

            obj.recomputeConicE();

            % optimizers initialization
            obj.solve_OCP = obj.init_OCP(var_convh);

            obj.directionalFeasRegion = obj.initDirFeasRegion(var_convh);
            obj.feasRegion = nan; % computation (slow) delayed, as only used in plots
        end


        function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
            % find the closest feasible initial state to the desired one

            f = obj.ccPoly.f;
            y = double2sdpvar(zeros(f,obj.N));

            z = sdpvar(obj.sys.nx, obj.N,'full'); 
            y_hat = sdpvar(f-length(obj.H_params),obj.N,'full');
            alpha = sdpvar(1,obj.N,'full');

            CC_params = setdiff(1:f,obj.H_params);
            y(CC_params,:) = y_hat; % Configuration-Constrained
            % y(obj.H_params,:) = obj.ccPoly.F(obj.H_params,:)*z + obj.ys(obj.H_params)*alpha; % Homothetic
            y_tmp = arrayfun(@(t) alpha(t)*obj.ys(obj.H_params), 1:obj.N, "UniformOutput",false);
            y(obj.H_params,:) = cat(2,y_tmp{:});


            % add the translation part
            y = y + obj.ccPoly.F*z;

            % OCP_optimizer YALMIP definitions
            x0=sdpvar(obj.sys.nx, 1,'full');
            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');


            switch class(costFun)
                case 'double'
                    cost = weighted2NormSquared(x0 - x0_des, costFun);
                case 'function_handle'
                    cost = costFun(x0, x0_des);
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(Stage_cost))
            end

            % % Constraints (look Equation (19) from paper)
            constr = [];

            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];

            % Stage + terminal constraints
            for k=1:obj.N-1
                constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1) )];
            end
            constr = [constr; obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys) )];

            % C-Set constraint
            % constr = [constr; alpha(:)>=0];
            % solve the OCP
            diagn = optimize(constr,cost,obj.sdp_opts);
            % diagn.solverinput.model
            feas_x0 = value(x0);
        end


        function feasRegion = computeFeasRegion(obj, n_facets, varargin)
            % compute a new feasible region given a template of #facets
            % complexity. varargin used to pass new (A_convh, B_convh)
            % H_feas = getTemplateMatrix(obj.sys.nx, n_facets);
            verts = polyGenEulerAngles(obj.sys.nx, n_facets);
            H_feas = Polyhedron(verts,ones(size(verts,1),1)).minHRep().A;
            h_feas = zeros(size(H_feas,1),1);

            x_convh = zeros(n_facets,obj.sys.nx);


            for i=1:size(H_feas,1)
                switch length(varargin)
                    case 0
                        x0_max = obj.directionalFeasRegion(H_feas(i,:)');
                    case 2
                        x0_max = obj.directionalFeasRegion(H_feas(i,:)',varargin{:});
                    otherwise % nargin-1 as obj is counted as argument
                        error("Wrong number of arguments: %d. Only 2 and 4 supported.", nargin-1);
                end
                h_feas(i) = H_feas(i,:)*x0_max;
                x_convh(i,:) = x0_max';

            end
            % TODO: which approx? Inner or outer?
            % feasRegion = Polyhedron(x_convh).minHRep().minVRep(); % inner
            feasRegion = Polyhedron(H_feas, h_feas).minHRep().minVRep(); % outer
        end
    end


    methods (Access = private)
        function ocpOptim = init_OCP(obj, var_convh)
            % define the YALMIP optimizer for the main QP
            f = obj.ccPoly.f;
            y = double2sdpvar(zeros(f,obj.N));

            z = sdpvar(obj.sys.nx, obj.N,'full'); 
            y_hat = sdpvar(f-length(obj.H_params),obj.N,'full');
            alpha = sdpvar(1,obj.N,'full');

            CC_params = setdiff(1:f,obj.H_params);
            y(CC_params,:) = y_hat; % Configuration-Constrained
            % y(obj.H_params,:) = obj.ccPoly.F(obj.H_params,:)*z + obj.ys(obj.H_params)*alpha; % Homothetic
            y_tmp = arrayfun(@(t) alpha(t)*obj.ys(obj.H_params), 1:obj.N, "UniformOutput",false);
            y(obj.H_params,:) = cat(2,y_tmp{:});


            % add the translation part
            y = y + obj.ccPoly.F*z;

            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

            x0=sdpvar(obj.sys.nx, 1,'full');

            if var_convh
                A_sdp = cell(1,obj.sys.nm);
                for i = 1:obj.sys.nm
                    A_sdp{i} = sdpvar(obj.sys.nx,obj.sys.nx,'full');
                end
                B_sdp = cell(1,obj.sys.nm);
                for i = 1:obj.sys.nm
                    B_sdp{i} = sdpvar(obj.sys.nx,obj.sys.nu,'full');
                end
            end

            % % Cost Function
            cost = 0;
            for t=1:obj.N-1
                cost = cost + obj.costFunMan.cost_L_k(y(:,t), u(:,:,t), obj.ys, obj.us);
            end
            cost = cost + obj.costFunMan.cost_L_N(y(:,obj.N), u(:,:,obj.N), obj.ys, obj.us);

            % % Constraints (look Equation (19) from paper)
            constr = [];

            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];

            % Stage + terminal constraints
            for k=1:obj.N-1
                if var_convh
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1), A_sdp,B_sdp )];
                else
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1))];
                end
            end
            if var_convh
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys), A_sdp,B_sdp)];
            else
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys))];
            end

            % C-Set constraint
            constr = [constr; alpha(:)>=0];
            % instantiate the Optimizer
            if var_convh
                ocpOptim = optimizer(constr,cost,...
                    obj.sdp_opts,...
                    [{x0},A_sdp(:)',B_sdp(:)'],...% input parameters
                    {y,u,cost});          % output parameters
            else
                ocpOptim = optimizer(constr,cost,...
                    obj.sdp_opts,...
                    {x0},...          % input parameters
                    {y,u,cost});  % output parameters
            end
        end


        function dirFeasRegion = initDirFeasRegion_bak(obj, var_convh)
            % define the YALMIP optimizer for maximization of the feasible
            % region with respect to a specific facet
            c_vec = sdpvar(obj.sys.nx,1,'full');
            x0 = sdpvar(obj.sys.nx,1,'full');

            % cost function
            cost = -c_vec'*x0; % maximize the direction

            f = obj.ccPoly.f; n_H = length(obj.H_params);
            alpha = double2sdpvar(zeros(f,obj.N));

            CC_params = setdiff(1:f,obj.H_params);
            alpha(CC_params,:) = sdpvar(f-n_H,obj.N,'full');
            alpha(obj.H_params,:) = repmat(sdpvar(1,obj.N,'full'),[n_H,1]);

            y = arrayfun(@(t) diag(alpha(:,t))*obj.ys,1:obj.N,"UniformOutput",false);
            y = [y{:}];

            % add the translation part
            z = sdpvar(obj.sys.nx, obj.N,'full');
            y = y + obj.ccPoly.F*z;

            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

            if var_convh
                A_sdp = cell(1,obj.sys.nm);
                for i = 1:obj.sys.nm
                    A_sdp{i} = sdpvar(obj.sys.nx, obj.sys.nx,'full');
                end
                B_sdp = cell(1,obj.sys.nm);
                for i = 1:obj.sys.nm
                    B_sdp{i} = sdpvar(obj.sys.nx, obj.sys.nu,'full');
                end
            end

            % % Constraints (look Equation (19) from paper)
            constr = [];

            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];

            % Stage + terminal constraints
            for k=1:obj.N-1
                if var_convh
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1), A_sdp,B_sdp )];
                else
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1))];
                end
            end
            if var_convh
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys), A_sdp,B_sdp)];
            else
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys) )];
            end

            % C-Set constraint
            constr = [constr; alpha(:)>=0];
            % instantiate the Optimizer
            if var_convh
                dirFeasRegion = optimizer(constr, cost,...
                    obj.sdp_opts,...
                    [{c_vec},A_sdp(:)',B_sdp(:)'],... % input parameters
                    {x0});      % output parameters
            else
                dirFeasRegion = optimizer(constr, cost,...
                    obj.sdp_opts,...
                    {c_vec},... % input parameters
                    {x0});      % output parameters
            end
        end

        function dirFeasRegion = initDirFeasRegion(obj, var_convh)
            % define the YALMIP optimizer for maximization of the feasible
            % region with respect to a specific facet
            c_vec = sdpvar(obj.sys.nx,1,'full');
            x0 = sdpvar(obj.sys.nx,1,'full');

            % cost function
            cost = -c_vec'*x0; % maximize the direction

            f = obj.ccPoly.f;
            y = double2sdpvar(zeros(f,obj.N));

            z = sdpvar(obj.sys.nx, obj.N,'full'); 
            y_hat = sdpvar(f-length(obj.H_params),obj.N,'full');
            alpha = sdpvar(1,obj.N,'full');

            CC_params = setdiff(1:f,obj.H_params);
            y(CC_params,:) = y_hat; % Configuration-Constrained
            % y(obj.H_params,:) = obj.ccPoly.F(obj.H_params,:)*z + obj.ys(obj.H_params)*alpha; % Homothetic
            y_tmp = arrayfun(@(t) alpha(t)*obj.ys(obj.H_params), 1:obj.N, "UniformOutput",false);
            y(obj.H_params,:) = cat(2,y_tmp{:});

            % % add the translation part
            y = y + obj.ccPoly.F*z;

            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

            if var_convh
                A_sdp = cell(1,obj.sys.nm);
                for i = 1:obj.sys.nm
                    A_sdp{i} = sdpvar(obj.sys.nx, obj.sys.nx,'full');
                end
                B_sdp = cell(1,obj.sys.nm);
                for i = 1:obj.sys.nm
                    B_sdp{i} = sdpvar(obj.sys.nx, obj.sys.nu,'full');
                end
            end

            % % Constraints (look Equation (19) from paper)
            constr = [];

            % Initial condition
            constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];

            % Stage + terminal constraints
            for k=1:obj.N-1
                if var_convh
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1), A_sdp,B_sdp )];
                else
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1))];
                end
            end
            if var_convh
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys), A_sdp,B_sdp)];
            else
                constr = [constr;
                    obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys) )];
            end

            % C-Set constraint
            constr = [constr; alpha(:)>=0];
            % instantiate the Optimizer
            if var_convh
                dirFeasRegion = optimizer(constr, cost,...
                    obj.sdp_opts,...
                    [{c_vec},A_sdp(:)',B_sdp(:)'],... % input parameters
                    {x0});      % output parameters
            else
                dirFeasRegion = optimizer(constr, cost,...
                    obj.sdp_opts,...
                    {c_vec},... % input parameters
                    {x0});      % output parameters
            end
        end


        function constr = constrSetS(obj,y,u,yp,A_convh,B_convh)
            % define the Configuration-Constrained 1-step propagating set S.

            if nargin == 4 % models won't be updated
                A_convh = obj.sys.A_convh;
                B_convh = obj.sys.B_convh;
            end

            constr = [];
            for j=1:obj.ccPoly.v
                for i=1:obj.sys.nm % #models
                    constr = [constr;
                        obj.ccPoly.F*(A_convh{i}*obj.ccPoly.Vi_s{j}*y + B_convh{i}*u(:,j) )+ obj.ccPoly.d <= yp];
                end
                constr = [constr;   obj.sys.X.A*obj.ccPoly.Vi_s{j}*y <= obj.sys.X.b;    obj.sys.U.A*u(:,j) <= obj.sys.U.b];

            end
            constr = [constr; obj.ccPoly.E*y <= 0];
        end


        function compute_mRCI_set(obj)
            y_m = sdpvar(obj.ccPoly.f,1);
            u_m = sdpvar(obj.sys.nu, obj.ccPoly.v,'full');

            % cost function
            cost = obj.costFunMan.cost_RCI(y_m,u_m);

            % constraints
            constr = obj.constrSetS(y_m,u_m,y_m);

            % solve the OCP
            optimize(constr,cost,obj.sdp_opts);

            obj.ys = value(y_m);
            obj.us = value(u_m);
            obj.RCI_setcost = value(cost);
        end

        function findIndices_yOpt(obj)
            y_v_indices = zeros(obj.ccPoly.f,1);

            % find for each column the non zero indices
            for i=1:length(obj.H_verts)
                y_v_indices = y_v_indices + any(obj.ccPoly.Vi_s{obj.H_verts(i)}, 1)';
            end
            % we now also have the multiplicity of each parameter. For our
            % purposes, if we find a nonzero, those are the parameters for
            % Homothetic transformation
            obj.H_params = find(y_v_indices);
        end

        function recomputeConicE(obj)


        end

    end

end














