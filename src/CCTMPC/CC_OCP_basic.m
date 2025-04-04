classdef CC_OCP_basic < handle
    % OCP as defined in Equation (19) from the paper. The main difference
    % from the OptimalController class is the use of "optimize" function
    % instead of "optimizer" from YALMIP. This leads to a simpler but less
    % efficient code, as DynSystem is an handle class whose system matrices
    % are updated across all the instances.

    properties (SetAccess = private)
        sys % dynamical system
        N % prediction horizon length
        ccPoly % Configuration-Constrained Polytope
        sdp_opts % YALMIP optimizer options (solver, verbosity...)
        costFunMan % support class for cost function definition
        gamma % user specified discount factor

        feasRegion % Polytope, (approximated) feasible region for the OCP

        ys, us % parameter vector and vertex control law defining the mRCI set
        RCI_setcost % base cost to be subtracted to have a Lyapunov function
    end

    methods % GETTER methods
        function feasReg = get.feasRegion(obj)
            if class(obj.feasRegion) == "Polyhedron"
                feasReg = obj.feasRegion;
            else
                % f=100, often enough to get a good approximation
                obj.feasRegion = obj.computeFeasRegion();
                feasReg = obj.feasRegion;
            end
        end
    end

    methods (Access = public)
        function obj = CC_OCP_basic(sys, ccPoly, costFunMan, N, gamma, sdp_opts)
            obj.sys = sys;
            obj.ccPoly = ccPoly;
            obj.costFunMan = costFunMan;
            obj.N = N;
            obj.gamma = gamma;
            obj.sdp_opts = sdp_opts;

            % pre-compute the mRCI to which we converge
            obj.compute_mRCI_set();

            obj.feasRegion = nan; % computation (slow) delayed, as only used in plots
        end

        function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
            % find the closest feasible initial state to the desired one

            % OCP_optimizer YALMIP definitions
            x0=sdpvar(obj.sys.nx, 1,'full');
            y=sdpvar(obj.ccPoly.f, obj.N,'full');
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

            % solve the OCP
            optimize(constr,cost,obj.sdp_opts);
            feas_x0 = value(x0);
        end


        function ocpSol = solve_OCP(obj, x0)
            % compute the main QP solution
            y=sdpvar(obj.ccPoly.f, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

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
                constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1) )];
            end
            constr = [constr; obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys) )];

            % solve the OCP
            optimize(constr,cost,obj.sdp_opts);

            % return optimization problem variables and cost
            ocpSol = {value(y),value(u),value(cost)};
        end


        function feasRegion = computeFeasRegion(obj, n_facets)
            % compute a new feasible region given a template of #facets
            % complexity.
            if nargin == 1
                n_facets = 100;
            end
            if isscalar(n_facets)
                verts = polyGenEulerAngles(obj.sys.nx, n_facets);
                H_feas = Polyhedron(verts,ones(size(verts,1),1)).minHRep().A;
            else
                assert(size(n_facets,2)==obj.sys.nx,"Wrong dimensional template")
                H_feas = n_facets; % we pass a template directly
                n_facets = size(H_feas,1);

            end

            h_feas = zeros(size(H_feas,1),1);

            x_convh = zeros(n_facets,obj.sys.nx);

            x0 = sdpvar(obj.sys.nx,1,'full');
            y=sdpvar(obj.ccPoly.f, obj.N,'full');
            u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');


            for i=1:size(H_feas,1)
                % cost function
                cost = -H_feas(i,:)*x0; % maximize the direction

                % % Constraints (look Equation (19) from paper)
                constr = [];

                % Initial condition
                constr = [constr; obj.ccPoly.F*x0 <= y(:,1)];
                % Stage + terminal constraints
                for k=1:obj.N-1
                    constr = [constr; obj.constrSetS(y(:,k), u(:,:,k), y(:,k+1) )];
                end
                constr = [constr; obj.constrSetS(y(:,obj.N), u(:,:,obj.N), obj.ys + obj.gamma*(y(:,obj.N)-obj.ys) )];

                % solve the OCP
                optimize(constr,cost,obj.sdp_opts);

                % return the height on the specific direction
                h_feas(i) = -value(cost);
                x_convh(i,:) = value(x0)';
            end
            % TODO: which approx? Inner or outer?
            % feasRegion = Polyhedron(x_convh).minHRep().minVRep(); % inner
            feasRegion = Polyhedron(H_feas, h_feas).minHRep().minVRep(); % outer
        end

    end

    methods (Access = private)

        function constr = constrSetS(obj,y,u,yp)
            % define the Configuration Constrained RFIT set.
            % NOTE: if the system convex hulls are updated, the constraints
            % do it accordingly
            constr = [];
            for j=1:obj.ccPoly.v
                for i=1:obj.sys.nm % #models
                    constr = [constr;
                        obj.ccPoly.F*(obj.sys.A_convh{i}*obj.ccPoly.Vi_s{j}*y + obj.sys.B_convh{i}*u(:,j) )+ obj.ccPoly.d <= yp];
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

    end

end