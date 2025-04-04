classdef H_SP_OCP_basic < handle
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

        ys, us % parameter vector defining the fixed mRCI set
        RCI_setcost % base cost to be subtracted to have a Lyapunov function

        y_bar_s % cell array of nm support functions
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
        function obj = H_SP_OCP_basic(sys, ccPoly, costFunMan, N, gamma, sdp_opts)
            obj.sys = sys;
            obj.ccPoly = ccPoly;
            obj.costFunMan = costFunMan;
            obj.N = N;
            obj.gamma = gamma;
            obj.sdp_opts = sdp_opts;

           % pre-compute the mRCI to which we converge
            obj.compute_mRCI_set();

            % pre-compute the support function for the parameter vector
            obj.computeSupportDynRPI();
            
            obj.feasRegion = nan; % computation (slow) delayed, as only used in plots
        end

        function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
            % find the closest feasible initial state to the desired one

            % OCP_optimizer YALMIP definitions
            x0=sdpvar(obj.sys.nx, 1,'full');
            z = sdpvar(obj.sys.nx, obj.N,'full');
            alpha = sdpvar(1, obj.N,'full');
            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            v = sdpvar(obj.sys.nu, obj.N,'full');
            % u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

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
            constr = [constr; obj.ccPoly.F*(x0 - z(:,1)) <= alpha(1)*obj.ys];

            % Stage + terminal constraints
            for k=1:obj.N-1
                constr = [constr;
                    obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)} )];
            end
            constr = [constr;
                obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1} )];

            % solve the OCP
            optimize(constr,cost,obj.sdp_opts);
            feas_x0 = value(x0);
        end


        function ocpSol = solve_OCP(obj, x0)
            % compute the main QP solution
            z = sdpvar(obj.sys.nx, obj.N,'full');
            alpha = sdpvar(1, obj.N,'full');
            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            v = sdpvar(obj.sys.nu, obj.N,'full');
            % u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

            % % Cost Function
            cost = 0;
            for t=1:obj.N-1
                y_t = obj.ccPoly.F*z(:,t) + alpha(t)*obj.ys;
                u_t = v(:,t) + alpha(t)*obj.us;
                cost = cost + obj.costFunMan.cost_L_k(y_t, u_t, obj.ys, obj.us);
            end
            y_N = obj.ccPoly.F*z(:,obj.N) + alpha(obj.N)*obj.ys;
            u_N = v(:,obj.N) + alpha(obj.N)*obj.us;
            cost = cost + obj.costFunMan.cost_L_N(y_N, u_N, obj.ys, obj.us);

            % % Constraints (look Equation (19) from paper)
            constr = [];

            % Initial condition
            constr = [constr; obj.ccPoly.F*(x0 - z(:,1)) <= alpha(1)*obj.ys];

            % Stage + terminal constraints
            for k=1:obj.N-1
                constr = [constr;
                    obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)} )];
            end
            constr = [constr;
                obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1} )];

            % solve the OCP
            optimize(constr,cost,obj.sdp_opts);

            % return optimization problem variables and cost
            ocpSol = {value(z),value(alpha),value(u),value(cost)};
        end


        function feasRegion = computeFeasRegion(obj, n_facets)
            % compute a new feasible region given a template of #facets
            % complexity.
            if nargin == 1
                n_facets = 100;
            end
            % H_feas = getTemplateMatrix(obj.sys.nx, n_facets);
            verts = polyGenEulerAngles(obj.sys.nx, n_facets);
            H_feas = Polyhedron(verts,ones(size(verts,1),1)).minHRep().A;
            h_feas = zeros(size(H_feas,1),1);

            x_convh = zeros(n_facets,obj.sys.nx);

            x0 = sdpvar(obj.sys.nx,1,'full');
            z = sdpvar(obj.sys.nx, obj.N,'full');
            alpha = sdpvar(1, obj.N,'full');
            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            v = sdpvar(obj.sys.nu, obj.N,'full');
            % u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');


            for i=1:size(H_feas,1)
                % cost function
                cost = -H_feas(i,:)*x0; % maximize the direction

                % % Constraints (look Equation (19) from paper)
                constr = [];

                % Initial condition
                constr = [constr; obj.ccPoly.F*(x0 - z(:,1)) <= alpha(1)*obj.ys];

                % Stage + terminal constraints
                for k=1:obj.N-1
                    constr = [constr;
                        obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)} )];
                end
                constr = [constr;
                    obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1} )];

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

            % compute the mRCI for the given template
            y_m = sdpvar(obj.ccPoly.f,1);
            u_m = sdpvar(obj.sys.nu, obj.ccPoly.v,'full');
            
            % cost function
            cost = obj.costFunMan.cost_RCI(y_m,u_m);

            % constraints
            constr = obj.constrSetS(y_m,u_m,y_m);

            % solve the OCP
            optimize(constr,cost, obj.sdp_opts);

            obj.ys = value(y_m); obj.us = value(u_m);
            obj.RCI_setcost = value(cost);
        end


        function constr = constrSetS_H_SP(obj, za, v, za_p)
            % define the Homothetic 1-step propagation set.
            z = za{1}; alpha = za{2};
            zp = za_p{1}; alpha_p = za_p{2};

            constr = [];
            for i=1:obj.sys.nm % #models
                constr = [constr;
                    obj.ccPoly.F*(obj.sys.A_convh{i}*z + obj.sys.B_convh{i}*v) + ...
                    alpha*obj.y_bar_s{i} + obj.ccPoly.d <= (obj.ccPoly.F*zp + alpha_p*obj.ys)];
            end

            for j=1:obj.ccPoly.v
                constr = [constr;   obj.sys.X.A*(z + alpha*obj.ccPoly.Vi_s{j}*obj.ys) <= obj.sys.X.b;
                    obj.sys.U.A*(v + alpha*obj.us(:,j)) <= obj.sys.U.b];
            end
            constr = [constr; alpha >= 0];
        end


        function computeSupportDynRPI(obj)
            % use MPT3 methods
            obj.y_bar_s = cell(1,obj.sys.nm);

            for i=1:obj.sys.nm
                Xs_p_verts = zeros(obj.sys.nx,obj.ccPoly.v);
                for j=1:obj.ccPoly.v
                    Xs_p_verts(:,j) = obj.sys.A_convh{i}*obj.ccPoly.Vi_s{j}*obj.ys + obj.sys.B_convh{i}*obj.us(:,j);
                end
                Xs_plus = Polyhedron(Xs_p_verts');
                obj.y_bar_s{i} = Xs_plus.support(obj.ccPoly.F');
            end
        end

    end

end