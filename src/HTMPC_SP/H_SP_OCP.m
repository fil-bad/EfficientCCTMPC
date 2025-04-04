classdef H_SP_OCP < handle
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
        gamma % user specified discount factor

        feasRegion % Polytope, (approximated) feasible region for the OCP

        ys, us % parameter vector defining the fixed mRCI set
        RCI_setcost % base cost to be subtracted to have a Lyapunov function

        y_bar_s % cell array of nm support functions
        hx_m, hu_m % precompute state and input constraints

        findFeasibleX0 % define Optimizer for feasible initial condition
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
        function obj = H_SP_OCP(sys, ccPoly, costFunMan, N, gamma, sdp_opts, var_convh)
            obj.sys = sys;
            obj.N = N;
            obj.ccPoly = ccPoly;
            obj.sdp_opts = sdp_opts;
            obj.costFunMan = costFunMan;
            obj.gamma = gamma;

            % pre-compute the mRCI to which we converge
            obj.compute_mRCI_set();

            % pre-compute the support functions for the parameter vector
            obj.precomputeSupportFun();

            % optimizers initialization
            obj.solve_OCP = obj.init_OCP(var_convh);

            obj.findFeasibleX0 = obj.initFeasibleX0();


            obj.directionalFeasRegion = obj.initDirFeasRegion(var_convh);
            obj.feasRegion = nan; % computation (slow) delayed, as only used in plots
        end


        % function feas_x0 = findFeasibleX0(obj, x0_des, costFun)
        %     % find the closest feasible initial state to the desired one
        % 
        %     % OCP_optimizer YALMIP definitions
        %     x0=sdpvar(obj.sys.nx, 1,'full');
        %     z = sdpvar(obj.sys.nx, obj.N,'full');
        %     alpha = sdpvar(1, obj.N,'full');
        %     % y=sdpvar(obj.ccPoly.f, obj.N,'full');
        %     v = sdpvar(obj.sys.nu, obj.N,'full');
        %     % u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');
        % 
        %     switch class(costFun)
        %         case 'double'
        %             cost = weighted2NormSquared(x0 - x0_des, costFun);
        %         case 'function_handle'
        %             cost = costFun(x0, x0_des);
        %         otherwise
        %             error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(Stage_cost))
        %     end
        % 
        %     % % Constraints (look Equation (19) from paper)
        %     constr = [];
        % 
        %     % Initial condition
        %     constr = [constr; obj.ccPoly.F*(x0 - z(:,1)) <= alpha(1)*obj.ys];
        % 
        %     % Stage + terminal constraints
        %     for k=1:obj.N-1
        %         constr = [constr;
        %             obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)} )];
        %     end
        %     constr = [constr;
        %         obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1} )];
        % 
        %     % solve the OCP
        %     optimize(constr,cost,obj.sdp_opts);
        %     feas_x0 = value(x0);
        % end


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
            feasRegion = Polyhedron(H_feas, h_feas); % outer
        end
    end


    methods (Access = private)
        function ocpOptim = init_OCP(obj, var_convh)
            % define the YALMIP optimizer for the main QP
             z = sdpvar(obj.sys.nx, obj.N,'full');
            alpha = sdpvar(1, obj.N,'full');
            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            v = sdpvar(obj.sys.nu, obj.N,'full');
            % u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

            x0 = sdpvar(obj.sys.nx,1,'full');

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
                y_t = obj.ccPoly.F*z(:,t) + alpha(t)*obj.ys;
                u_t = repmat(v(:,t),1,obj.ccPoly.v) + alpha(t)*obj.us;
                cost = cost + obj.costFunMan.cost_L_k(y_t, u_t, obj.ys, obj.us);
            end
            y_N = obj.ccPoly.F*z(:,obj.N) + alpha(obj.N)*obj.ys;
            u_N = repmat(v(:,obj.N),1,obj.ccPoly.v) + alpha(obj.N)*obj.us;
            cost = cost + obj.costFunMan.cost_L_N(y_N, u_N, obj.ys, obj.us);

            % % Constraints (look Equation (19) from paper)
            constr = [];

            % Initial condition
            constr = [constr; obj.ccPoly.F*(x0 - z(:,1)) <= alpha(1)*obj.ys];

            % Stage + terminal constraints
            for k=1:obj.N-1
                if var_convh
                    constr = [constr;
                        obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)}, A_sdp,B_sdp)];
                else
                    constr = [constr;
                        obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)} )];
                end
            end
            if var_convh
                constr = [constr;
                    obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1}, A_sdp,B_sdp)];
            else
                constr = [constr;
                    obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1} )];
            end

            % instantiate the Optimizer
            if var_convh
                ocpOptim = optimizer(constr,cost,...
                    obj.sdp_opts,...
                    [{x0},A_sdp(:)',B_sdp(:)'],...% input parameters
                    {z,alpha,v,cost});          % output parameters
            else
                ocpOptim = optimizer(constr,cost,...
                    obj.sdp_opts,...
                    {x0},...          % input parameters
                    {z,alpha,v,cost});  % output parameters
            end
        end


        function dirFeasRegion = initDirFeasRegion(obj, var_convh)
            % define the YALMIP optimizer for maximization of the feasible
            % region with respect to a specific facet
            c_vec = sdpvar(obj.sys.nx,1,'full');
            x0 = sdpvar(obj.sys.nx,1,'full');

            % cost function
            cost = -c_vec'*x0; % maximize the direction

            z = sdpvar(obj.sys.nx, obj.N,'full');
            alpha = sdpvar(1, obj.N,'full');
            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            v = sdpvar(obj.sys.nu, obj.N,'full');
            % u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

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
            constr = [constr; obj.ccPoly.F*(x0 - z(:,1)) <= alpha(1)*obj.ys];

             % Stage + terminal constraints
            for k=1:obj.N-1
                if var_convh
                    constr = [constr;
                        obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)}, A_sdp,B_sdp)];
                else
                    constr = [constr;
                        obj.constrSetS_H_SP({z(:,k),alpha(k)}, v(:,k), {z(:,k+1),alpha(k+1)} )];
                end
            end
            if var_convh
                constr = [constr;
                    obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1}, A_sdp,B_sdp)];
            else
                constr = [constr;
                    obj.constrSetS_H_SP({z(:,obj.N),alpha(obj.N)}, v(:,obj.N), {obj.gamma*z(:,obj.N),obj.gamma*(alpha(obj.N)-1)+1} )];
            end

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


        function feasX0 = initFeasibleX0(obj)
            x0=sdpvar(obj.sys.nx, 1,'full');
            z = sdpvar(obj.sys.nx, obj.N,'full');
            alpha = sdpvar(1, obj.N,'full');
            % y=sdpvar(obj.ccPoly.f, obj.N,'full');
            v = sdpvar(obj.sys.nu, obj.N,'full');
            % u=sdpvar(obj.sys.nu, obj.ccPoly.v, obj.N,'full');

            % input variables
            x0_des = sdpvar(obj.sys.nx, 1,'full');
            costMatrix = sdpvar(obj.sys.nx, obj.sys.nx,'full');

            % cost function
            cost = weighted2NormSquared(x0 - x0_des, costMatrix);

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

            % instantiate the Optimizer
            feasX0 = optimizer(constr, cost,...
                obj.sdp_opts,...
                {x0_des, costMatrix},... % input parameters
                {x0, cost});      % output parameters
        end



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

            constr = [constr; obj.sys.X.A*z + alpha*obj.hx_m <= obj.sys.X.b;
                            obj.sys.U.A*v + alpha*obj.hu_m <= obj.sys.U.b];
            % for j=1:obj.ccPoly.v
            %     constr = [constr;   obj.sys.X.A*(z + alpha*obj.ccPoly.Vi_s{j}*obj.y_m) <= obj.sys.X.b;
            %         obj.sys.U.A*(v + alpha*obj.u_m(:,j)) <= obj.sys.U.b];
            % end
            constr = [constr; alpha >= 0];
        end


        function precomputeSupportFun(obj)
            % use MPT3 methods

            % disturbance vector
            obj.y_bar_s = cell(1,obj.sys.nm);

            for i=1:obj.sys.nm
                Xs_p_verts = zeros(obj.sys.nx,obj.ccPoly.v);
                for j=1:obj.ccPoly.v
                    Xs_p_verts(:,j) = obj.sys.A_convh{i}*obj.ccPoly.Vi_s{j}*obj.ys + obj.sys.B_convh{i}*obj.us(:,j);
                end
                Xs_plus = Polyhedron(Xs_p_verts');
                obj.y_bar_s{i} = Xs_plus.support(obj.ccPoly.F');
            end

            % state and input constraints            
            HX = obj.sys.X.A; HU = obj.sys.U.A;
            % TODO: improve the below initial method!
            obj.hx_m = HX*obj.sys.X.chebyCenter.x;
            obj.hu_m = HU*obj.sys.U.chebyCenter.x;

            for j=1:obj.ccPoly.v
                obj.hx_m = max(obj.hx_m, HX*obj.ccPoly.Vi_s{j}*obj.ys);
                obj.hu_m = max(obj.hu_m, HU*obj.us(:,j));
            end
            


        end

    end

end