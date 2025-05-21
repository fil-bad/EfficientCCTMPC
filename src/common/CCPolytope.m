classdef CCPolytope
    % Definition of polytopes for a given template matrix. The goal is to
    % get polytopes which are entirely simple for their use in the control
    % scheme. A method from [https://arxiv.org/abs/2309.02384,Section IV-D]
    % is also included, to compute an appropriate transformation matrix
    % such that the resulting polytope is a CCPolytope

    properties (SetAccess = private)
        sys % dynamical system associated to the polytope

        f % #facets
        v % #vertices
        F % fixed polytope template
        E % conic constraint in parameter space (Ey<=0)
        Vi_s % cell array of isolated vertices of the polytope

        Ui_s % cell array of vertex control input selectors
        d % support (polytopic enclosure) of the uncertainty set
    end


    methods (Access = public)
        function obj = CCPolytope(sys, template, templateOpts)

            obj.sys = sys;
            if isa(template,'cell') % passing a Triple {F,E,V}
                obj.F = template{1}; obj.E = template{2}; obj.Vi_s = template{3};
                obj.f = size(obj.F,1); obj.v = numel(obj.Vi_s);

            elseif isa(template,'double') % passing only F
                % templateOpts:
                % - computeRCI (def=false) use NLP to get F an RCI template
                % - justFeasibleRCI (def=false) disable NLP maximization cost
                if (nargin < 3)
                    templateOpts = struct('computeRCI',true,'justFeasibleRCI',false,...
                        'y_sigma', ones(size(template,1),1));
                end

                if templateOpts.computeRCI
                    obj.F = obj.computeRCITemplate(template, templateOpts);
                else
                    obj.F = template;
                end
                obj.f = size(obj.F,1);

                % we require that Poly(F,y_sigma) is simple
                obj.Vi_s = obj.getVi_s(templateOpts.y_sigma);
                obj.v = numel(obj.Vi_s);

                obj.E = obj.computeEMatrix(); 
            end

            % reduce the number of redundant rows
            obj.E = Polyhedron(obj.E,zeros(size(obj.E,1),1)).minHRep().A;
            obj.E = sparse(obj.E); % usually very sparse

            obj.Ui_s = arrayfun(@(i) kron(eye(sys.nu), double(1:obj.v == i)),...
                1:obj.v,"UniformOutput",false);

            % compute the polytopic enclosure of disturbance set.
            W_dist_ext = (obj.sys.Bw*obj.sys.W_dist);
            obj.d = W_dist_ext.support(obj.F');
        end
    end

    methods (Access = private)

        function E = computeEMatrix(obj)
            % conic constraint for the chosen face configuration
            I = eye(obj.f);
            E = zeros(obj.v*obj.f, obj.f);
            for i = 1:obj.v
                E((i-1)*obj.f+1:i*obj.f,:) = obj.F * obj.Vi_s{i} - I;
            end
        end

        function V = getVi_s(obj, y_sigma)
            simplePoly_V = Polyhedron(obj.F, y_sigma).V;
            % ensure minimal V-representation. Note that we use directly 
            % ccdlib to avoid convhulln in minVRep() method. If the polytope
            % is too complex, that method just softlock MATLAB
            [verts_min] = cddmex('reduce_v',struct('V',simplePoly_V));
            verts = num2cell(verts_min.V',1);

            % simplePoly = Polyhedron(obj.F, y_sigma).minVRep; 
            % verts_min = num2cell(simplePoly.V',1);

            % computing Vi_s such that Vi*y0 is an isolated vertex of the CCPolytope
            V = cell(size(verts));

            for i=1:length(verts)
                % check for each individual vertex which constraint is active
                [~, idx_nx_active] = mink(abs(obj.F*verts{i} - y_sigma), obj.sys.nx);
                Vi_mask = false(obj.f,1);
                Vi_mask(idx_nx_active) = true;

                % compute 1_Vi to get Vi = inv(F)*1_Vi
                one_mat = zeros(obj.sys.nx, obj.f);  % nx x f
                one_mat(:, Vi_mask) = eye(obj.sys.nx);

                V{i} = obj.F(Vi_mask,:) \ one_mat;
            end
        end

        function F = computeRCITemplate(obj, F_hat, templOpts)
            % from [EfficientCCTMPC paper, Section IV-B]: solve an NLP to
            % find a linear (invertible) transformation of F_hat, so that
            % F = F_hat*T_inv is such that Fx<=1 defines a RCI set.

            nx = obj.sys.nx; nu = obj.sys.nu;
            HW = obj.sys.W_dist.A; hW = obj.sys.W_dist.b; f_W = size(hW,1);

            % ensure Z is in minimal representation
            Z = Polyhedron(F_hat, templOpts.y_sigma).minHRep();

            f_Z = size(Z.A,1);
            % enforce RCI condition on vertices
            zj_s = Z.V';

            % optimization problem
            opti = casadi.Opti();

            T = opti.variable(nx,nx);
            T_inv = opti.variable(nx,nx);
            uj_s = opti.variable(nu,size(zj_s,2));

            if templOpts.justFeasibleRCI
                costFun = 0; %feasible solution
            else
                HX = obj.sys.X.A; hX = obj.sys.X.b;

                % use strong duality to reduce the number of constraints
                Lam_X = opti.variable(f_Z,size(HX,1));
                epsilon = opti.variable(f_Z,1);

                opti.subject_to( Lam_X(:) >= 0);
                opti.subject_to(epsilon >=0);

                opti.subject_to( Lam_X*hX <= ones(f_Z,1) + epsilon );
                opti.subject_to( Lam_X*HX == F_hat*T_inv);

                costFun = sum(epsilon); % 1-norm
            end

            % RCI condition
            for j=1:size(zj_s,2) % #vertices
                for i=1:obj.sys.nm % #models
                    A_i = obj.sys.A_convh{i}; B_i = obj.sys.B_convh{i};

                    % use strong duality on each vertex
                    Lam_W_k = opti.variable(f_Z,f_W);
                    opti.subject_to( Lam_W_k(:) >= 0);
                    opti.subject_to( Lam_W_k*HW == F_hat*T_inv*obj.sys.Bw);
                    w_hat = Lam_W_k*hW;
                    opti.subject_to(...
                        F_hat*T_inv*(A_i*T*zj_s(:,j)+B_i*uj_s(:,j)) + w_hat <= Z.b)
                end
                opti.subject_to( obj.sys.X.A*T*zj_s(:,j) <= obj.sys.X.b );
                opti.subject_to( obj.sys.U.A*uj_s(:,j) <= obj.sys.U.b );
            end
            opti.subject_to( T*T_inv == eye(nx) );

            opti.minimize(costFun);

            % set print_level=0 and print_time=0 to remove any debug output
            % 0, 3 (no iteration), 5 (all details)
            ipopt_opts = struct('print_level',5,'sb','yes','mu_strategy','adaptive');
            opti_opts = struct('ipopt', ipopt_opts, 'print_time', 1,'expand',true);
            opti.solver('ipopt', opti_opts)

            % identity matrix as initial condition (improve NLP convergence)
            opti.set_initial(T, eye(nx));
            opti.set_initial(T_inv,	eye(nx));

            sol = opti.solve();

            F = F_hat*sol.value(T_inv);
        end

    end

end
