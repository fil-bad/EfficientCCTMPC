classdef CCPolytope_test 
    % Definition of polytopes for a given template matrix. The goal is to
    % get polytopes which are entirely simple for their use in the control
    % scheme. A method from [https://arxiv.org/abs/2309.02384,Section IV-D]
    % is also included, to compute an appropriate transformation matrix
    % such that the resulting polytope is a CCPolytope

    properties (SetAccess = private)
        sys % dynamical system associated to the polytope
        ccpolytope % Configuration-Constrained polytope
        f % #facets
        v % #vertices
        Vi_s % cell array of isolated vertices of the polytope
        Ui_s % cell array of vertex control input selectors
        E % conic constraint in parameter space (Ey<=0)

        d % support (polytopic enclosure) of the uncertainty set

        uj_templ
        dual_gap
        y_sigma
        F
    end

    % properties(Dependent)
    %     F % in H-rep (Ax<=b) for the CCPolytope, F=A. Defined for brevity
    % end
    % 
    % methods % GETTER methods
    %     function F = get.F(obj)
    %         F = obj.ccpolytope.A;
    %     end
    % end

    methods (Access = public)
        function obj = CCPolytope_test(sys, F, getCCPoly,justFeasible)
            % Two admissible constructor signatures:
            % {sys, F}: use a provided CC-template directly (useful for nx == 2)
            % {sys, C_tilde, setDist}: compute a valid CC-template (suggested for nx >= 3)

            obj.sys = sys;

            if nargin <= 3
                justFeasible = false;
            end

            if isa(F,"Polyhedron")
                obj.ccpolytope = F;
                % obj.ccpolytope = Polyhedron( F.A./F.b,ones(size(F.b)) ).minHRep().minVRep();
            else
                if getCCPoly % modify template F to get a CCPolytope
                    [F,obj.uj_templ,obj.dual_gap] = obj.getCCTemplateMatrixCASADI(F,justFeasible);
                    % F(abs(F)<5e-8) = 0;
                end
                % using ones() as initial parameter, the CCPolytope is entirely simple
                % obj.ccpolytope = Polyhedron(F,ones(size(F,1),1)).minHRep().minVRep();
            end

            obj.F = F;

            obj.f = size(obj.F,1);
            obj.y_sigma = ones(obj.f,1);

            obj.Vi_s = obj.getVi_s();
            obj.v = numel(obj.Vi_s);

            e = @(i) double(1:obj.v == i);
            obj.Ui_s = arrayfun(@(i) kron(eye(sys.nu),e(i)), 1:obj.v,"UniformOutput",false);


            obj.E = sparse(obj.computeEMatrix()); % usually very sparse
            % obj.E(abs(obj.E)<5e-8) = 0;
            % reduce the number of redundant rows
            obj.E = Polyhedron(obj.E,sparse(size(obj.E,1),1)).minHRep().A;
            obj.E = sparse(obj.E);

            % compute the polytopic enclosure of disturbance set.
            W_dist_ext = (obj.sys.Bw*obj.sys.W_dist);
            obj.d = W_dist_ext.support(obj.F');
        end
    end

    methods (Access = private)

        function E = computeEMatrix(obj)
            % conic constraint defining the vertex configuration domain explicitly
            E = zeros(obj.v*obj.f, obj.f);
            for i=1:obj.v
                E((i-1)*obj.f+1:i*obj.f,:) = obj.F*obj.Vi_s{i}-eye(obj.f);
            end
        end


        function V = getVi_s(obj)
            % computing Vi_s such that Vi*y0 is an isolated vertex of the CCPolytope
            % vertices = num2cell(obj.ccpolytope.V',1);
            vertices = num2cell(obj.findPolytopeVertices(obj.F, obj.y_sigma)',1);
            V = cell(size(vertices));

            for i=1:length(vertices)
                % check for each individual vertex which constraint is active
                % Vi_mask = (obj.F*vertices{i}-obj.ccpolytope.b).^2 <= 1e-12;

                [~,idx_active] = sort(abs(obj.F*vertices{i}-obj.y_sigma));
                Vi_mask = false(obj.f,1);
                Vi_mask(idx_active(1:obj.sys.nx)) = true(obj.sys.nx,1);

                % Vi_mask = abs(obj.F*vertices{i}-obj.y_sigma) <= 1e-13;
                % assert(nnz(Vi_mask)==obj.sys.nx, "The vertex cannot be defined uniquely")

                % compute 1_Vi to get Vi = inv(F)*1_Vi
                one_mat = zeros(obj.sys.nx,obj.f); % nx x f
                one_mat(sub2ind(size(one_mat),1:obj.sys.nx,find(Vi_mask)')) = 1;

                V{i} = obj.F(Vi_mask,:) \ one_mat;
            end
        end

        function vertices = findPolytopeVertices(obj, A, b)
% findPolytopeVertices - Computes all vertices of a polytope given by Ax <= b.
%
% Syntax: vertices = findPolytopeVertices(A, b)
%
% Inputs:
%    A - m-by-n matrix of coefficients in the linear inequalities.
%    b - m-by-1 vector.
%
% Outputs:
%    vertices - p-by-n matrix. Each row corresponds to a vertex of the polytope.
%
% Example:
%    % Define a square in R^2:
%    A = [ 1  0;
%          0  1;
%         -1  0;
%          0 -1];
%    b = [1; 1; 1; 1];
%    verts = findPolytopeVertices(A, b);
%    disp(verts);
%
% Note:
%    It is assumed that the polytope is 'simple' so that each vertex is the 
%    intersection of exactly n constraints.
%
% Author: Your Name
% Date: Today's Date

    tol = 1e-11;   % Tolerance for numerical comparisons

    [m, n] = size(A);
    vertices = [];

    % Generate all combinations of m constraints taken n at a time.
    combs = nchoosek(1:m, n);
    num_combs = size(combs, 1);

    for idx = 1:num_combs
        I = combs(idx, :);
        A_active = A(I, :);
        b_active = b(I);
        
        % Only consider combinations where the active constraints define a unique point.
        if rank(A_active) < n
            continue;
        end
        
        % Solve the system A_active * x = b_active
        x_candidate = A_active \ b_active;
        
        % Check feasibility: the candidate must satisfy all the inequalities.
        if all(A * x_candidate <= b + tol)
            vertices = [vertices; x_candidate(:)']; %#ok<AGROW>
        end
    end

    % Remove duplicate vertices (accounting for numerical accuracy issues)
    vertices = unique(round(vertices/tol)*tol, 'rows');
end



        function [F,uj_templ,dual_gap] = getCCTemplateMatrixCASADI(obj, F_hat, justFeasible)
            % from [https://doi.org/10.1109/TAC.2024.3454528,Section III-D]:
            % - firstly, define a C_tilde matrix (f,nx), thus fixing the
            %   polytope complexity)
            % - then, solve an NLP to find the proper transformation
            %   matrix W_inv such that F = C_tilde*W_inv is a template for
            %   a CCPolytope.

            % here we simplify the problem as the minkowski sum of two
            % polytopes with the same template is the sum of their RHSes.

            nx = obj.sys.nx; nu = obj.sys.nu; nw = size(obj.sys.W_dist.b,1);

            HW = obj.sys.W_dist.A; hW = obj.sys.W_dist.b;
            

            Z_set = Polyhedron(F_hat, ones(size(F_hat,1),1)).minHRep();

            F_hat_f = size(Z_set.A,1);
            % we enforce RCI condition on these vertices zj_s
            zj_s = num2cell(Z_set.V',1);

            % optimization problem
            opti = casadi.Opti();

            W = opti.variable(nx,nx);
            M = opti.variable(nx,nx);
            uj_s = opti.variable(nu,length(zj_s));

            if justFeasible
                costFun = 0; %feasible solution
            else
                % % % % maximize volume of the Affine (Linear) Transformation
                % W_sx = casadi.SX.sym('W',nx, nx);
                % [~, R] = qr(W_sx);
                % % qr_fun = casadi.Function('qr_fun',{W_sx},{trace(log(R))});
                % qr_fun = casadi.Function('qr_fun',{W_sx},{diag(log(R))'*[100*ones(3,1);1*ones(7,1)]});
                % costFun = -qr_fun(W);


                HX = obj.sys.X.A; hX = obj.sys.X.b;

                epsilon = opti.variable(F_hat_f,1);
                opti.subject_to(epsilon >=0);

                % costFun = epsilon'*epsilon; % 2-norm
                costFun = sum(epsilon); % 1-norm
                % using strong duality to reduce the number of constraints
                Lam_X = opti.variable(F_hat_f,size(HX,1));
                opti.subject_to( Lam_X(:) >= 0);

                opti.subject_to( Lam_X*hX <= ones(F_hat_f,1) + epsilon );
                opti.subject_to( Lam_X*HX == F_hat*M);
            end

            % using strong duality to reduce the number of constraints
            w_hat = opti.variable(F_hat_f,1);
            Lambda = opti.variable(F_hat_f,nw);

            opti.subject_to( w_hat >= 0);
            opti.subject_to( Lambda(:) >= 0);

            opti.subject_to( Lambda*hW <= w_hat);
            % w_hat = Lambda*hW;
            opti.subject_to( Lambda*HW == F_hat*M*obj.sys.Bw);

            constr_dual = cell(obj.sys.nm, length(zj_s));
            for j=1:length(zj_s) % vertices
                for i=1:obj.sys.nm % #models
                    constr_dual{i,j} =  F_hat*M*(obj.sys.A_convh{i}*W*zj_s{j}+obj.sys.B_convh{i}*uj_s(:,j)) + w_hat <=ones(size(F_hat,1),1);
                    opti.subject_to(constr_dual{i,j});
                end
                opti.subject_to( obj.sys.X.A*W*zj_s{j} <= obj.sys.X.b );
                opti.subject_to( obj.sys.U.A*uj_s(:,j) <= obj.sys.U.b );
            end
            opti.subject_to( W*M == eye(nx) );


            opti.minimize(costFun);

            % set print_level=0 and print_time=0 to remove any debug output
            % 0, 3 (no iteration), 5 (all details)
            ipopt = struct('print_level',5,'sb','yes');%, 'mu_strategy', 'adaptive');
            opti_opts = struct('ipopt', ipopt, 'print_time', 1, 'expand',true);
            opti.solver('ipopt',opti_opts)

            % identity matrix as initial condition (improve NLP convergence)
            opti.set_initial(W, eye(nx));
            opti.set_initial(M,	eye(nx));

            % opti.set_initial(epsilon,1);

            sol = opti.solve();

            F = F_hat*sol.value(M);

            uj_templ = sol.value(uj_s);


            dual_gap = cell(obj.sys.nm, length(zj_s));
            for j=1:length(zj_s) % vertices
                for i=1:obj.sys.nm % #models
                    dual_gap{i,j} = sol.value(opti.dual(constr_dual{i,j}));
                end
            end

            % disp("cost: ")
            % disp(sol.value(W))
            
            % disp("cost: " + sum(sol.value(epsilon)))
            % nlp_M = opti.to_function('nlp',{},{M},{},{'M'});
            % M_struct = nlp_M();
            % F = F_hat*full(M_struct.M);
        end

    end

end