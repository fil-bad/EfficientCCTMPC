classdef CostFunctionManagerCasadi
    
    properties (SetAccess = private)
        cost_RCI    % casadi.Function with signature: (ys, us) -> scalar
        cost_L_k    % casadi.Function with signature: (y_k, u_k, ys, us) -> scalar
        cost_L_N    % casadi.Function with signature: (y_N, u_N, ys, us) -> scalar
        cost_L_k_H  % function for HTMPC cost, signature: (z_k,v_k,a_k) -> scalar
        cost_L_N_H  % function for HTMPC cost, signature: (z_N,v_N,a_N) -> scalar
    end

    methods (Access = public)
        function obj = CostFunctionManagerCasadi(ccPoly, RCI_cost, Stage_cost, Term_cost, mRCIset)
            % In this CasADi version, Stage_cost and Term_cost must be double matrices.
            % RCI_cost is expected to be a struct or a casadi.Function.
            
            if nargin == 5 % compute the 

            end



            if isa(Stage_cost, 'double')
                obj.cost_L_k = obj.buildStageCost(ccPoly, Stage_cost);
            else
                error('Stage_cost must be provided as a double matrix for the CasADi version.');
            end

            if isa(Term_cost, 'double')
                obj.cost_L_N = obj.buildTerminalCost(ccPoly, Term_cost);
            else
                error('Term_cost must be provided as a double matrix for the CasADi version.');
            end

            if isstruct(RCI_cost)
                obj.cost_RCI = obj.buildRCICost(ccPoly, RCI_cost);
            elseif isa(RCI_cost,'casadi.Function')
                obj.cost_RCI = RCI_cost;
            else
                error('Wrong type for RCI_cost. Only struct and CasADi function are supported.');
            end
        end
    end

    methods (Access = private)
        function cost_RCI_func = buildRCICost(~, ccPoly, RCI_cost)
            % Build the RCI cost function as a casadi.Function.
            ys = casadi.SX.sym('ys', ccPoly.f, 1);
            us = casadi.SX.sym('us', ccPoly.sys.nu, ccPoly.v);

            if isfield(RCI_cost, "Qy")
                % Cost defined in the parameter space:
                cost_expr = 0.5 * ( ys' * RCI_cost.Qy * ys + us(:)' * RCI_cost.Qu * us(:) );
            else
                % Cost defined in state space.
                V_bar = mean(cat(3, ccPoly.Vi_s{:}), 3);
                U_bar = mean(cat(3, ccPoly.Ui_s{:}), 3);
                L_vol = 0;
                for i = 1:ccPoly.v
                    diffV = V_bar - ccPoly.Vi_s{i};
                    diffU = U_bar - ccPoly.Ui_s{i};
                    L_vol = L_vol + squared2norm([diffV * ys; diffU*us(:)], RCI_cost.Qv);
                end
                L_cen = squared2norm( blkdiag(V_bar, U_bar) * [ys; us(:)], RCI_cost.Qc );
                cost_expr = L_vol + L_cen;
            end

            cost_RCI_func = casadi.Function('cost_RCI', {ys, us}, {cost_expr}).expand();
        end

        function stage_cost_func = buildStageCost(~, ccPoly, Stage_cost)
            % Build the stage cost function as a casadi.Function.
            % Expected inputs: y_k, u_k, ys, us.
            y_k = casadi.SX.sym('y_k', ccPoly.f, 1);
            u_k = casadi.SX.sym('u_k', ccPoly.sys.nu, ccPoly.v);
            ys = casadi.SX.sym('ys', ccPoly.f, 1);
            us = casadi.SX.sym('us', ccPoly.sys.nu, ccPoly.v);

            expr = squared2norm([y_k - ys; u_k(:) - us(:)], Stage_cost);
            stage_cost_func = casadi.Function('stage_cost', {y_k, u_k, ys, us}, {expr}).expand();
        end

        function term_cost_func = buildTerminalCost(~, ccPoly, Term_cost)
            % Build the terminal cost function as a casadi.Function.
            % Expected inputs: y_N, u_N, ys, us.
            y_N = casadi.SX.sym('y_N', ccPoly.f, 1);
            u_N = casadi.SX.sym('u_N', ccPoly.sys.nu, ccPoly.v);
            ys = casadi.SX.sym('ys', ccPoly.f, 1);
            us = casadi.SX.sym('us', ccPoly.sys.nu, ccPoly.v);

            expr = squared2norm([y_N - ys; u_N(:) - us(:)], Term_cost);
            term_cost_func = casadi.Function('term_cost', {y_N, u_N, ys, us}, {expr}).expand();
        end
    end
end
