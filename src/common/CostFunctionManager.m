classdef CostFunctionManager
    % Collector class to handle objective function definition, that must be
    % later passed to the OCP classes.
    
    properties (SetAccess = private)
        cost_RCI % must be a convex function (not necessarily strict)
        cost_L_k % strictly convex function (for quadratic cost, Q > 0)
        cost_L_N % strictly convex function (for quadratic cost, P > 0 && Q <= (1-gamma^2)*P)

        LQR_cost % cell array containing (Q,R) for computing the state feedback gain
    end
    
    methods (Access = public)
        function obj = CostFunctionManager(ccPoly, Stage_cost, Term_cost, RCI_cost)
            % if matrices are passed, construct a quadratic cost function;
            % otherwise, set the function handle accordingly.
            % NOTE: in the second case, function signature must be:
            % Stage_cost:   @(y_k,u_k,ys,us)
            % Term_cost:    @(y_N,u_N,ys,us)
            % RCI_cost:     @(ys,us)
                        
            switch class(Stage_cost)
                case 'double'
                    obj.cost_L_k = @(y_k,u_k,ys,us) weighted2NormSquared([y_k-ys; u_k(:)-us(:)],Stage_cost);
                case 'function_handle'
                    obj.cost_L_k = Stage_cost;
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(Stage_cost))
            end
            
            switch class(Term_cost)
                case 'double'
                    obj.cost_L_N = @(y_N,u_N,ys,us) weighted2NormSquared([y_N-ys; u_N(:)-us(:)],Term_cost);
                case 'function_handle'
                    obj.cost_L_N = Term_cost;
                otherwise
                    error("Unexpected DataType: %s. Only 'double' and 'function_handle' accepted.",class(RCI_cost))
            end

            switch class(RCI_cost)
                case 'struct'
                    obj.cost_RCI = obj.buildRCICost(ccPoly, RCI_cost);
                    nx = ccPoly.sys.nx;
                    if isfield(RCI_cost,"Qy")
                        obj.LQR_cost = {RCI_cost.Qy*eye(nx),RCI_cost.Qu*eye(ccPoly.sys.nu)};
                    else
                        obj.LQR_cost = {RCI_cost.Qv(1:nx,1:nx),RCI_cost.Qv(nx+1:end,nx+1:end)};
                    end
                case 'function_handle'
                    obj.cost_RCI = RCI_cost;
                    obj.LQR_cost = {eye(ccPoly.sys.nx),eye(ccPoly.sys.nu)};
                otherwise
                    error("Unexpected DataType: %s. Only 'struct' and 'function_handle' accepted.",class(RCI_cost))
            end
        end
    end
    
    methods (Access = private)

        function cost_RCI = buildRCICost(~,ccpoly,RCI_cost)
            % refer to Section II-D from the paper.

            if isfield(RCI_cost,"Qy") % we're defining the cost in parameter space
                cost_RCI = @(ys,us) 0.5*(ys'*RCI_cost.Qy*ys + us(:)'*RCI_cost.Qu*us(:));
            
            else % the cost is in state space
                V_bar = mean(cat(3, ccpoly.Vi_s{:}), 3); % each Vi_s{i} is a 2D matrix, averaging along 3rd dim.
                U_bar = mean(cat(3, ccpoly.Ui_s{:}), 3); % nu x (nu*v)

                % set volume
                L1 = @(ys, us) sum(arrayfun(@(i) weighted2NormSquared(  [(V_bar-ccpoly.Vi_s{i})*ys ; (U_bar-ccpoly.Ui_s{i})*us(:)], RCI_cost.Qv), 1:ccpoly.v ));
                % set center from zero
                L2 = @(ys, us) weighted2NormSquared( blkdiag(V_bar,U_bar)*[ys;us(:)], RCI_cost.Qc);

                cost_RCI = @(ys, us) L1(ys, us) + L2(ys, us);
            end
        end
    end
end













