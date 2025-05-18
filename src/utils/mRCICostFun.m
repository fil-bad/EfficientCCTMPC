function cost_RCI_func = mRCICostFun(ccPoly, RCI_cost)
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

cost_RCI_func = casadi.Function('cost_RCI', {ys, us}, {cost_expr},...
    {'y_m','u_m'},{'\ell(y,u)'});
end