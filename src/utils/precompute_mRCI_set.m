function [ys,us] = precompute_mRCI_set(ccPoly,cost_RCI)
opti = casadi.Opti('conic');

% Define decision variables with CasADi
y_m = opti.variable(ccPoly.f, 1);
u_m = opti.variable(ccPoly.sys.nu, ccPoly.v);

% Define cost function
cost = cost_RCI(y_m, u_m);

% Build constraints
opti.subject_to(constrSetS_csd(ccPoly, y_m, u_m, y_m));

% Define the optimization problem
opti.minimize(cost);

% it's often convenient use gurobi when #vertices is large
opti.solver('gurobi',struct(),...
    struct('outputflag',1))
% opti.solver('daqp',struct(),struct('iter_limit',3000))

sol = opti.solve();
ys = sol.value(y_m); us = sol.value(u_m);
end


function constr = constrSetS_csd(ccPoly, y, u, yp)
% This function replicates the Configuration-Constrained RFIT set constraints in CasADi.
% It returns a cell array of constraint expressions.
constr = {};
for j = 1:ccPoly.v   % loop over vertices
    for i = 1:ccPoly.sys.nm   % loop over models
        constr{end+1} = ccPoly.F * (ccPoly.sys.A_convh{i} * ccPoly.Vi_s{j} * y + ...
            ccPoly.sys.B_convh{i} * u(:,j)) + ccPoly.d <= yp;
    end
    constr{end+1} = ccPoly.sys.X.A * ccPoly.Vi_s{j} * y <= ccPoly.sys.X.b;
    constr{end+1} = ccPoly.sys.U.A * u(:,j) <= ccPoly.sys.U.b;
end
constr{end+1} = ccPoly.E * y <= 0;
end