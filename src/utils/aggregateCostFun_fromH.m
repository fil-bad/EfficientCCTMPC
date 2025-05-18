function cfm_csd = aggregateCostFun_fromH(ccPoly, mRCI_costFun, Q_H, gamma, y_m,u_m)
% compute the cost functions handles and aggregate them inside a struct. we
% have: cost_RCI, cost_L_k_H, cost_L_N_H, cost_L_k_CC, cost_L_N_CC
cfm_csd = struct;

% RCI cost
cfm_csd.cost_RCI = mRCI_costFun;

% Tracking cost (HTMPC)
z = casadi.SX.sym('z_k',ccPoly.sys.nx,1);
v = casadi.SX.sym('v_k',ccPoly.sys.nu,1);
a = casadi.SX.sym('a_k',1,1);

L_k_H = squared2norm([z; v; a-1],Q_H);
cfm_csd.cost_L_k_H = casadi.Function('stage_cost_H',{z,v,a},{L_k_H},...
    {'z_k','v_k','a_k'},{'L_k_H'}).expand();

L_N_H = squared2norm([z; v; a-1],(1/(1-gamma^2))*Q_H);
cfm_csd.cost_L_N_H = casadi.Function('term_cost_H',{z,v,a},{L_N_H},...
    {'z_N','v_N','a_N'},{'L_N_H'}).expand();

% Tracking cost for CCTMPC
y = casadi.SX.sym('y_k',size(y_m));
u = casadi.SX.sym('u_k',size(u_m));

% adapt cost matrix
I_tilde = kron(ones(ccPoly.v,1),eye(ccPoly.sys.nu));
F_tilde = [blkdiag(ccPoly.F,I_tilde), [y_m;u_m(:)] ];

pinv_F_tilde = pinv(F_tilde);
Q_CC = pinv_F_tilde'*Q_H*pinv_F_tilde;

% P = F_tilde * pinv_F_tilde;      
% I_n = eye(size(P,1));
% Z = 1e-4 * (I_n - P);
% 
% Q_CC = Q_CC + Z;

L_k_CC = squared2norm([y-y_m; u(:)-u_m(:)], Q_CC );
cfm_csd.cost_L_k_CC = casadi.Function('stage_cost_CC',{y,u},{L_k_CC},...
    {'y_k','u_k'},{'L_k_CC'}).expand();

L_N_CC = squared2norm([y-y_m; u(:)-u_m(:)], (1/(1-gamma^2))*Q_CC );
cfm_csd.cost_L_N_CC = casadi.Function('term_cost_CC',{y,u},{L_N_CC},...
    {'y_N','u_N'},{'L_N_CC'}).expand();
end