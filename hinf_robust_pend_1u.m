%%setup_rotpen

Lr = [0.16 0.2159];

clear A B C D

for i=1:2
    A{i} = [0 0 1 0;
        0 0 0 1;
        0  Mp^2*(Lp/2)^2*Lr(i)*g/Jt -Dr*(Jp+Mp*(Lp/2)^2)/Jt -Mp*(Lp/2)*Lr(i)*Dp/Jt;
        0  Mp*g*(Lp/2)*(Jr+Mp*Lr(i)^2)/Jt -Mp*(Lp/2)*Lr(i)*Dr/Jt -Dp*(Jr+Mp*Lr(i)^2)/Jt];
    
    B{i} = [0; 0; (Jp+Mp*(Lp/2)^2)/Jt; Mp*(Lp/2)*Lr(i)/Jt];
    
    
    B_u{i} = Kg * kt * B{i} / Rm;
    
    A{i}(3,3) = A{i}(3,3) - Kg^2*kt*km/Rm*B{i}(3);
    A{i}(4,3) = A{i}(4,3) - Kg^2*kt*km/Rm*B{i}(4);
    
    B_w{i} = zeros(4,1);
    
    C{i} = [1 0 0 0; 0 1 0 0];
    
    D_u{i} = zeros(2,1);
    
    c_row = size(C{i}, 1);
    b_col = size(B_u{i}, 2);
    D_w{i} = zeros(c_row, b_col);
end


%%
a_row = size(A{1,1}, 1);
a_col = size(A{1,1}, 2);
b_row = size(B_u{1,1}, 1);
b_col = size(B_u{1,1}, 2);
c_row = size(C{1,1}, 1);
c_col = size(C{1,1}, 2);
bw_col = size(B_w{1,1}, 2); % input size

Q = sdpvar(a_row, a_col,'symmetric'); %% Positive definite matrix
Y = sdpvar(b_col, b_row, 'full'); %% Y = KQ
W = sdpvar(c_row, c_row, 'symmetric'); %  NxN com N igual ao número de linhas de C

gama = sdpvar(1); % Hinf gamma - objective




% limit = 20.0;
% gama = limit*limit;
epsilon = 0.375;

% If exists a symmetric matrix Q, a matrix P and a scalar gamma > 0 such
% that
LMIs = set([]);
for i=1:2
    
    H2_1 = [ W               (C{i}*Q+D_u{i}*Y);
            (Q*C{i}'+Y'*D_u{i}')    Q];
        
   
    H_inf = [(A{i}*Q+Q*A{i}'-B_u{i}*Y-Y'*B_u{i}')       B_w{i}              (Q*C{i}'-Y'*D_u{i}');
              B_w{i}'                                  -gama*eye(bw_col)     D_w{i}';
              (C{i}*Q-D_u{i}*Y)                         D_w{i}              -eye(c_row)];
    
    %control_limit = [Q Y'; Y 5^2];
          
    LMIs = [LMIs, Q>= epsilon*eye(a_row), H_inf<=0, H2_1 >= 0];
end

obj = gama; % Objective function to be minimized

% obj = [];
% Using the SeDuMi solver
opt = sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-8, ...
                'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                'sedumi.stepdif', 2);
sol = optimize(LMIs, obj);

Q = double(Q);
Y = double(Y);
gama = double(gama);

% then the system is stabilizable by the control law u = Y(Q^-1)x
K = Y/Q % K h_inf
K_Q = [-5.26 28.16 -2.76 3.22];


% and the Hinf norm of the system is |G_wz| = sqrt(gamma) 
H_inf = sqrt(gama)
A_hinf = (A{2}-B_u{2}*K); eig(A_hinf)