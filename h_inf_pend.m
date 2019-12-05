%% Controle H_inf

% dx = Ax +B_w*w + B_u*u
% z = Cx + D_w*w + D_u*u
%%

clear A B C D

A = [0 0 1 0;
     0 0 0 1;
     0 Mp^2*(Lp/2)^2*Lr*g/Jt -Dr*(Jp+Mp*(Lp/2)^2)/Jt -Mp*(Lp/2)*Lr*Dp/Jt;
     0 Mp*g*(Lp/2)*(Jr+Mp*Lr^2)/Jt -Mp*(Lp/2)*Lr*Dr/Jt -Dp*(Jr+Mp*Lr^2)/Jt]

B = [0; 
     0; 
     (Jp+Mp*(Lp/2)^2)/Jt; 
     Mp*(Lp/2)*Lr/Jt]

C = [1 0 0 0
     0 1 0 0]
 
D_u = zeros(2,1)



%Add actuator dynamics
A(3,3) = A(3,3) - Kg^2*kt*km/Rm*B(3);
A(4,3) = A(4,3) - Kg^2*kt*km/Rm*B(4);
B_u = Kg * kt * B / Rm;

B_w = zeros(4,1);

%% Exemplo apostila robusto
% A = [-0.2  2   2    0;
%       0   -1   0    0;
%       0    0  -2    0;
%       0    0   0  -10];
%  
% B_w = [ 2  0;
%         0 20;
%         0 30;
%         0  0];
% 
% B_u = [0; -20; -30; -3];
% 
% C = [1 0 0 0;
%      0 0 0 -3];
%  
%  D_u = [0; -1];
%  
%  D_w = [0 0; 0 0];
 
 
% Problema: Como deterinar a matriz de ganhos K que estabiliza o sistema e
% minimiza sua norma H_inf.

% Q: Matriz simétrica, definida positiva
% Y = KQ - para linearização pelo complemento de schur
a_row = size(A, 1);
a_col = size(A, 2);
b_row = size(B_u, 1);
b_col = size(B_u, 2);
c_row = size(C, 1); % output size
bw_col = size(B_w, 2); % input size
D_w = zeros(c_row, b_col)

Q = sdpvar(a_row, a_col,'symmetric'); %% Matriz simétrica, definida positiva da dimensão de A
Y = sdpvar(b_col, b_row, 'full'); %% Y = KQ com as dimensões de B'
gama = sdpvar(1); % gamma da norma H_inf - função objetivo
limit = 1.1;
gama = limit^2;

% Se existem matrizes Q, Y, e um escalar gamma > 0, satisfazendo o seguinte
% problema de otimização
LMIs = set([]);
% H_inf = [(A*Q+Q*A'+B_u*Y+Y'*B_u')   B_w                 (Q*C'+Y'*D_u');
%           B_w'                      -gama*eye(bw_col)    D_w';
%           (C*Q+D_u*Y)               D_w                 -eye(c_row)];


H_inf = [(A*Q+Q*A'-B_u*Y-Y'*B_u')   B_w                 (Q*C'-Y'*D_u');
          B_w'                      -gama*eye(bw_col)    D_w';
          (C*Q-D_u*Y)               D_w                 -eye(c_row)];

LMIs = [LMIs, Q>0, H_inf<0];
obj = gama; % Função objetivo a ser minimizada
obj = [];
opt = sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-8, ...
                'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                'sedumi.stepdif', 2);
sol = optimize(LMIs, obj);

Q = double(Q);
Y = double(Y);
gama = double(gama);

% O sistema é estabilizável pela lei de controle u = Y(Q^-1)x
K = Y/Q % K h_inf
K_Q = [-5.26 28.16 -2.76 3.22]


% A norma H_inf do sistema é dada por |G_wz|_inf < gamma
H_inf = sqrt(gama)

%%
% % Plot
% [reference,t] = gensig('square',5,10,0.1);
% reference = deg2rad(20*reference);
% 
% figure(2)
% sys = ss(A-B_u*K_Q, B_u, C, D_u);
% lsim(sys, reference, t)
% %step(sys2)
% 
% figure(1)
% sys2 = ss(A-B_u*K, B_u, C, D_u);
% lsim(sys2, reference, t)
% %step(sys2)