Lr = [0.15 0.2159];
Lp = [0.3365 0.60];

clear A B C D

for i=1:2
    for j=1:2
        A{2*(i-1)+j} = [0 0 1 0;
                        0 0 0 1;
                        0  Mp^2*(Lp(j)/2)^2*Lr(i)*g/Jt -Dr*(Jp+Mp*(Lp(j)/2)^2)/Jt -Mp*(Lp(j)/2)*Lr(i)*Dp/Jt;
                        0  Mp*g*(Lp(j)/2)*(Jr+Mp*Lr(i)^2)/Jt -Mp*(Lp(j)/2)*Lr(i)*Dr/Jt -Dp*(Jr+Mp*Lr(i)^2)/Jt];
                    
        B{2*(i-1)+j} = [0; 0; (Jp+Mp*(Lp(j)/2)^2)/Jt; Mp*(Lp(j)/2)*Lr(i)/Jt];
        
        B_u{2*(i-1)+j} = Kg * kt * B{2*(i-1)+j} / Rm;
        
        A{2*(i-1)+j}(3,3) = A{2*(i-1)+j}(3,3) - Kg^2*kt*km/Rm*B{2*(i-1)+j}(3);
        A{2*(i-1)+j}(4,3) = A{2*(i-1)+j}(4,3) - Kg^2*kt*km/Rm*B{2*(i-1)+j}(4);
        
        B_w{2*(i-1)+j} = zeros(4,1);
        
        C{2*(i-1)+j} = [1 0 0 0; 0 1 0 0];
        
        D_u{2*(i-1)+j} = zeros(2,1);
        
        c_row = size(C{2*(i-1)+j}, 1);
        b_col = size(B_u{2*(i-1)+j}, 2);
        D_w{2*(i-1)+j} = zeros(c_row, b_col);
    end
end


%%
a_row = size(A{1,1}, 1);
a_col = size(A{1,1}, 2);
b_row = size(B_u{1,1}, 1);
b_col = size(B_u{1,1}, 2);
c_row = size(C{1,1}, 1);
c_col = size(C{1,1}, 2);
bw_col = size(B_w{1,1}, 2); % input size

Q = sdpvar(a_row, a_col,'symmetric'); %% Matriz simétrica, definida positiva da dimensão de A
Y = sdpvar(b_col, b_row, 'full'); %% Y = KQ com as dimensões de B'
gama = sdpvar(1); % gamma da norma H_inf - função objetivo
limit = 20.0;
gama = limit*limit;

LMIs = set([]);
for i=1:4
    H_inf = [(A{i}*Q+Q*A{i}'-B_u{i}*Y-Y'*B_u{i}')       B_w{i}              (Q*C{i}'-Y'*D_u{i}');
              B_w{i}'                                  -gama*eye(bw_col)     D_w{i}';
              (C{i}*Q-D_u{i}*Y)                         D_w{i}              -eye(c_row)];
    
    LMIs = [LMIs, Q>=0, H_inf<=0];
end

LMIs = [LMIs, Q>=0, H_inf<=0];
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
A_hinf = (A{4}-B_u{4}*K); eig(A_hinf)
%%
% Plot
% [reference,t] = gensig('square',5,10,0.1);
% reference = deg2rad(20*reference);
% 
% figure(1)
% sys = ss(A{4}-B_u{4}*K_Q, B_u{4}, C{4}, D_u{4});
% lsim(sys, reference, t)
% %step(sys2)
% 
% figure(2)
% sys2 = ss(A{4}-B_u{4}*K, B_u{4}, C{4}, D_u{4});
% lsim(sys2, reference, t)
%step(sys2)