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



% %%
% clc
% clear A B C D D_u D_w
% for i=1:2
%     for j=1:2
%         A{2*(i-1)+j} = [0 0 0 0;
%                         0 0 0 0;
%                         2 0 1 2;
%                             0 2 2 1];
% 
%         B_u{2*(i-1)+j} = [1 0; 0 1; 0 0; 0 0];
%         B_w{2*(i-1)+j} = [0 0; 0 0; 0 0; 0 0];
% 
%         C{2*(i-1)+j} = [4 -3 3 0; -2 2 -1 2];
% 
%         D_u{2*(i-1)+j} = zeros(c_row,b_col);
%         D_w = D_u;
%         
%     end
% end
% 
% Y = sdpvar(2, 2, 'full'); %% Y = KQ com as dimensões de B'
% a_row = size(A{1,1}, 1);
% a_col = size(A{1,1}, 2);
% b_row = size(B_u{1,1}, 1);
% b_col = size(B_u{1,1}, 2);
% c_row = size(C{1,1}, 1);
% c_col = size(C{1,1}, 2);
% bw_col = size(B_w{1,1}, 2); % input size

% D_w = zeros(c_row, b_col);
% D_u = zeros(c_row,b_col);

Q = sdpvar(a_row, a_col,'symmetric'); %% Matriz simétrica, definida positiva da dimensão de A
Y = sdpvar(2, 1, 'full'); %% Y = KQ com as dimensões de B'
gama = sdpvar(1); % gamma da norma H_inf - função objetivo

LMIs = set([]);
for i=1:4
    %Y = N
    %W = Q
    phi = Q*A{i}' + A{i}*Q - B_u{i}*Y*C{i} - C{i}'*Y'*B_u{i}';
    a12 = B_w{i};
    a13 = Q*C{i}'-C{i}'*Y'*D_u{i}';
    a21 = B_w{i}';
    a22 = -gama*eye(bw_col);
    a23 = D_w{i}';
    a31 = (C{i}*Q-D_u{i}*Y*C{i});
    a32 = D_w{i};
    a33 = -eye(c_row);
    
    H_inf = [phi    a12    a13;
             a21    a22    a23;
             a31    a32    a33];
         
    LMIs = [LMIs, Q>=0, H_inf <= 0];
end


LMIs = [LMIs, Q>=0, H_inf <= 0];
obj = gama; % Função objetivo a ser minimizada
% obj = [];
opt = sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-8, ...
                'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                'sedumi.stepdif', 2);
sol = optimize(LMIs, gama, opt);

Q = double(Q);
Y = double(Y);
% gama = double(gama);

% O sistema é estabilizável pela lei de controle u = Y(Q^-1)x
M = C{4}*Q/C{4};
K = Y/M % K h_inf
K_Q = [-5.26 28.16 -2.76 3.22];
K = [4.6497 6.2001; 5.7551 8.2946];

% A norma H_inf do sistema é dada por |G_wz|_inf < gamma
H_inf = sqrt(gama)
A_hinf = (A{4}-B_u{4}*K*C{4}); eig(A_hinf)
sys = ss(A_hinf, B_u{4}, C{4}, D_u{4});
step(sys)
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