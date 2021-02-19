% Numerical Example Oya2006
clc
clear

% Define A and B matrix
A = [1 4; 0 -1];
B = [0; 1];


P = sdpvar(2,2, 'symmetric');
Q = eye(2);
Y = sdpvar(2,2);
W = sdpvar(1, 1, 'symmetric');
a = 3;

% Equation 7 in LMI form
LMIs = set([]);
riccati_constraint = [Y*A'+A*Y-B*W'*B'-B*W*B'+2*a*Y Y; 
                              Y -1\Q];
                                                    
% Equation 7 = 0
const = [LMIs, riccati_constraint <= 0, Y>=0, W>=0];

% You can choose one of the solvers below
opt = sdpsettings('solver', 'sedumi'); % SeDuMi
%opt = sdpsettings('solver', 'sdpt3'); % SDPT-3

optimize(const,[],opt)

% Get value of X
Y = double(Y)
W = double(W)
P = 1\Y
K = -W*B'*P
meta = [1.0855 -1.5356; -1.5356 4.5318]
w_meta = 0.21708

K = [-0.13017 -9.2008];
Y = [1.0855 -1.5356; -1.5356 4.5318];
Ak = A+B*K;
y = sdpvar(1);
sig = 3.1612;
epsilon = sdpvar(1);

LMIs_2 = set([]);
eq_25 = [Ak'*P+P*A+y*P^2 eye(2); eye(2) (-y/sig)*eye(2)];
eq_26 = [y 1; 1 epsilon];

consts = [LMIs_2, eq_25 <=0, eq_26>=0, y>=0];
obj = epsilon;

sol = optimize(consts, obj, opt);
y = double(y)
epsilon = double(epsilon)
                          
                          
