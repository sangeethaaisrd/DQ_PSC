%%PSC for 6dof rigid body using DQ 
%%Rotation sequence - 3,2,1
%% Using Lagrange ploynomial to approximate states and control 
%% Collocation points are generated as the CGL nodes and clenshaw curtis weights are used for approximating the integration.
clc; clear; close all;
set(0, 'defaultFigureWindowState', 'maximized');
set(0, 'defaultAxesFontSize', 20);
set(0, 'defaultAxesLineWidth', 1.5);
set(groot, 'defaultAxesFontName','Century');
set(groot, 'defaultLegendFontName','Century');
set(groot, 'defaultTextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultColorbarTickLabelinterpreter','latex'); 
set(groot, 'defaultLineLineWidth', 1.5);


%% initial conds and final conds and bounds
tau0 = 0;
tau_f = 650; %10 minutes
w_x_ini = 0;
w_y_ini = 0;
w_z_ini = 0.0012; % orbital velocity in rad/s 
v_x_ini = 0.0;
v_y_ini = 0.0;
v_z_ini = 0.1;
mc = 20;
alpha = pi/3;   %FOV

%initial condition 

w_b_i_b = [w_x_ini w_y_ini w_z_ini];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
v_b_i_b = [v_x_ini v_y_ini v_z_ini];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
w_q_b_i_b = quaternion(0,w_b_i_b);
v_q_b_i_b = quaternion(0,v_b_i_b);
w_dq_b_i_b = dualquaternion(w_q_b_i_b,v_q_b_i_b);

 
q_b_i = quaternion(1,[0 0 0]); % initial quaternion of body frame w.r.t inertial frame
q_b_i_conj = conj(q_b_i);

q_dq_b_i = dualquaternion();
q_dq_b_i.qr = q_b_i; 
r_b = quaternion();
q_dq_b_i.qd =  0.5*q_dq_b_i.qr*r_b; % rotation first followed by translation

w_dq_b_i_b_dot =dualquaternion();
q_dq_b_i_dot = dualquaternion();

%final condition

w_x_fin = 0;
w_y_fin = 0;
w_z_fin = 0.0012; % orbital velocity in rad/s 
v_x_fin = 0.0;
v_y_fin = 0.0;
v_z_fin = 0.1;

w_b_i_b_fin = [w_x_fin  w_y_fin  w_z_fin ];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
v_b_i_b_fin = [v_x_fin  v_y_fin  v_z_fin ];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
w_q_b_i_b_fin = quaternion(0,w_b_i_b_fin);
v_q_b_i_b_fin = quaternion(0,v_b_i_b_fin);
w_dq_b_i_b_fin = dualquaternion(w_q_b_i_b_fin,v_q_b_i_b_fin);

 
q_b_i_fin = quaternion(1,[0 0 0]); % final quaternion of body frame w.r.t inertial frame
q_b_i_conj_fin = conj(q_b_i_fin);

q_dq_b_i_fin = dualquaternion();
q_dq_b_i_fin.qr = q_b_i_fin; 
r_b_fin = quaternion(1,[0 0 10]);% final translation distance along z axis is 10m
q_dq_b_i_fin.qd =  0.5*q_dq_b_i_fin.qr*r_b_fin; % rotation first followed by translation

w_dq_b_i_b_dot_fin =dualquaternion();
q_dq_b_i_dot_fin = dualquaternion();

% boundary condition

bc = [w_dq_b_i_b.qr.s  (w_dq_b_i_b.qr.v)' w_dq_b_i_b.qd.s (w_dq_b_i_b.qd.v)' q_dq_b_i.qr.s  (q_dq_b_i.qr.v)' q_dq_b_i.qd.s (q_dq_b_i.qd.v)'  w_dq_b_i_b_fin.qr.s  (w_dq_b_i_b_fin.qr.v)' w_dq_b_i_b_fin.qd.s (w_dq_b_i_b_fin.qd.v)' q_dq_b_i_fin.qr.s  (q_dq_b_i_fin.qr.v)' q_dq_b_i_fin.qd.s (q_dq_b_i_fin.qd.v)'];

% bounds 

%minima

w_x_min = 0;
w_y_min = 0;
w_z_min = 0; % orbital velocity in rad/s 
v_x_min = 0.0;
v_y_min = 0.0;
v_z_min = 0.1;

w_b_i_b_min = [w_x_min  w_y_min  w_z_min ];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
v_b_i_b_min = [v_x_min  v_y_min  v_z_min ];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
w_q_b_i_b_min = quaternion(0,w_b_i_b_min);
v_q_b_i_b_min = quaternion(0,v_b_i_b_min);
w_dq_b_i_b_min = dualquaternion(w_q_b_i_b_min,v_q_b_i_b_min);

 
q_b_i_min = quaternion(1,[0 0 0]); % final quaternion of body frame w.r.t inertial frame
q_b_i_conj_min = conj(q_b_i_min);

q_dq_b_i_min = dualquaternion();
q_dq_b_i_min.qr = q_b_i_min; 
r_b_min = quaternion(1,[0 0 0]);% final translation distance along z axis is 10m
q_dq_b_i_min.qd =  0.5*q_dq_b_i_min.qr*r_b_min; % rotation first followed by translation

%maxima

w_x_max = 0;
w_y_max = 0;
w_z_max = 0; % orbital velocity in rad/s 
v_x_max = 0.0;
v_y_max = 0.0;
v_z_max = 0.1;

w_b_i_b_max = [w_x_max  w_y_max  w_z_max ];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
v_b_i_b_max = [v_x_max  v_y_max  v_z_max ];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
w_q_b_i_b_max = quaternion(0,w_b_i_b_max);
v_q_b_i_b_max = quaternion(0,v_b_i_b_max);
w_dq_b_i_b_max = dualquaternion(w_q_b_i_b_max,v_q_b_i_b_max);

 
q_b_i_max = quaternion(1,[0 0 0]); % final quaternion of body frame w.r.t inertial frame
q_b_i_conj_max = conj(q_b_i_max);

q_dq_b_i_max = dualquaternion();
q_dq_b_i_max.qr = q_b_i_max; 
r_b_max = quaternion(1,[0 0 0]);% final translation distance along z axis is 10m
q_dq_b_i_max.qd =  0.5*q_dq_b_i_max.qr*r_b_max; % rotation first followed by translation




Txmin = -0.3;
Txmax = 0.3;
Tymin = -0.3;
Tymax = 0.3;
Tzmin = -0.3;
Tzmax = 0.3;

Fxmin = -0.3;
Fxmax = 0.3;
Fymin = -0.3;
Fymax = 0.3;
Fzmin = -0.3;
Fzmax = 0.3;

%% Node distribution, Clenshaw Curtis weights and D matrix
a = -1;
b = 1;
N = 70;
tk = (((b - a) / 2) .* cos(linspace(0, pi, N + 1)) + (b + a) / 2)';
D = -Dmatrix_CGL(tk);
w = flip(cc_quad_weights(N));
tk = flip(tk);
tau = ((tau_f-tau0).*tk+(tau_f+tau0))./2;

%% TODO:bounds

lb = [w_dq_b_i_b_min.qr.s.*ones(N+1, 1);  (w_dq_b_i_b_min.qr.v)'.*ones(N+1, 1);  w_dq_b_i_b_min.qd.s*ones(N+1, 1); (w_dq_b_i_b_min.qd.v)'.*ones(N+1, 1);  q_dq_b_i_min.qr.s*ones(N+1, 1);  (q_dq_b_i_min.qr.v)'.*ones(N+1, 1);  q_dq_b_i_min.qd.s*ones(N+1, 1); (q_dq_b_i_min.qd.v)'.*ones(N+1, 1); ];
ub = [w_dq_b_i_b_max.qr.s*ones(N+1, 1);  (w_dq_b_i_b_max.qr.v)'.*ones(N+1, 1);  w_dq_b_i_b_max.qd.s*ones(N+1, 1); (w_dq_b_i_b_max.qd.v)'.*ones(N+1, 1);  q_dq_b_i_max.qr.s*ones(N+1, 1);  (q_dq_b_i_max.qr.v)'.*ones(N+1, 1);  q_dq_b_i_max.qd.s*ones(N+1, 1); (q_dq_b_i_max.qd.v)'.*ones(N+1, 1); ];

%% equality constraints and linear constraints
A = []; B = []; Aeq = []; Beq = [];
% TODO : Here introdcue the constraint - norm of real quaternion =1 . qd.qr = 0

%% initial guess for decision vector - [qr_0(0)...qr_0(N) qr_1(0)...qr_1(N) qr_2(0)...qr_2(N)...qr_3(0)...qr_3(N) 
% qd_0(0)     ...... Ydot(0)...Ydot(N) Zdot(0)...Zdot(N) 
% Tx(0)...Tx(N) Ty(0)...Ty(N) Tz(0)...Tz(N)]#

DV0 = awgn([xguess; xdotguess; yguess; ydotguess; zguess; zdotguess; Txguess; Tyguess; Tzguess],650);
%DV0 = [xguess; xdotguess; yguess; ydotguess; zguess; zdotguess; Txguess; Tyguess; Tzguess];

%% Optimization options
options =  optimoptions ('fmincon','Display','Iter','OptimalityTolerance',...
1e-4, 'ConstraintTolerance', 1e-1, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
500000,'Algorithm','sqp');

[DV, costval, exitflag, output] = fmincon(@(DV)costfunc(DV, w, tau0, tau_f), DV0, A, B,...
    Aeq, Beq, lb, ub, @(DV)nonlcon(DV, D, bc, w, tau0, tau_f, omega, mc, alpha),options);

exitflag
output