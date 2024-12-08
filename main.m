%%PSC for 6dof rigid body using DQ 
%%Rotation sequence - 3,2,1
%% Using Lagrange ploynomial to approximate states and control are approximated by Lagrange ploynomial
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

addpath("RK/","dq_fcn/") 
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

 
q_b_i_fin = quaternion(1,[0 0 0]); % initial quaternion of body frame w.r.t inertial frame
q_b_i_conj_fin = conj(q_b_i_fin);

q_dq_b_i_fin = dualquaternion();
q_dq_b_i_fin.qr = q_b_i_fin; 
r_b_fin = quaternion(1,[0 0 10]);% final translation distance along z axis is 10m
q_dq_b_i_fin.qd =  0.5*q_dq_b_i_fin.qr*r_b_fin; % rotation first followed by translation

w_dq_b_i_b_dot_fin =dualquaternion();
q_dq_b_i_dot_fin = dualquaternion();

% boundary condition

bc = [w_dq_b_i_b.qr.s  (w_dq_b_i_b.qr.v)' w_dq_b_i_b.qd.s (w_dq_b_i_b.qd.v)' q_dq_b_i.qr.s  (q_dq_b_i.qr.v)' q_dq_b_i.qd.s (q_dq_b_i.qd.v)'  w_dq_b_i_b_fin.qr.s  (w_dq_b_i_b_fin.qr.v)' w_dq_b_i_b_fin.qd.s (w_dq_b_i_b_fin.qd.v)' q_dq_b_i_fin.qr.s  (q_dq_b_i_fin.qr.v)' q_dq_b_i_fin.qd.s (q_dq_b_i_fin.qd.v)' w_dq_b_i_b_dot.qr.s  (w_dq_b_i_b_dot.qr.v)' w_dq_b_i_b_dot.qd.s (w_dq_b_i_b_dot.qd.v)' q_dq_b_i_dot.qr.s  (q_dq_b_i_dot.qr.v)' q_dq_b_i_dot.qd.s (q_dq_b_i_dot.qd.v)'  w_dq_b_i_b_dot_fin.qr.s  (w_dq_b_i_b_dot_fin.qr.v)' w_dq_b_i_b_dot_fin.qd.s (w_dq_b_i_b_dot_fin.qd.v)' q_dq_b_i_dot_fin.qr.s  (q_dq_b_i_dot_fin.qr.v)' q_dq_b_i_dot_fin.qd.s (q_dq_b_i_fin_dot.qd.v)'];


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

%% bounds
lb = [xmin.*ones(N+1, 1); xdotmin.*ones(N+1, 1); ymin.*ones(N+1, 1); ydotmin.*ones(N+1, 1); zmin.*ones(N+1, 1); zdotmin.*ones(N+1, 1); Txmin.*ones(N+1, 1); Tymin.*ones(N+1, 1); Tzmin.*ones(N+1, 1)];
ub = [xmax.*ones(N+1, 1); xdotmax.*ones(N+1, 1); ymax.*ones(N+1, 1); ydotmax.*ones(N+1, 1); zmax.*ones(N+1, 1); zdotmax.*ones(N+1, 1); Txmax.*ones(N+1, 1); Tymax.*ones(N+1, 1); Tzmax.*ones(N+1, 1)];

%% equality constraints and linear constraints
A = []; B = []; Aeq = []; Beq = [];

%% initial guess for decision vector - [X(0)...X(N) Y(0)...Y(N) Z(0)...Z(N) 
% Xdot(0)...X(N) Ydot(0)...Ydot(N) Zdot(0)...Zdot(N) 
% Tx(0)...Tx(N) Ty(0)...Ty(N) Tz(0)...Tz(N)]

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