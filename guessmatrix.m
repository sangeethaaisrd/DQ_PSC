clc
clear all
close all


%% initial conds and final conds and bounds
%% initial conds and final conds and bounds
tau0 = 0;
tau_f = 650;%650; %10 minutes
w_x_ini = 0;
w_y_ini = 0;
w_z_ini = 0.0; % orbital velocity in rad/s 
v_x_ini = 0.0;
v_y_ini = 0.0;
v_z_ini = 0.0;% in m/s
m = 20;

% xini = 0.05;
% xfin = 0.02;
% yini = 0.05;
% yfin = 0.02;
% zini = -10.05;
% zfin = -0.1;

% xdotini = 0.0;
% xdotfin = 0.005;
% ydotini = 0.0;
% ydotfin = 0.005;
% zdotini = 0.0;
% zdotfin = 0.005;
%node distribution 
a = -1;
b = 1;
N = 40;
tk = (((b - a) / 2) .* cos(linspace(0, pi, N + 1)) + (b + a) / 2)';
tau = ((tau_f-tau0).*tk+(tau_f+tau0))./2;


%initial condition 

w_b_i_b_ini = [w_x_ini w_y_ini w_z_ini];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
v_b_i_b_ini = [v_x_ini v_y_ini v_z_ini];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
w_q_b_i_b_ini = quaternion(0,w_b_i_b_ini);
v_q_b_i_b_ini = quaternion(0,v_b_i_b_ini);
w_dq_b_i_b_ini = dualquaternion(w_q_b_i_b_ini,v_q_b_i_b_ini);

 
q_b_i_ini = quaternion(1,[0 0 0]); % initial quaternion of body frame w.r.t inertial frame
% q_b_i_conj = conj(q_b_i);

q_dq_b_i_ini = dualquaternion();
q_dq_b_i_ini.qr = q_b_i_ini; 
r_b_ini = quaternion(0,[0.05,0.05,-10.05]);
q_dq_b_i_ini.qd =  0.5*q_dq_b_i_ini.qr*r_b_ini; % rotation first followed by translation

w_dq_b_i_b_dot_ini =dualquaternion();
q_dq_b_i_dot_ini = dualquaternion();

%final condition

w_x_fin = 0.0;
w_y_fin = 0.0;
w_z_fin = 0.0; % orbital velocity in rad/s 
v_x_fin = 0.005;
v_y_fin = 0.005;
v_z_fin = 0.005;

w_b_i_b_fin = [w_x_fin  w_y_fin  w_z_fin ];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
v_b_i_b_fin = [v_x_fin  v_y_fin  v_z_fin ];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
w_q_b_i_b_fin = quaternion(0,w_b_i_b_fin);
v_q_b_i_b_fin = quaternion(0,v_b_i_b_fin);
w_dq_b_i_b_fin = dualquaternion(w_q_b_i_b_fin,v_q_b_i_b_fin);

 
q_b_i_fin = normalize(quaternion(0.7866,[0.4330 0.0795 -0.4330])); % final quaternion of body frame w.r.t inertial frame(-45,30,45)euler angle
           
% q_b_i_conj_fin = conj(q_b_i_fin);

q_dq_b_i_fin = dualquaternion();
q_dq_b_i_fin.qr = q_b_i_fin; 
r_b_fin = quaternion(0,[0.02 0.02 -0.1]);% final translation distance along z axis is 10m
q_dq_b_i_fin.qd =  0.5*q_dq_b_i_fin.qr*r_b_fin; % rotation first followed by translation


load("guess.mat"); %

Th = resample(Thrust, tau);

i=1;
t0=0;
tf =1;
% N = 20;
incr = (tf-t0)/N;
%%
for t = t0 :incr:tf

q_dq_b_i_guess(i,1) = dqinterp(t,q_dq_b_i_ini,q_dq_b_i_fin); % norm of qr is verified to be one for all guessvalues

% A(i) = dot(dq_interp(i).qr,dq_interp(i).qd)% ensured that dot product of real and dual part of DQ is zero.
% norm_qr(i) = norm(dq_interp(i).qr) % ensured that norm of real quaternion = 1 
r_b(i,1) = rb_from_dq(q_dq_b_i_guess(i,1));
w_b_i_b_guess = interpolate_3d(w_b_i_b_ini,w_b_i_b_fin,t);
v_b_i_b_guess = interpolate_3d(v_b_i_b_ini,v_b_i_b_fin,t);

w_q_b_i_b_guess(i,1) = quaternion(0,w_b_i_b_guess);
v_q_b_i_b_guess(i,1) = quaternion(0,v_b_i_b_guess);

w_dq_b_i_b_guess(i,1) = dualquaternion(w_q_b_i_b_guess(i,1),v_q_b_i_b_guess(i,1));
F_q_guess(i,1) = quaternion(0,[Th.Data(i, 1),Th.Data(i, 2),Th.Data(i, 3)]);
% T_q_guess(i,1) =quaternion(0,[Th.Data(i, 1),Th.Data(i, 2),Th.Data(i, 3)]);
T_q_guess(i,1) = quaternion();
F_dq_guess(i,1) = dualquaternion(F_q_guess(i,1),T_q_guess(i,1));

guess_1(:,i) = dq2vec(w_dq_b_i_b_guess(i,1));
guess_2(:,i) = dq2vec(q_dq_b_i_guess(i,1));
guess_3(:,i) = dq2vec(F_dq_guess(i,1));
%F
%T

i=i+1;

end
size(guess_1)



guess_1 =guess_1(:);
guess_2 =guess_2(:);
guess_3 =guess_3(:);

guess_DV= cat(1,guess_1,guess_2,guess_3);



function interpolated_point = interpolate_3d(point1, point2, t)

    interpolated_point = (1 - t) .* point1 + t .* point2;

end
