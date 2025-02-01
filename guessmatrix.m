%clc
clear all
close all


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

w_x_fin = 0.001;
w_y_fin = 0;
w_z_fin = 0.0012; % orbital velocity in rad/s 
v_x_fin = 0.0;
v_y_fin = 0.0;
v_z_fin = 0.0001;

w_b_i_b_fin = [w_x_fin  w_y_fin  w_z_fin ];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
v_b_i_b_fin = [v_x_fin  v_y_fin  v_z_fin ];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
w_q_b_i_b_fin = quaternion(0,w_b_i_b_fin);
v_q_b_i_b_fin = quaternion(0,v_b_i_b_fin);
w_dq_b_i_b_fin = dualquaternion(w_q_b_i_b_fin,v_q_b_i_b_fin);

 
q_b_i_fin = normalize(quaternion(1,[2 3 4])); % final quaternion of body frame w.r.t inertial frame
q_b_i_conj_fin = conj(q_b_i_fin);

q_dq_b_i_fin = dualquaternion();
q_dq_b_i_fin.qr = q_b_i_fin; 
r_b_fin = quaternion(0,[0 0 10]);% final translation distance along z axis is 10m
q_dq_b_i_fin.qd =  0.5*q_dq_b_i_fin.qr*r_b_fin; % rotation first followed by translation

w_q_b_i_b = quaternion(1,[2 3 4]);
w_q_b_i_b_fin = quaternion(5,[6 7 8]);

i=1;
%%
for t =0:0.01:1

q_dq_b_i_guess(i,1) = dqinterp(t,q_dq_b_i,q_dq_b_i_fin);
% A(i) = dot(dq_interp(i).qr,dq_interp(i).qd)% ensured that dot product of real and dual part of DQ is zero.
% norm_qr(i) = norm(dq_interp(i).qr) % ensured that norm of real quaternion = 1 

w_q_b_i_b_guess(i,1)= quat_interp(w_q_b_i_b,w_q_b_i_b_fin,t);
v_q_b_i_b_guess(i,1) = quat_interp(v_q_b_i_b,v_q_b_i_b_fin,t);

%F
%T

i=i+1;

end

