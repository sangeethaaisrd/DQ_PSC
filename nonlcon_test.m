clc
clear all
close all

DV= ones(1,80);
N = (length(DV)-16)/16; % in this case N =4 
m = 20;
J = [3.0514 0 0;0 2.6628 0; 0 0 2.1879]; % Inertia in Kgm^2
% Dual Inertia matrix 8x8 
J_dq = dualInertia(m,J);

w_dq_b_i_b.qr.s = DV(1:N+1);%1-5
w_dq_b_i_b.qr.v(1,:) = DV(N+2:2*N+2); %6-10
w_dq_b_i_b.qr.v(2,:) = DV(2*N+3:3*N+3);%11-15
w_dq_b_i_b.qr.v(3,:) = DV(3*N+4:4*N+4);%16-20

w_dq_b_i_b.qd.s = DV(5*N+1:6*N+1);%21-25
w_dq_b_i_b.qd.v(1,:) = DV(6*N+2:7*N+2);%26-30
w_dq_b_i_b.qd.v(2,:) = DV(7*N+3:8*N+3);%31-35
w_dq_b_i_b.qd.v(3,:) = DV(9*N:10*N);%36-40


q_dq_b_i.qr.s = DV(10*N+1:11*N+1);%41-45
q_dq_b_i.qr.v(1,:) = DV(11*N+2:12*N+2);%46-50
q_dq_b_i.qr.v(2,:) = DV(12*N+3:13*N+3);%51-55
q_dq_b_i.qr.v(3,:) = DV(14*N:15*N) ; %56-60 

q_dq_b_i.qd.s = DV(15*N+1:16*N+1);%61-65
q_dq_b_i.qd.v(1,:) = DV(16*N+2:17*N+2);%66-70
q_dq_b_i.qd.v(2,:) = DV(17*N+3:18*N+3);%71-75
q_dq_b_i.qd.v(3,:) = DV(19*N:20*N);%76-80

F_x = 0;
F_y = 0.01;
F_z = 0.001;
T_x = 0.001;
T_y = 0.1;
T_z = 1;

F = [F_x  F_y  F_z ];  %  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
T = [T_x  T_y  T_z];%  angular velocity of body frame(b) w.r.t inertial frame expressed in body frame(b)
F_q = quaternion(0,F);
T_q = quaternion(0,T);
F_dq = dualquaternion(F_q,T_q)

q_dq_b_i= dq_mat(q_dq_b_i)
w_dq_b_i_b = dq_mat(w_dq_b_i_b)
pdt = q_dq_b_i.*w_dq_b_i_b
% A = J_dq.*w_dq_b_i_b
A = cross_times((w_dq_b_i_b),J_dq.*(w_dq_b_i_b))
B = repmat(F_dq,[1,5])
D = inv(J_dq).*(minus_times(B,A))


function dq = dq_mat(dq_struct)

   n = size(dq_struct.qr.s)
   
   for i=1:1:n(2)
            qr(i) = quaternion(dq_struct.qr.s(:,i),[dq_struct.qr.v(:,i)]);
            qd(i) = quaternion(dq_struct.qd.s(:,i),[dq_struct.qd.v(:,i)]); 
            dq(1,i) = dualquaternion(qr(i),qd(i));
   end
end


function J_dq = dualInertia(m,J)
% Dual inertia matrix. Assumes quaternion of the form [s;v]
    J_dq = [zeros(4,4) [1 zeros(1,3); zeros(3,1) m*eye(3)]; ...
           [1 zeros(1,3); zeros(3,1) J] zeros(4,4)];
end