%% 6DOF rigid body dynamics equation in DQ form 
function[c,ceq] = nonlcon(DV,J_dq,F)


N = (length(DV)-16)/16; % in this case N =100 

w_dq_b_i_b.qr.s = DV(1:N+1);
(w_dq_b_i_b.qr.v) = DV(N+2:4*N+4);
w_dq_b_i_b.qd.s = DV(4*N+5:5*N+5);
(w_dq_b_i_b.qd.v)' = DV(5*N+6:8*N+8);
q_dq_b_i.qr.s = DV(8*N+9:9*N+9);
(q_dq_b_i.qr.v)' = DV(9*N+10:12*N+12);
q_dq_b_i.qd.s = DV(12*N+13:13*N+13);
(q_dq_b_i.qd.v)' = DV(13*N+14:16*N+16);

q_dq_b_i.qr = quaternion(q_dq_b_i.qr.s,(q_dq_b_i.qr.v)');
q_dq_b_i.qd = quaternion(q_dq_b_i.qd.s,(q_dq_b_i.qd.v)');
q_dq_b_i = dualquaternion(q_dq_b_i.qr,q_dq_b_i.qd);


w_dq_b_i_b.qr = quaternion(w_dq_b_i_b.qr.s,(w_dq_b_i_b.qr.v)');
w_dq_b_i_b.qd = quaternion(w_dq_b_i_b.qd.s,(w_dq_b_i_b.qd.v)');
w_dq_b_i_b = dualquaternion(w_dq_b_i_b.qr,w_dq_b_i_b.qd);

A = cross((w_dq_b_i_b),J_dq*(w_dq_b_i_b));

q_dq_b_i_dot_con = 2/(tau_f-tau0)).(D*q_dq_b_i-0.5*q_dq_b_i*w_dq_b_i_b);
w_dq_b_i_b_dot_con = (2/(tau_f-tau0)).*(D*w_dq_b_i_b - inv(J_dq)*(F - A));

ceq = [q_dq_b_i_dot_con;w_dq_b_i_b_dot_con;X(1)-bc(1);Xdot(1)-bc(2);Y(1)-bc(3);Ydot(1)-bc(4);Z(1)-bc(5);Zdot(1)-bc(6)];
c = [abs(X(end))-bc(7); abs(Xdot(end))-bc(8);abs(Y(end))-bc(9);abs(Ydot(end))-bc(10);-(Z(end)-bc(11));Z(end);abs(Zdot(end))-bc(12);fov_con];%X.^2+Y.^2-(Z.^2).*(tan(alpha/2).^2)];%X(ind).^2+Y(ind).^2-(Z(ind).^2).*(tan(alpha/2).^2)];

end