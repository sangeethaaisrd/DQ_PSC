%% 6DOF rigid body dynamics equation in DQ form - constraint 
function[c,ceq] = nonlcon(DV,J_dq,F)


N = (length(DV)-22)/22; % in this case N =100 

w_dq_b_i_b.qr.s = DV(1:N+1);
w_dq_b_i_b.qr.v(1,:) = DV(N+2:2*N+2);
w_dq_b_i_b.qr.v(2,:) = DV(2*N+3:3*N+3);
w_dq_b_i_b.qr.v(3,:) = DV(3*N+4:4*N+4);

w_dq_b_i_b.qd.s = DV(4*N+5:5*N+5);
w_dq_b_i_b.qd.v(1,:) = DV(5*N+6:6*N+6);
w_dq_b_i_b.qd.v(2,:) = DV(6*N+7:7*N+7);
w_dq_b_i_b.qd.v(3,:) = DV(7*N+8:8*N+8);

q_dq_b_i.qr.s = DV(8*N+9:9*N+9);
q_dq_b_i.qr.v(1,:) = DV(9*N+10:10*N+10);
q_dq_b_i.qr.v(2,:) = DV(10*N+11:11*N+11);
q_dq_b_i.qr.v(3,:) = DV(11*N+12:12*N+12);

q_dq_b_i.qd.s = DV(12*N+13:13*N+13);
q_dq_b_i.qd.v(1,:) = DV(13*N+14:14*N+14);
q_dq_b_i.qd.v(2,:) = DV(14*N+15:15*N+15);
q_dq_b_i.qd.v(3,:) = DV(15*N+16:16*N+16);


q_dq_b_i.qr = quaternion(q_dq_b_i.qr.s,(q_dq_b_i.qr.v)');
q_dq_b_i.qd = quaternion(q_dq_b_i.qd.s,(q_dq_b_i.qd.v)');
q_dq_b_i = dualquaternion(q_dq_b_i.qr,q_dq_b_i.qd);


w_dq_b_i_b.qr = quaternion(w_dq_b_i_b.qr.s,(w_dq_b_i_b.qr.v)');
w_dq_b_i_b.qd = quaternion(w_dq_b_i_b.qd.s,(w_dq_b_i_b.qd.v)');
w_dq_b_i_b = dualquaternion(w_dq_b_i_b.qr,w_dq_b_i_b.qd);

A = cross((w_dq_b_i_b),J_dq*(w_dq_b_i_b));

q_dq_b_i_dot_con = (2/(tau_f-tau0)).*(D*q_dq_b_i-0.5*q_dq_b_i*w_dq_b_i_b);
w_dq_b_i_b_dot_con = (2/(tau_f-tau0)).*(D*w_dq_b_i_b - inv(J_dq)*(F - A));


ceq = [q_dq_b_i_dot_con.qr.s;(w_dq_b_i_b_dot_con.qr.v)';q_dq_b_i_dot_con.qd.s,(w_dq_b_i_b_dot_con.qd.v)',w_dq_b_i_b.qr.s(1,1)-bc(1);w_dq_b_i_b.qr.v(1,1)-bc(2);w_dq_b_i_b.qr.v(2,1)-bc(3);w_dq_b_i_b.qr.v(3,1)-bc(4);q_dq_b_i.qr.s(1,1)-bc(5);q_dq_b_i.qr.v(1,1)-bc(6);q_dq_b_i.qr.v(2,1)-bc(7);q_dq_b_i.qr.v(3,1)-bc(8)];
c = [w_dq_b_i_b.qr.s(1,end)-bc(9); w_dq_b_i_b.qr.v(1,end)-bc(10);w_dq_b_i_b.qr.v(2,end)-bc(11);w_dq_b_i_b.qr.v(3,end)-bc(12);w_dq_b_i_b.qd.s(1,end)-bc(13);w_dq_b_i_b.qd.v(1,end)-bc(14);w_dq_b_i_b.qd.v(2,end)-bc(15);w_dq_b_i_b.qd.v(3,end)-bc(16)];
end