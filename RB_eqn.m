%% 6DOF rigid body dynamics equation in DQ form 
function[result] = RB_eqn (DV,J_dq,F)


N = (length(DV)-16)/16; % in this case N =100 

w_dq_b_i_b.qr.s = DV(1:N+1);
(w_dq_b_i_b.qr.v)' = DV(N+2:4*N+4);
w_dq_b_i_b.qd.s = DV(4*N+5:5*N+5);
(w_dq_b_i_b.qd.v)' = DV(5*N+6:8*N+8);
q_dq_b_i.qr.s = DV(8*N+9:9*N+9);
(q_dq_b_i.qr.v)' = DV(9*N+10:12*N+12);
q_dq_b_i.qd.s = DV(12*N+13:13*N+13);
(q_dq_b_i.qd.v)' = DV(13*N+14:16*N+16);

w_dq_b_i_b.qr = quaternion(w_dq_b_i_b.qr.s,(w_dq_b_i_b.qr.v)');
w_dq_b_i_b.qd = quaternion(w_dq_b_i_b.qd.s,(w_dq_b_i_b.qd.v)');
w_dq_b_i_b = dualquaternion(w_dq_b_i_b.qr,w_dq_b_i_b.qd);

A = cross((w_dq_b_i_b),J_dq*(w_dq_b_i_b));
result = inv(J_dq)*(F - A);
end