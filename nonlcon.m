 %% 6DOF rigid body dynamics equation in DQ form - constraint 
function[c,ceq] = nonlcon(DV,J_dq)

%DV= ones(1,2424);
N = (length(DV)-24)/24 % in this case N =100 

w_dq_b_i_b.qr.s = DV(1:N+1);%1-101 % This will be all zero always 
w_dq_b_i_b.qr.v(1,:) = DV(N+2:2*N+2); %102-202
w_dq_b_i_b.qr.v(2,:) = DV(2*N+3:3*N+3);%203-303
w_dq_b_i_b.qr.v(3,:) = DV(3*N+4:4*N+4);%303-404

w_dq_b_i_b.qd.s = DV(4*N+5:5*N+5);%405-505 % This will be all zero always 
w_dq_b_i_b.qd.v(1,:) = DV(5*N+6:6*N+6);%506-606
w_dq_b_i_b.qd.v(2,:) = DV(6*N+7:7*N+7);%607-707
w_dq_b_i_b.qd.v(3,:) = DV(7*N+8:8*N+8);%708-808


q_dq_b_i.qr.s = DV(8*N+9:9*N+9);%809-909
q_dq_b_i.qr.v(1,:) = DV(9*N+10:10*N+10);%910-1010
q_dq_b_i.qr.v(2,:) = DV(10*N+11:11*N+11);%1011-1111
q_dq_b_i.qr.v(3,:) = DV(11*N+12:12*N+12) ; %1112-1212

q_dq_b_i.qd.s = DV(12*N+13:13*N+13);%1213-1313 
q_dq_b_i.qd.v(1,:) = DV(13*N+14:14*N+14);%1314-1414
q_dq_b_i.qd.v(2,:) = DV(14*N+15:15*N+15);%1415-1515
q_dq_b_i.qd.v(3,:) = DV(15*N+16:16*N+16);%1516-1616

F_dq.qr.s = DV(16*N+17:17*N+17);% This will be all zero always 
F_dq.qr.v(1,:) = DV(17*N+18:18*N+18);%
F_dq.qr.v(2,:) = DV(18*N+19:19*N+19);%
F_dq.qr.v(3,:) = DV(19*N+20:20*N+20) ; % 

F_dq.qd.s = DV(20*N+21:21*N+21);%
F_dq.qd.v(1,:) = DV(21*N+22:22*N+22);%
F_dq.qd.v(2,:) = DV(22*N+23:23*N+23);%
F_dq.qd.v(3,:) = DV(23*N+24:24*N+24);%

%%
% q_dq_b_i.qr = quaternion(q_dq_b_i.qr.s,(q_dq_b_i.qr.v)');
% q_dq_b_i.qd = quaternion(q_dq_b_i.qd.s,(q_dq_b_i.qd.v)');
% q_dq_b_i = dualquaternion(q_dq_b_i.qr,q_dq_b_i.qd);
% 
% 
% w_dq_b_i_b.qr = quaternion(w_dq_b_i_b.qr.s,(w_dq_b_i_b.qr.v)');
% w_dq_b_i_b.qd = quaternion(w_dq_b_i_b.qd.s,(w_dq_b_i_b.qd.v)');
% w_dq_b_i_b = dualquaternion(w_dq_b_i_b.qr,w_dq_b_i_b.qd);
% 
% A = cross((w_dq_b_i_b),J_dq*(w_dq_b_i_b));

q_dq_b_i
q_dq_b_i= dq_mat(q_dq_b_i)
w_dq_b_i_b = dq_mat(w_dq_b_i_b)
F_dq = dq_mat(F_dq)
pdt = q_dq_b_i.*w_dq_b_i_b
A = cross_times((w_dq_b_i_b),J_dq.*(w_dq_b_i_b))
D = inv(J_dq).*(minus_times(F_dq,A))

% q_dq_b_i_dot_con = (2/(tau_f-tau0)).*(D*q_dq_b_i-0.5*q_dq_b_i.*w_dq_b_i_b);% 100x8 
% w_dq_b_i_b_dot_con = (2/(tau_f-tau0)).*(D*w_dq_b_i_b - inv(J_dq)*(F_dq - A));



ceq = [q_dq_b_i_dot_con.qr.s;(w_dq_b_i_b_dot_con.qr.v)';q_dq_b_i_dot_con.qd.s,(w_dq_b_i_b_dot_con.qd.v)',w_dq_b_i_b.qr.s(1,1)-bc(1);w_dq_b_i_b.qr.v(1,1)-bc(2);w_dq_b_i_b.qr.v(2,1)-bc(3);w_dq_b_i_b.qr.v(3,1)-bc(4);q_dq_b_i.qr.s(1,1)-bc(5);q_dq_b_i.qr.v(1,1)-bc(6);q_dq_b_i.qr.v(2,1)-bc(7);q_dq_b_i.qr.v(3,1)-bc(8);norm(q_dq_b_i_dot_con.qr)-1; norm(w_dq_b_i_b_dot_con)-1;dot(q_dq_b_i_dot_con.qr,q_dq_b_i_dot_con.qd),dot(w_dq_b_i_b_dot_con.qr,w_dq_b_i_b_dot_con.qd)];
c = [w_dq_b_i_b.qr.s(1,end)-bc(9); w_dq_b_i_b.qr.v(1,end)-bc(10);w_dq_b_i_b.qr.v(2,end)-bc(11);w_dq_b_i_b.qr.v(3,end)-bc(12);w_dq_b_i_b.qd.s(1,end)-bc(13);w_dq_b_i_b.qd.v(1,end)-bc(14);w_dq_b_i_b.qd.v(2,end)-bc(15);w_dq_b_i_b.qd.v(3,end)-bc(16)];
end




function dq = dq_mat(dq_struct)

   n = size(dq_struct.qr.s);
   
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