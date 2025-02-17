 %% 6DOF rigid body dynamics equation in DQ form - constraint 
function[c,ceq] = nonlcon(DV,J_dq,tau0,tau_f,D,bc)

DV =DV';

N = (length(DV)-24)/24 ;% in this case N =100 

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



q_dq_b_i= dq_mat(q_dq_b_i);

w_dq_b_i_b = dq_mat(w_dq_b_i_b);
F_dq = dq_mat(F_dq);
pdt = q_dq_b_i.*w_dq_b_i_b;

pdt1 = 0.5*(dqarray_vec(pdt))';

A = cross_times((w_dq_b_i_b),J_dq.*(w_dq_b_i_b));
Dy = inv(J_dq).*(minus_times(F_dq,A));
Dy = dqarray_vec(Dy);

q_dq_b_i =dqarray_vec(q_dq_b_i);
w_dq_b_i_b =dqarray_vec(w_dq_b_i_b);

A1 = (2/(tau_f-tau0)).*(D*q_dq_b_i'- pdt1);

B1 = (2/(tau_f-tau0)).*(D*w_dq_b_i_b' - Dy');


q_dq_b_i_dot_con.qr.s = A1(:,1);
q_dq_b_i_dot_con.qr.v = A1(:,2:4);
q_dq_b_i_dot_con.qd.s = A1(:,5);
q_dq_b_i_dot_con.qd.v = A1(:,6:8);

w_dq_b_i_b_dot_con.qr.s = B1(:,1);
w_dq_b_i_b_dot_con.qr.v = B1(:,2:4);
w_dq_b_i_b_dot_con.qd.s = B1(:,5);
w_dq_b_i_b_dot_con.qd.v = B1(:,6:8);

[norm_q_dq_b_i_dot_con,dot_q_dq_b_i_dot_con] = norm_dot(q_dq_b_i_dot_con,N);
[norm_w_dq_b_i_b_dot_con,dot_w_dq_b_i_b_dot_con] = norm_dot(w_dq_b_i_b_dot_con,N);
norm_q_dq_b_i_dot_con = norm_q_dq_b_i_dot_con';
norm_w_dq_b_i_b_dot_con = norm_w_dq_b_i_b_dot_con';
dot_q_dq_b_i_dot_con =dot_q_dq_b_i_dot_con';
dot_w_dq_b_i_b_dot_con = dot_w_dq_b_i_b_dot_con';
 

 ceq = [q_dq_b_i_dot_con.qr.s;q_dq_b_i_dot_con.qr.v(:,1);q_dq_b_i_dot_con.qr.v(:,2);q_dq_b_i_dot_con.qr.v(:,3);q_dq_b_i_dot_con.qd.s;q_dq_b_i_dot_con.qd.v(:,1);q_dq_b_i_dot_con.qd.v(:,2);q_dq_b_i_dot_con.qd.v(:,3);w_dq_b_i_b_dot_con.qr.s;w_dq_b_i_b_dot_con.qr.v(:,1);w_dq_b_i_b_dot_con.qr.v(:,2);w_dq_b_i_b_dot_con.qr.v(:,3);w_dq_b_i_b_dot_con.qd.s;w_dq_b_i_b_dot_con.qd.v(:,1);w_dq_b_i_b_dot_con.qd.v(:,2);w_dq_b_i_b_dot_con.qd.v(:,3);w_dq_b_i_b(1,1)-bc(1);w_dq_b_i_b(2,1)-bc(2);w_dq_b_i_b(3,1)-bc(3);w_dq_b_i_b(4,1)-bc(4);q_dq_b_i(1,1)-bc(5);q_dq_b_i(2,1)-bc(6);q_dq_b_i(3,1)-bc(7);q_dq_b_i(4,1)-bc(8);norm_q_dq_b_i_dot_con-1; norm_w_dq_b_i_b_dot_con-1;dot_q_dq_b_i_dot_con;dot_w_dq_b_i_b_dot_con];
 c = [w_dq_b_i_b(1,end)-bc(17,1); w_dq_b_i_b(2,end)-bc(18);w_dq_b_i_b(3,end)-bc(19);w_dq_b_i_b(4,end)-bc(20);w_dq_b_i_b(5,end)-bc(21);w_dq_b_i_b(6,end)-bc(22);w_dq_b_i_b(7,end)-bc(23);w_dq_b_i_b(8,end)-bc(24)];
end

function [abs_q,dot_dq] = norm_dot(dq_mat,N)

for i = 1:1:N+1
    abs_q(i)= norm(quaternion(dq_mat.qr.s(i,1),[dq_mat.qr.v(i,1:3)]));
    qr(i) = quaternion(dq_mat.qr.s(i,1),[dq_mat.qr.v(i,1:3)]);
    qd(i) = quaternion(dq_mat.qd.s(i,1),[dq_mat.qd.v(i,1:3)]);
    dot_dq(i) = dot(qr(i),qd(i));
end
end 


function dq= dq_mat(dq_struct)

   n = size(dq_struct.qr.s);
   
   for i=1:1:n(2)
            qr(i) = quaternion(dq_struct.qr.s(:,i),[dq_struct.qr.v(:,i)]);
            qd(i) = quaternion(dq_struct.qd.s(:,i),[dq_struct.qd.v(:,i)]); 
            dq(1,i) = dualquaternion(qr(i),qd(i));
   end
end

function dqvec_array = dqarray_vec(arr)

   n = size(arr);
   
   for i=1:1:n(2)
           dqvec_array(:,i) =  dq2vec(arr(1,i));
   end
end
