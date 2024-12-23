%% 6DOF rigid body dynamics equation in DQ form 
function [c,ceq] = mycon(x)
q1 = quaternion(x(1),[x(2),x(3),x(4)]);
q2 = quaternion(x(5),[x(6),x(7),x(8)]);

ceq = [norm(q1)-1 norm(q2)-1];
c =[];
end 