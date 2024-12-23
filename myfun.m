%% 6DOF rigid body dynamics equation in DQ form 
function f = myfun(x,a1)
q1 = quaternion(x(1),[x(2),x(3),x(4)]);
q2 = quaternion(x(5),[x(6),x(7),x(8)]);
func = a1*(mtimes(q1, q1)+ mtimes(q2, q2));
f = sqrt(func.s^2+func.v(1)^2+func.v(2)^2+func.v(3)^2);
end 