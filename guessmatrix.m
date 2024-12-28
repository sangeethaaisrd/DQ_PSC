%clc
clear all
close all
i=0;
qi = normalize(quaternion(1,[0,0,0]));
qn  = normalize(quaternion(1,[2,3,4]));
eps = 0.01;
i=1;
for t =0:0.01:1

quat(i)= slerp (qi, qn, t, eps);
quat(i) = normalize(quat(i));
i=i+1;

end

 

