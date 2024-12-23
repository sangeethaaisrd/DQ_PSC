%% Trial quaternion program with fmincon

a1 = 2; % define parameters first
options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
%X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
q1 = quaternion(1,[0,0,0]);
q2 = quaternion(1,[0,0,0]);

x = [q1.s q1.v' q2.s q2.v'];
x0 = [-0.8733 0 0 -0.4872 -0.8733 -0.4872 0 0];
x = fmincon(@(x) myfun(x,a1),x0,[],[],[],[],[],[],@(x) mycon(x),options);