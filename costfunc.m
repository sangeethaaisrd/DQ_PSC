function cost = costfunc(DV, w, tau0, tau_f)
%COSTFUNC Cost function to minimize control effort 
%   DV - Decision Vector
N = (length(DV)-22)/22;
Fx = DV(16*N+17:17*N+17);
Fy = DV(17*N+18:18*N+18);
Fz = DV(18*N+19:19*N+19);
Fx = DV(19*N+20:20*N+20);
Fy = DV(20*N+21:21*N+21);
Fz = DV(21*N+22:22*N+22);
Tnormsq = Tx.^2+Ty.^2+Tz.^2;
Fnormsq = Fx.^2+Fy.^2+Fz.^2;
tot_norm = Tnormsq + Fnormsq;
cost = (tot_norm'*w)*(tau_f-tau0)/2;  %integration weights using clenshaw curtis - w
end

