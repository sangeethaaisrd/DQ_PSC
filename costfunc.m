function cost = costfunc(DV, w, tau0, tau_f)
%COSTFUNC Cost function to minimize control effort 
%   DV - Decision Vector

N = (length(DV)-24)/24;
Fx = DV(17*N+18:18*N+18);
Fy = DV(18*N+19:19*N+19);
Fz = DV(19*N+20:20*N+20) ;
Tx = DV(21*N+22:22*N+22);
Ty = DV(22*N+23:23*N+23);
Tz =  DV(23*N+24:24*N+24);
Tnormsq = Tx.^2+Ty.^2+Tz.^2;
Fnormsq = Fx.^2+Fy.^2+Fz.^2;
tot_norm = Tnormsq + Fnormsq;
cost = (tot_norm'*w)*(tau_f-tau0)/2;  %integration weights using clenshaw curtis - w
end

