function[lb,ub] = bounds(w_dq_min,w_dq_max,q_dq_min,q_dq_max,F_dq_max,F_dq_min,N)
upto = 3*8;
lb_ = [dq2vec(w_dq_min) ; dq2vec(q_dq_min);dq2vec(F_dq_min)];
ub_ = [dq2vec(w_dq_max) ; dq2vec(q_dq_max);dq2vec(F_dq_max)];
size(lb_);
size(ub_)

lb = 0;
ub = 0;
for i =1:1:upto
    lb =  cat (1,lb,lb_(i,1).*ones(N+1,1))
    ub = cat(1,ub,ub_(i,1).*ones(N+1,1))
end
lb = lb(2:end)
ub = ub(2:end)
size(lb)
size(ub)
end
