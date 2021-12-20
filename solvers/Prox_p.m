%% proximal mapping of the non-overlapped group lasso penalty function 
%% P(x) = c1 \|x\|_1 + c2 * sum_i w_i ||x_{G_i}|| . 
%% June 26,2017 Zhang Yangjing

function [u] = Prox_p(v,c,P)
c1 = c(1); c2 = c(2);

utmp = ProxL1(v,c1);
u = P.times(utmp); 

[u] = P.ProjL2(u,c2);
u = utmp - P.trans(u);

