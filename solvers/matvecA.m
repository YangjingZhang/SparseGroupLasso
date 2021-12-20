%% generate the matrix of the Newton system
%% V = I + sigma A M(u) A^*
%% A is an operator

function [Vx,Proj2_Pv,v] = matvecA(u,Amap,ATmap,c,P,sig,x)
c1 = c(1); c2 =  c(2);
v = ProxL1(u,c1);
Pv = P.times(v);
n = length(x);
[Proj2_Pv,grp_nrm] = P.ProjL2(Pv,c2);
supp_v2 = (v~=0);% 
if all(~supp_v2)
    Vx = x;
    return;
end
if isempty(find(grp_nrm ~= 0,1))
    Vx = x;
    return;
end
Atx = ATmap(x);
TAtx = supp_v2.*Atx;
PTAtx = P.matrix*TAtx;
Vx2 = PTAtx;
for k = find(grp_nrm ~= 0)'
    kstart = P.ind(1,k);
    kend = P.ind(2,k);
    G_k = P.G(kstart:kend);
    v_k = v(G_k);
    
    cw = c2*P.ind(3,k);
    par1 = cw/grp_nrm(k);
    par1 = sig*par1;
    par2 = par1/grp_nrm(k)^2;   

    PTAtx_k = PTAtx(kstart:kend);
    Vx2(kstart:kend) = par1*PTAtx_k - (par2*(v_k'*PTAtx_k))*v_k;
end
Vx3 = P.matrix'*Vx2;
Vx3 = sig*TAtx - Vx3;
Vx3 = Amap(Vx3);
Vx = x + Vx3;
end