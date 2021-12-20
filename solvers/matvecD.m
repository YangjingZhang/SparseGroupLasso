%% matrix 
%% D = [B C] woodbury formula
function [D,sp_dim,id_yes] = matvecD(u,A,c,P,sig)
c1 = c(1); c2 = c(2);
v = ProxL1(u,c1);

Pv = P.times(v);
[n,p] = size(A);
id_yes = 0;
[~,grp_nrm] = P.ProjL2(Pv,c2);
I_grp = find(grp_nrm ~= 0);
supp_v = (v~=0);
r = nnz(supp_v);
r2 = length(I_grp);
D = []; sp_dim = 0;
if ~(r && r2) 
    id_yes = 1;
    return;
end
B = zeros(n,r + r2);
C = zeros(n,r2);
s_start = 1; i = 1;
for k = I_grp'
    G_k = P.G(P.ind(1,k):P.ind(2,k));
    v_k = v(G_k);  
    indvk = (v_k~=0);    
    cw = c2*P.ind(3,k);
    par1 = cw/grp_nrm(k);
    par1 = sig*par1;
    par2 = par1/grp_nrm(k)^2;     

    Al = A(:,G_k);%
    Al = Al(:,indvk);%
    
    Bl = sqrt(sig-par1)*Al;
    lenind1 = size(Bl,2);
    s_end = s_start + lenind1;
    B(:,s_start:s_end-1) = Bl;
    s_start = s_end;  
        
    if nnz(v_k)
        pv = v_k(indvk);%
        Asl = Al*pv;
        cl = sqrt(par2)*Asl;
        C(:,i) = cl;
        i = i + 1;       
    end

end
sp_dim = s_start + r2- 1 ;
B(:,s_start:sp_dim) = C;
D = B(:,1:sp_dim);
end
