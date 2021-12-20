%% generate the matrix of the Newton system
%% V = I + sigma A M(u) A^*

function [V,Proj2_Pv,v] = mat_ssn(u,A,c,P,sig)
c1 = c(1); c2 =  c(2);
v = ProxL1(u,c1);
Pv = P.times(v);
[n,p] = size(A);
[Proj2_Pv,grp_nrm] = P.ProjL2(Pv,c2);
supp_v2 = (v~=0);% 
if all(~supp_v2)
    V = speye(n);
    return;
end
if isempty(find(grp_nrm ~= 0,1))
    V = speye(n);
    return;
end
V = zeros(n,n);

for k = find(grp_nrm ~= 0)'
    G_k = P.G(P.ind(1,k):P.ind(2,k));
    v_k = v(G_k);
    indvk = (v_k~=0);
    
    cw = c2*P.ind(3,k);
    par1 = cw/grp_nrm(k);
    par1 = sig*par1;
    par2 = par1/grp_nrm(k)^2;   

    Al = A(:,G_k);%
    Al = Al(:,indvk);%
    
    M1 = Al*Al';      
    V = V + (sig-par1)*M1;
    
        
    if nnz(v_k)
        pv = v_k(indvk);%
        Asl = Al*pv;
        M2 = Asl*Asl';
        V = V + par2*M2;        
    end


end

for i = 1:n; V(i,i) = V(i,i) + 1;end
end