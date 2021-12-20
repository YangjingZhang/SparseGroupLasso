%% Define P map. 
%% Zhang Yangjing, 05 May 2017
%
%
% Input - prescribed group G_i, i=1,...,g
%       - G\in R^s, a column vector containing the indices of all the 
%         overlapping groups. G =[ G_1; G_2; ...; G_g]
%       - ind \in R^{3*g}, indices, G(ind(1,i):ind(2,i))denotes the indices
%         for the i-th group,
%         ind(3,i)= wi denotes the weight for the i-th group, 
%         e.g., we can take (wi)^2 = |G_i|.
%
% Output - P
%
% P.matrix   - matrix representation of P  
%
% P.grpNUM   - number of groups
% P.grpLEN   - a column vector [|G_1|, |G_2|, ... ,|G_grpNUM|]'
% P.ntotal   - |G_1|+|G_2|+...+|G_grpNUM|
%
%
% P.times    - P.times(x) computes P(x)%active
% P.trans    - P.trans(y) computes P'(y)%active
%
% P.Pi       - P.Pi(i) matrix representation of Pi, Pi(x)=x_{G_i}
%
%
% P.ProjL2   - P.ProjL2(z,c2) computes projection onto B_2 = B1*B2*...
%              Bi={zi \in R^|G_i|: norm(zi,2)<=c2*wi},i=1,2,...,grpNUM %active
%
% P.Lasso_fx - P.Lasso_fx(x) computes sum_i {wi*norm(P_i(x),2)}
% P.Lasso_fz - P.Lasso_fz(z) computes sum_i {wi*norm(zi,2)}

function P = Def_P(n,G,ind)

grpNUM = size(ind,2);
P.grpNUM = grpNUM;
grpLEN = (ind(2,:) - ind(1,:))' + 1;
P.grpLEN = grpLEN;

ntotal = sum(grpLEN);
P.ntotal = ntotal;
P.ind = ind;
P.G = G;

I = [1:ntotal];
J = G';
V = ones(1,ntotal);
Pma = sparse(I,J,V);
P.matrix = Pma;
 
P.Pi = @(i) P_i(i,G,ind,n,grpLEN);

P.times = @(x) Pma*x;
P.trans = @(y) Pma'*y;

P.ProjL2 = @(z,c1) mexProjL2(z,c1,ind,grpNUM);
P.ProxL2 = @(z,c1) mexProxL2(z,c1,ind,grpNUM);

P.Lasso_fx = @(x)mexfz(Pma*x,ind,grpNUM);
P.Lasso_fz = @(z)mexfz(z,ind,grpNUM);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pi] = P_i(i,G,ind,n,grpLEN)
tmp = grpLEN(i);
I = [1:tmp];
J = G(ind(1,i):ind(2,i));
V = ones(1,tmp);
Pi = sparse(I,J,V,tmp,n);
end


