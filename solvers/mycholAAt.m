function [L] = mycholAAt(ATT,m)
   
   pertdiag = 1e-15*ones(m,1); 
   ATT = ATT + spdiags(pertdiag,0,m,m);
%%
   if (nnz(ATT) < 0.2*m*m); 
       use_spchol=1; 
   else
       use_spchol=0; 
   end
   if (use_spchol)
      [L.R,L.p,L.perm] = chol(sparse(ATT),'vector'); 
      L.Rt = L.R'; 
      L.matfct_options = 'spcholmatlab'; 
   else
      if issparse(ATT); ATT = full(ATT); end;           
      L.matfct_options = 'chol';
      L.perm = [1:m]; 
      [L.R,indef] = chol(ATT); 
   end
%%***********************************************************