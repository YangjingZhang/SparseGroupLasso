% v=argmin_v  kappa |v|_1+0.5*|v-w|^2

function v = ProxL1(w, kappa)
if kappa < 0; 
   error('kappa needs to be nonegative');
end
%v = sign(w).*max(abs(w) - kappa,0); %
v = max(0, w-kappa)-max( 0, -w-kappa);
