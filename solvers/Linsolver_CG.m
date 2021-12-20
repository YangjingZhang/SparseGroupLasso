function [dy,resnrm,solve_ok] = Linsolver_CG(Ainput,rhs,u,c,P,par)
n = length(rhs);
Ayes = isfield(Ainput,'A');
solver = 3;%1:p direct %2:d direct woodbury-formula %3:d pcg 
density = par.nnz; %l_1 sparsity
dn = 10000;
sig = par.sigma;

if n <= 1000 && Ayes
    solver = 1;
end

if density <= n && Ayes && density <= dn
    solver = 2;
end

if (n > 5e3 && density >= 1000) || (n > 2000 && density > 5000) || (n > 100 && density > 8000)
   solver = 3;
end


if solver == 1 
    V = mat_ssn(u,Ainput.A,c,P,sig);
    if n <= 1000
        dy = V\rhs;
    else
        LAAt = mycholAAt(V,n);
        dy = mylinsysolve(LAAt,rhs);
    end
    resnrm = 0; solve_ok = 1;
end
if solver == 2 
    [V2,D,sp_dim,id_yes] = mat2_ssn(u,Ainput.A,c,P,sig);
    if id_yes
        dy = rhs;
    elseif sp_dim <= 1000
        rhstmp = D'*rhs;
        dy = V2\rhstmp;
        dy = D*dy;
        dy = rhs - dy;        
    else
        rhstmp = D'*rhs;
        LAAt2 = mycholAAt(V2,sp_dim);
        dy = mylinsysolve(LAAt2,rhstmp);
        dy = D*dy;
        dy = rhs - dy;
    end

    resnrm = 0; solve_ok = 1;
end

if solver == 3 
    if Ayes
        [D,~,id_yes] = matvecD(u,Ainput.A,c,P,sig);
        if id_yes
            dy = rhs;
        else
            rhstmp = D'*rhs;
            Afun = @(x) (x+D'*(D*x));
            [dy,~,resnrm,solve_ok] = psqmrGL(Afun,rhstmp,par);
            dy = D*dy;
            dy = rhs - dy;
        end
    else
        Afun = @(x) matvecA(u,Ainput.Amap,Ainput.ATmap,c,P,sig,x);
        [dy,~,resnrm,solve_ok] = psqmrGL(Afun,rhs,par);
    end
end
end
    