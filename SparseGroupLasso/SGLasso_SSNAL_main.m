function[obj,y,z,x,info,runhist] = SGLasso_SSNAL_main(Amap0,ATmap0,b,c,P,parmain,y,z,x)
tstart = parmain.tstart;
tstart_cpu = parmain.tstart_cpu;
Lip = parmain.Lip;
maxit = parmain.maxit;
printyes = parmain.printyes;
stoptol = parmain.stoptol;
stopopt = parmain.stopopt;
p = parmain.p; n = parmain.n;
existA = parmain.existA;
if existA; A = parmain.A; end  
c1 = c(1); c2 = c(2);
stop = 0;

sigmaLip = 1/Lip;

Aty = ATmap0(y);  Ax = Amap0(x); Px = P.matrix*(x);

obj(1) = 0.5*norm(Ax - b)^2 + c2*P.Lasso_fz(Px) + c1*sum(abs(x));
obj(2) = 0;
normb = 1 + norm(b);
Ainput_nal.Amap = Amap0;
Ainput_nal.ATmap = ATmap0;
if existA; Ainput_nal.A = A; end
sig = max(1/sqrt(Lip),min([1,sigmaLip,c1,1/c1,1/c2]));
Rp = Ax - b - y;
Rd = Aty + z;
primfeas = norm(Rp)/normb;
dualfeas = norm(Rd)/(1+norm(z));
maxfeas = max(primfeas,dualfeas);
relgap = abs(obj(1) - obj(2))/(1+abs(obj(1))+abs(obj(2)));

if printyes
%     fprintf('\n ***********************************************');
%     fprintf('******************************************');
%     fprintf('\n \t\t GROUP LASSO: SSNAL      ');
%     fprintf('\n *********************************************');
%     fprintf('*******************************************\n');
    fprintf('\n n=%d, m=%d, ',p,n);
    fprintf('tol=%1.1e, parameters:c1=%4.3f, c2=%4.3f\n',stoptol,c1,c2);
    fprintf('\n ---------------------------------------------------');
    fprintf('\n iter|  [pinfeas  dinfeas]    relgap|    pobj          dobj       |');
    fprintf(' time | sigma |');
    fprintf('\n*********************************************');
    fprintf('*******************************************************');
    fprintf('\n #%3.1d|  %3.2e %3.2e  %- 3.2e %- 8.7e %- 8.7e  %5.1f',...
       0,primfeas,dualfeas,relgap,obj(1),obj(2),etime(clock,tstart)); 
    fprintf('  %3.2e ',sig);
end

%% ssncg
parNCG.tolconst = 0.5;
parNCG.p = p;

maxitersub = 10;
prim_win = 0;
dual_win = 0;

ssncgop.existA = existA;
ssncgop.tol = stoptol;
ssncgop.precond = 0;

sigmamax = 1e6;%1e6;%1e6;%1
sigmamin = 1e-4;%1e-5;%1e-5; %1


%% main loop
for iter = 1:maxit
    parNCG.sigma = sig;
    
    if dualfeas < 1e-5
        maxitersub = max(maxitersub,30);
    elseif dualfeas < 1e-3
        maxitersub = max(maxitersub,30);
    elseif dualfeas < 1e-1
        maxitersub = max(maxitersub,20);
    end
    ssncgop.maxitersub = maxitersub;
    %% Newton direction
    [y,z,Aty,Prox_u,x,Ax,parNCG,runhist_NCG,info_NCG] = ...
        SGLasso_SSNCG(y,Aty,x,Ax,Ainput_nal,b,c,P,parNCG,ssncgop);
    if info_NCG.breakyes < 0
        parNCG.tolconst = max(parNCG.tolconst/1.06,1e-3);
    end
    
    %% compute things: pobj,dobj,infeas,gap,...
    Rd = Aty + z;
    Px = P.matrix*(x);
    normRd = norm(Rd);
    dualfeas = normRd/(1+norm(z));
    Axb = Ax - b;
    Rp = Axb - y;
    normRp = norm(Rp);
    primfeas = normRp/normb;

    lasso = c2*P.Lasso_fz(Px) + c1*sum(abs(x));
    dualobj = -norm(y)^2/2 - b'*y;
    primobj = norm(Axb)^2/2 + lasso;

    gap = primobj - dualobj;
    relgap = abs(gap)/(1+abs(primobj)+abs(dualobj));
    
    if stopopt == 1
        stop = max(dualfeas,relgap) < stoptol;
    elseif stopopt == 2
        eta_tmp = max(dualfeas,primfeas);
        if eta_tmp < stoptol
            eta_1 = norm(x - Prox_p(x+z,c,P))/(1+norm(x));
            stop = eta_1 < stoptol;
        end
    elseif stopopt == 3
        stop = (dualfeas < stoptol) && (relgap < stoptol_gap);
    elseif stopopt ==4
        stop = (dualfeas < stoptol) && (gap < stoptol*norm(b,2)^2);

    end
    ttime = etime(clock,tstart);
    runhist.primfeas(iter)    = primfeas; 
    runhist.dualfeas(iter)    = dualfeas; 
    runhist.sigma(iter)       = sig; 
    runhist.primobj(iter)     = primobj;
    runhist.dualobj(iter)     = dualobj;
    runhist.gap(iter)         = gap;
    runhist.relgap(iter)      = relgap; 
    runhist.ttime(iter)       = ttime;
    runhist.itersub(iter)     = info_NCG.itersub;
    
    if (printyes)
        fprintf('\n %5.0d| [%3.2e %3.2e]  %- 3.2e| %- 5.4e %- 5.4e |',...
               iter,primfeas,dualfeas,relgap,primobj,dualobj); 
        fprintf(' %5.1f| %3.2e|',ttime, sig);   
    end
    %% check termination
    if stop || (iter == maxit) 
        termination = 'converged';
        if iter == maxit, termination = 'maxiter reached'; end
        runhist.termination = termination;        
        runhist.iter = iter;
        runhist.nnz = cardcal(x,0.999);
        obj(1) = primobj;
        obj(2) = dualobj;              
        break;
    end
    %% update sigma
    ratio = primfeas / (dualfeas + eps);
    runhist.ratio_seq(iter) = ratio;

    if ratio < 1
        prim_win = prim_win + 1;
    else
        dual_win = dual_win + 1;
    end
    sigma_update_iter = sigma_fun(iter); 
    if primfeas > 100*stoptol;
        if runhist_NCG.av_findstep > 5
            sigmascale = sqrt(3);
        else
            sigmascale = 3;
        end
    else
        if runhist_NCG.av_findstep > 5
            sigmascale = sqrt(5); %sqrt(5);
        else
            sigmascale = 5;
        end
    end
    
    if (rem(iter,sigma_update_iter)==0) &&  info_NCG.breakyes < 0 
        if (prim_win > max(1,1.2*dual_win)) 
            prim_win = 0;
            sig = min(sigmamax,sig*sigmascale);
        elseif(dual_win > max(1,1.2*prim_win))
            dual_win = 0;
            sig = max(sigmamin,sig/sigmascale);
        end
        if info_NCG.breakyes >=0 && iter >= 10
            sig = max(sigmamin,2*sig/sigmascale);
        end
    end    
end
if stopopt == 2
    kktres = max(eta_tmp,eta_1); 
    runhist.kktres = kktres;
    info.kktres = kktres;
end 
info.maxfeas = maxfeas;
info.iter = iter;
info.ttime = ttime;
info.termination = termination;
info.relgap = relgap;
info.msg = termination;
end
function sigma_update_iter = sigma_fun(iter)   
  if (iter < 10)
     sigma_update_iter = 2; 
  elseif (iter < 20) 
     sigma_update_iter = 3; 
  elseif (iter < 200)
     sigma_update_iter = 3; %5
  elseif (iter < 500)
     sigma_update_iter = 10; 
  end
end
  
