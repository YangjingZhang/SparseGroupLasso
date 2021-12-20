function [y,z,Aty,Prox_u,x,Ax,par,runhist,info] = SGLasso_SSNCG(y0,Aty0,x0,Ax0,Ainput,b,c,P,par,options)

printsub = 1;
breakyes = 0;
maxitersub = 50;
tiny = 1e-10;
tol = 1e-6;
maxitpsqmr = 500;
precond = 0;
% Ascaleyes = 0;
if isfield(options,'printsub'); printsub = options.printsub; end
if isfield(options,'maxitersub'); maxitersub = options.maxitersub; end
if isfield(options,'tol'); tol = options.tol; end

existA = options.existA;
sig = par.sigma;
normb = 1+norm(b);
%%
Amap = @(x) Ainput.Amap(x);
ATmap = @(x) Ainput.ATmap(x);
y = y0; Aty = Aty0;
% Ax = Ax0;
% Rp = Ax - b - y;
% normRp = norm(Rp);
u = x0/sig - Aty;
Prox_u = Prox_p(u,c,P);
sigProx_u = sig*Prox_u;
z = u - Prox_u;
psi_y = -(b'*y + 0.5*norm(y)^2 + 0.5*sig*norm(Prox_u)^2);

%% main Newton iteration
for itersub = 1:maxitersub
    x = sigProx_u;
    Ax = Amap(x);
    Grad = - y - b + Ax;
    normGrad = norm(Grad)/normb;    
    priminf_sub = normGrad; 
    Rd = Aty + z;
    normRd = norm(Rd);
    dualinf_sub = normRd/(1+norm(z));
    if max(priminf_sub,dualinf_sub) < tol
       tolsubconst = 0.9;
    else
       tolsubconst = 1e-2;%0.05
    end
    tolsub = max(min(1,par.tolconst*dualinf_sub),tolsubconst*tol);
    runhist.priminf(itersub) = priminf_sub;
    runhist.dualinf(itersub) = dualinf_sub;
    runhist.psi_y(itersub)   = psi_y;

    if printsub
        fprintf('\n      %2.0d  %- 11.10e  %3.2e   %3.2e %3.2e',...
                 itersub,psi_y,priminf_sub,dualinf_sub,par.tolconst);
    end

    if normGrad < tolsub && itersub > 1
        msg = 'good termination in subproblem:';
        if printsub
            fprintf('\n       %s  ',msg);
            fprintf(' dualinfes = %3.2e, normGrad = %3.2e, tolsub = %3.2e',...
                          dualinf_sub,normGrad,tolsub);
        end
        breakyes = -1;
        break;
    end

    
    %% compute Newton direction 
    %% precond = 0, 
    par.epsilon = min(1e-3,0.1*normGrad);
    par.precond = precond;
    if precond == 1
        par.invdiagM = 1/(1+sig);
    end
    if (dualinf_sub > 1e-3) || (itersub <= 5)
        maxitpsqmr = max(maxitpsqmr,200); 
    elseif (dualinf_sub > 1e-4)	 
        maxitpsqmr = max(maxitpsqmr,300); 
    elseif (dualinf_sub > 1e-5)	 
        maxitpsqmr = max(maxitpsqmr,400); 
    elseif (dualinf_sub > 5e-6)
        maxitpsqmr = max(maxitpsqmr,500); 
    end    
    if (itersub > 1) 
        prim_ratio = priminf_sub/runhist.priminf(itersub-1); 
        dual_ratio = dualinf_sub/runhist.dualinf(itersub-1); 
    else
        prim_ratio = 0; dual_ratio = 0;
    end
    rhs = Grad;
    tolpsqmr = min(5e-3,0.001*norm(rhs));%1e-8;%% smaller tolpsqmr often give better primfeasorg.
    const2 = 1;
    if itersub > 1 && (prim_ratio > 0.5 || priminf_sub > 0.1*runhist.priminf(1))
        const2 = 0.5*const2;
    end
	if (dual_ratio > 1.1); const2 = 0.5*const2; end;    
    tolpsqmr = const2*tolpsqmr;
    par.tol = tolpsqmr; par.maxit = maxitpsqmr;
    nnz = length(find(abs(x)>tol));
    par.nnz = nnz;
    [dy,resnrm,solve_ok] = Linsolver_CG(Ainput,rhs,u,c,P,par);
    
    Atdy = ATmap(dy);
    iterpsqmr = length(resnrm)-1;
    if (printsub)
          fprintf('| %3.1e %3.1e %3.0d %-3d',par.tol,resnrm(end),iterpsqmr);
          fprintf(' %2.1f %2.0d',const2,nnz);
    end    
    par.iter = itersub;
    %% line search
    if (itersub <= 3) && (dualinf_sub > 1e-4) || (itersub <3)
        stepop = 1;
    else
        stepop = 2;
    end
    steptol = 1e-5; step_opt.stepop = stepop;

    [psi_y,u,Prox_u,sigProx_u,z,y,Aty,alp,iterstep] = ...
        findstep(b,sig,psi_y,u,Prox_u,sigProx_u,z,y,Aty,dy,Atdy,c,P,steptol,step_opt);

    runhist.solve_ok(itersub) = solve_ok;
    runhist.psqmr(itersub)    = iterpsqmr;
    runhist.findstep(itersub) = iterstep; 
    runhist.av_findstep       = mean(runhist.findstep);
    if alp < tiny; breakyes = 11; end
    psiy_ratio = 1; 
    if itersub > 1
         psiy_ratio = (psi_y - runhist.psi_y(itersub-1))/(abs(psi_y)+eps);
    end
    if printsub
         fprintf(' %3.2e %2.0f',alp,iterstep);
         if (psiy_ratio < 0); fprintf('-'); end
    end
     %% check for stagnation
     if (itersub > 4)
         idx = max(1,itersub-3):itersub; 
         tmp = runhist.priminf(idx); 
         ratio = min(tmp)/max(tmp);
         if (all(runhist.solve_ok(idx) <= -1)) && (ratio > 0.9) ... 
              && (min(runhist.psqmr(idx)) == max(runhist.psqmr(idx))) ...
              && (max(tmp) < 5*tol)
              fprintf('#')
              breakyes = 1; 
         end
         const3 = 0.7;
         priminf_1half  = min(runhist.priminf(1:ceil(itersub*const3))); 
         priminf_2half  = min(runhist.priminf(ceil(itersub*const3)+1:itersub));
         priminf_best   = min(runhist.priminf(1:itersub-1)); 
         priminf_ratio  = runhist.priminf(itersub)/runhist.priminf(itersub-1);
         dualinf_ratio  = runhist.dualinf(itersub)/runhist.dualinf(itersub-1);  
         stagnate_idx   = find(runhist.solve_ok(1:itersub) <= -1);
         stagnate_count = length(stagnate_idx); 
         idx2 = [max(1,itersub-7):itersub]; 
         if (itersub >= 10) && all(runhist.solve_ok(idx2) == -1) ... 
              && (priminf_best < 1e-2) && (dualinf_sub < 1e-3)                    
              tmp = runhist.priminf(idx2); 
              ratio = min(tmp)/max(tmp); 
              if (ratio > 0.5) 
                  if (printsub); fprintf('##'); end
                  breakyes = 2; 
              end
         end
          if (itersub >= 15) && (priminf_1half < min(2e-3,priminf_2half)) ...
               && (dualinf_sub < 0.8*runhist.dualinf(1)) && (dualinf_sub < 1e-3) ...
               && (stagnate_count >= 3) 
               if (printsub); fprintf('###'); end
               breakyes = 3;
          end
          if (itersub >= 15) && (priminf_ratio < 0.1) ...
               && (priminf_sub < 0.8*priminf_1half) ...
               && (dualinf_sub < min(1e-3,2*priminf_sub)) ...
               && ((priminf_sub < 2e-3) || (dualinf_sub < 1e-5 && priminf_sub < 5e-3)) ...
               && (stagnate_count >= 3) 
               if (printsub); fprintf(' $$'); end
               breakyes = 4;  
          end
          if (itersub >=10) && (dualinf_sub > 5*min(runhist.dualinf)) ...
               && (priminf_sub > 2*min(runhist.priminf)) %% add: 08-Apr-2008
               if (printsub); fprintf('$$$'); end
               breakyes = 5;              
          end
          if (itersub >= 20)
             %% add: 12-May-2010
             dualinf_ratioall = runhist.dualinf(2:itersub)./runhist.dualinf(1:itersub-1);
             idx = find(dualinf_ratioall > 1); 
             if (length(idx) >= 3)
                  dualinf_increment = mean(dualinf_ratioall(idx)); 
                  if (dualinf_increment > 1.25)
                      if (printsub); fprintf('^^'); end
                      breakyes = 6;     
                  end                    
             end              
          end 
          if breakyes > 0
              x = sig*Prox_u;
              Ax = Amap(x);
              break
          end
          
     end       
end
info.maxCG = max(runhist.psqmr);
info.avgCG = sum(runhist.psqmr)/itersub;
info.breakyes = breakyes;
info.itersub = itersub;
info.tolconst = par.tolconst;


function [psi_y,u,Prox_u,sigProx_u,z,y,Aty,alp,iter] = ...
        findstep(b,sig,psi_y0,u0,Prox_u0,sigProx_u0,z0,y0,Aty0,dy,Atdy,c,P,tol,options)
    op = options.stepop;
    printlevel = 1;
    maxit = ceil(log(1/(tol+eps))/log(2));
    c1 = 1e-4; c2 = 0.9;
%%    
    tmp1 = dy'*(y0+b);
    tmp2 = norm(dy)^2;
    
    g0 = sig*Atdy'*Prox_u0 - tmp1;
    if g0 <= 0
        alp = 0; iter = 0;
        if printlevel
            fprintf('\n Need an ascent direction, %2.1e  ',g0);
        end
        u = u0;
        y = y0;
        z = z0;
        Aty = Aty0;
        Prox_u = Prox_u0;
        sigProx_u = sigProx_u0;
        psi_y = psi_y0;
        return;
    end
    
   alp = 1; alpconst = 0.5; 
   for iter = 1:maxit
      if iter == 1          
         alp = 1; LB = 0; UB = 1; 
      else
         alp = alpconst*(LB+UB);
      end
      y = y0 + alp*dy;
      u = u0 - alp*Atdy;
      Prox_u = Prox_p(u,c,P);
      z = u - Prox_u;
      galp = tmp1 + alp*tmp2;
      sigProx_u = sig*Prox_u;
      galp = Atdy'*sigProx_u - galp;
      psi_y = -(b'*y + 0.5*norm(y)^2 + 0.5*sig*norm(Prox_u)^2);

      if printlevel > 1
          fprintf('\n --------------------------------------- \n');
          fprintf('\n alp = %2.2f, psi_y = %11.10e, psi_y0 = %11.10e',alp,psi_y,psi_y0);
          fprintf('\n --------------------------------------- \n');
      end      
      if iter == 1
          gLB = g0; gUB = galp;
          if (sign(gLB)*sign(gUB) > 0)%gUB > 0
              if printlevel; fprintf('|');  end
              Aty = Aty0 + alp*Atdy;
              return;
          end
      end
      if (abs(galp) < c2*abs(g0)) && (psi_y - psi_y0 - c1*alp*g0 > 1e-8/max(1,abs(psi_y0)))
         if (op == 1) || ((op == 2) && (abs(galp) < tol))
            if (printlevel); fprintf(':'); end
            Aty = Aty0 + alp*Atdy;
            return;
         end
      end
      if (sign(galp)*sign(gUB) < 0)
         LB = alp; gLB = galp;
      elseif (sign(galp)*sign(gLB) < 0) 
         UB = alp; gUB = galp; 
      end
   end %end maxit
   if iter == maxit
      Aty = Aty0 + alp*Atdy;
   end
   if (printlevel); fprintf('m'); end

