%%*************************************************************************
%% SSNAL:
%% A semismooth Newton augmented Lagrangian method for solving 
%% sparse group lasso problems with non-overlapped groups:
%% (P) min   (1/2)|| A x - b||^2 + p(x),
%% where p(x) = c1 \|x\|_1 + c2 * sum_i w_i ||x_{G_i}||, c1>0, c2>0, w_i>0
%% x \in R^p, b \in R^n
%%
%% (D) - min (1/2)||y||^2 + <b,y> + P^*(z) s.t. A^* y + z = 0.
%%
%%*************************************************************************
%% SSNAL: 
%% Copyright (c) 2019 by
%% Defeng Sun, Kim-Chuan Toh, Ning Zhang, and Yangjing Zhang
%%*************************************************************************
%% Input - prescribed group G_i, i = 1,...,g,
%%       - G \in R^s, a column vector containing the indices of all the 
%%         groups, i.e., G =[G_1; G_2; ...; G_g],
%%       - ind \in R^{3*g}, indices, G(ind(1,i):ind(2,i)) denotes the indices
%%         for the i-th group,
%%         ind(3,i)= w_i denotes the weight for the i-th group, 
%%         e.g., we can take (w_i)^2 = |G_i|.
%%

function[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,p,c,G,ind,options,y0,z0,x0)

rng('default')
stoptol = 1e-6;  
printyes = 1; 
maxit = 5000;  
stopopt = 2; %1:relgap+feas 2:kkt
if isfield(options,'stoptol');     stoptol = options.stoptol; end
if isfield(options,'stopopt');     stopopt = options.stopopt; end
if isfield(options,'printyes');    printyes = options.printyes; end
if isfield(options,'maxit');       maxit = options.maxit; end
if isfield(options,'Lip');         Lip = options.Lip; end
if isfield(options,'stoptol_gap'); stoptol_gap  = options.stoptol_gap ; end

if (stopopt == 3) && ~exist('stoptol_gap','var')
    stoptol_gap = stoptol;
end
%%
fprintf('\n*************************************************************************************');
fprintf('\n SparseGroupLasso');
fprintf('\n Authors: Defeng Sun, Kim-Chuan Toh, Ning Zhang, Yangjing Zhang                      ');
fprintf('\n*************************************************************************************');

%% Amap and ATmap
tstart = clock;
tstart_cpu = cputime;
n = length(b);
normb = 1 + norm(b);
if isstruct(Ainput)
    if isfield(Ainput,'A');
        A = Ainput.A;
        Amap0 = @(x) A*x;
    elseif isfield(Ainput,'Amap')
        Amap0 = Ainput.Amap;
    else
        error(' Ainput.Amap not defined')
    end
    if exist('A','var')
        ATmap0 = @(y) A'*y;
    elseif isfield(Ainput,'ATmap')
        ATmap0 = Ainput.ATmap;
    else
        error(' Ainput.ATmap not defined')
    end
else
    A = Ainput;
    Amap0 = @(x) A*x;
    ATmap0 = @(y) A'*y;
end
existA = exist('A','var');
if ~existA
  fprintf('\n Remark: better performance can be achieved if the'); 
  fprintf('\n matrix representation of the linear operator A is available');
end

AATmap0 = @(x) Amap0(ATmap0(x));
if ~exist('Lip','var')
    eigsopt.issym = 1;
    tstartLip = clock;
    Lip = eigs(AATmap0,m,1,'LA',eigsopt);
    fprintf('\n Lip = %3.2e, time = %3.2f', Lip, etime(clock, tstartLip));
end
if ~exist('y','var') || ~exist('z','var') || ~exist('x','var')
    y = zeros(n,1); z = zeros(p,1); x = z;
else
    y = y0; z = z0; x = x0;
end
%% group map
P = Def_P(p,G,ind);
%%
parmain.tstart = tstart;
parmain.tstart_cpu = tstart_cpu;
parmain.stoptol = stoptol;
parmain.Lip = Lip;
if isfield(options,'sigma'); parmain.sigma = options.sigma; end
parmain.maxit = maxit;
parmain.stopopt = stopopt;
parmain.printyes = printyes;
parmain.existA = existA;
if existA; parmain.A = A; end
parmain.p = p;
parmain.n = n;
[obj_main,y,z,x,info_main,runhist_main] = ...
    SGLasso_SSNAL_main(Amap0,ATmap0,b,c,P,parmain,y,z,x);
iter = info_main.iter;
ttime = etime(clock,tstart);
ttime_cpu = cputime - tstart_cpu;
msg = info_main.msg;
if (iter == maxit)
    msg = ' maximum iteration reached';
end
Aty = ATmap0(y); Ax = Amap0(x);
Rd = Aty + z;
Px = P.matrix*(x);
normRd = norm(Rd);
dualfeas = normRd/(1+norm(z));
Axb = Ax - b;
Rp = Axb - y;
normRp = norm(Rp);
primfeas = normRp/normb;
maxfeas = max(primfeas,dualfeas);

eta = norm(x - Prox_p(x+z,c,P))/(1+norm(x));

lasso = c(2)*P.Lasso_fz(Px) + c(1)*sum(abs(x));
dualobj = -norm(y)^2/2 - b'*y;
primobj = norm(Axb)^2/2 + lasso;
obj = [primobj;dualobj];
gap = primobj - dualobj;
relgap = abs(gap)/(1+abs(primobj)+abs(dualobj));
    


runhist.totaltime = ttime;
runhist.primobj = primobj;
runhist.dualobj = dualobj;
runhist.maxfeas = maxfeas;
runhist.eta = eta;
info.relgap = relgap;
info.iter = iter;
info.time = ttime;
info.time_cpu = ttime_cpu;
info.eta = eta;
info.obj = obj;
info.maxfeas = maxfeas;
if printyes
    fprintf('\n****************************************\n');
    fprintf([' SSNAL       : ',msg,'\n']);
    fprintf(' iteration   : %d\n',iter);
    fprintf(' time        : %3.2f\n',ttime);
    fprintf(' cputime     : %3.2f\n',ttime_cpu);
    fprintf(' prim_obj    : %4.8e\n',primobj);
    fprintf(' dual_obj    : %4.8e\n',dualobj);
    fprintf(' relgap      : %4.5e\n',relgap);
    fprintf(' primfeas    : %3.2e\n',primfeas);
    fprintf(' dualfeas    : %3.2e\n',dualfeas);
    fprintf(' eta         : %3.2e\n',eta);
    fprintf(' nnz         : %d\n',cardcal(x,0.999));
end

