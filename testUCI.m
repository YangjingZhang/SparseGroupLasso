%% Test on UCI data set
clear; rng('default');
addpath(genpath(pwd));
probname = ['UCIdata',filesep,'pyrim_scale_expanded5.mat'];
load(probname);
%% dimension, options
[n,p] = size(A);
if exist('mexMatvec')
    Amap = @(x) mexMatvec(A,x,0);
    ATmap = @(y) mexMatvec(A,y,1);
else
    AT = A';
    Amap  = @(x) A*x;
    ATmap = @(x) AT*x;
end
AATmap = @(x) Amap(ATmap(x));
eigsopt.issym = 1;
Lip = eigs(AATmap,length(b),1,'LA',eigsopt);
fprintf('\n Lip const = %3.2e, nomrb = %3.2e ', Lip, norm(b));
stoptol = 1e-6;
opts.stoptol = stoptol;
opts.Lip = Lip;
Ainput.A = A;
Ainput.Amap = @(x) Amap(x);
Ainput.ATmap = @(x) ATmap(x);
%% group structure
grpsize = 300; grpNUM = round(p/grpsize);
[A,G,ind,~] = getGroup(A,grpNUM);
%% solver
lam1_max = norm(A'*b,Inf);
lambda1 = lam1_max*(1e-3); lambda2 = lambda1;
c = [lambda1;lambda2];
[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,p,c,G,ind,opts);



