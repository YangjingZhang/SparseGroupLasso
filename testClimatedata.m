%% Test on climate data set
clear; rng('default');
addpath(genpath(pwd));
probname = ['climatedata',filesep,'Climatedata.mat'];
load(probname);
%% dimension, group structure, options
[n,p] = size(A);
G = [1:p];
ind = zeros(3,10512);
ind(1,:) = 7*[0,1:10512-1] + 1;
ind(2,:) = 7*[1:10512];
ind(3,:) = sqrt(ind(2,:) - ind(1,:));
AATmap = @(x) A*(A'*x);
eigsopt.issym = 1;
Lip = eigs(AATmap,length(b),1,'LA',eigsopt);
fprintf('\n Lip const = %3.2e, nomrb = %3.2e ', Lip, norm(b));
tol = 1e-8;
opts.stoptol = tol;
opts.stopopt = 4; %stop if dualfeas<tol and pobj-dobj<tol*norm(b)^2
opts.Lip = Lip;
opts.printyes = 1;
Ainput.A = A;
Ainput.Amap = @(x) Amap(x);
Ainput.ATmap = @(x) ATmap(x);
%% solver
lambda1 = 100; lambda2 = 100;
c = [lambda1;lambda2];
[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,p,c,G,ind,opts);
%% plot active groups  in predicting Dakar's air temperature
xx = reshape(x,7,10512);
groupnorm = sum(xx.*xx).^0.5;
CurrentSubs = [31;138];
geoplot_climate(CurrentSubs,groupnorm);