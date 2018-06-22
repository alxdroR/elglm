% Example script for fitting an "LNP" model using the maximum expected log-likelihood estimator (MELE)
% or maximum penalized expected log-likelihood estimator (MPELE). The demo also
% shows how to use these two estimators as initializors for a quick
% approximation to the maximum likelihood estimator of the model
%
% adr 
% adr2110@gmail.com


%% set parameters and simulate from LNP model
[X,y,thtaM] = LNPsimulation;

%% Part 1 : Maximum likelihood 

%% infer parameters using maximum likelhood estimator (MLE)
mlet0 = tic; bmle = glmfit(X,y,'poisson','constant','off'); mlet = toc(mlet0);

%% infer parameters using fast maximum likelihood approximation using PCG initialized at the maximum expected log-liklihood (MELE)
% approximate MLE solution by initializing using MELE and then solving 
% the maximum likelihood optimization using preconditioned conjugate gradient (PCG) ascent with 
% the expected log-likelihood informing the precondioner (see paper for details)


maxNumIters = 3; % maximum number of PCG iterations
wnCov = @(z) z; % stimulus covariance function for white noise
melet0 = tic; [bapp,~,bmele]=MLApprox(X,y,wnCov,maxNumIters);melet = toc(melet0);
%% plot results and compare computation time with and without MELE initialization
plotEstimatorsLNP(thtaM,bmle,bmele,bapp,mlet,melet);
set(gcf,'Name','No likelihood penalty function');


%% Part 2 : Penalized Maximum Likelihood 
%% infer parameters using L2-penalized maximum likelhood estimator (MLE)
lambda = 0.5;
mlet0 = tic; bmle = lassoglm(X,y,'poisson','Alpha',1e-8,'Lambda',lambda,'Standardize',false); mlet = toc(mlet0);
%% infer parameters using a fast L2-penalized MLE approximation via PCG initialized at the maximum penalized expected log-liklihood (MPELE)

pwnCov = @(z,n,lambda) z/(n+lambda); % penalized covariance function for white noise
melet0 = tic;[bapp,~,bmele]=MLApprox(X,y,pwnCov,maxNumIters,'L2Lambda',lambda);melet = toc(melet0);
%% plot results and compare computation time with and without MELE initialization
figure;plotEstimatorsLNP(thtaM,bmle,bmele,bapp,mlet,melet);
set(gcf,'Name','L2-penalized likelihood function');


