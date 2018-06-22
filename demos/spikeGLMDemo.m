% Example script for fitting an "spiking" GLM model using the 
% maximum penalized expected log-likelihood estimator (MPELE) 
% to inform the initialization and optimization of an approximation to the 
% maximum likelihood estimator of the model. 
%
% adr 
% adr2110@gmail.com
%% set parameters and simulate from model
[Stim,tsp,ggsim] = spikeGLMsimulation;
% the hardcoded value in spikeGLMsimulation is to simulate spiking
% responses from 10 cells. 
cell2fit=1 ; % cell number we will fit
%% infer parameters using penalized maximum likelhood estimator (L1 penalty applied to coupling coefficients)

% create a fitting structure to work with functions from Johnathan Pillow's GLMspiketools package
gg0 = prepareSpikeGLMFit(Stim,tsp,size(ggsim.k,1),ggsim.dt,cell2fit);
mlet0 = tic;thetaMLE = MLGLMNET(Stim,gg0);mlet = toc(mlet0);
%% infer parameters using a fast penalized maximum likelhood estimator approximation by keeping stimlulus filter fixed at the maximum penalized expected log-liklihood (MPELE)

wnCov = @(z) z; % stimulus covariance function for white noise
melet0 = tic;thetaMLEApprox = MLApproxGLMNET(Stim,gg0,wnCov);melet = toc(melet0);

%% view exact solutions at best values 
plotEstimatorsSpikeGLM(ggsim,cell2fit,thetaMLE,thetaMLEApprox,gg0,mlet,melet)


