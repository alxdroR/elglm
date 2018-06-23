function THTA = MLApproxGLMNET(Stim,gg,meleCov,varargin)
% THTA = MLApproxGLMNET(Stim,gg,meleCov,varargin)
% MLApproxPCG
% Fast approximation of the maximum likelihood estimate of a GLM spike
% train model using the maximum expected log-likelihood to initialize the "stimulus" parameters
%


options = struct('L2Lambda',[]);
options = parseNameValueoptions(options,varargin{:});

global RefreshRate

% Step 1: Construct the marginalized stimulus filter estimator for an LNP model.  
% -------------------------------------------------------------------------
if isempty(options.L2Lambda)
    thetaMELE = mele(Stim,gg,meleCov,'poisson');
else
    thetaMELE = mpele(Stim,gg,meleCov,options.L2Lambda,'poisson');
end


% reshape mpele
nkt = size(gg.k,1);stimP = size(Stim,2);
gg.k = reshape(thetaMELE,nkt,stimP);
gg.kt = pinv(gg.ktbas)*gg.k;

% ------ construct a design matrix that holds stimulus parameters fixed to MELE
[X,Y]=setup_GLMNET(gg,Stim,1);


% glmnet parameters -- set penalty for stimulus scalar multiple and spike self-history
% terms to 0 (no regularization of these terms)
opts = glmnetSet;opts.penalty_factor = ones(1,size(X,2));opts.penalty_factor(1:1+numel(gg.ih)) = 0;
fit = glmnet(X,Y,'binomial',opts);

% extract parameters for the fit GLMNET structure
pMLEApprox = fit.dim(1)+1; % number fit + offset 
THTA = zeros(length(gg.kt(:))-1+pMLEApprox,fit.dim(2));
for k=1:fit.dim(2)
    THTA(:,k) = [gg.kt(:)*fit.beta(1,k)*-1;fit.a0(k)*-1-log(gg.dt/RefreshRate);fit.beta(2:end,k)*-1];
end

