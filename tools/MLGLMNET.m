function THTA = MLGLMNET(Stim,gg)

global RefreshRate

% create the design matrix to work iwth GLMNet
[X,Y]=setup_GLMNET(gg,Stim,0);

% glmnet parameters -- set penalty for stimulus and spike self-history
% terms to 0 (no regularization of these terms)
opts = glmnetSet;opts.penalty_factor = ones(1,size(X,2));opts.penalty_factor(1:numel(gg.kt)+numel(gg.ih)) = 0;options.maxit = 1e3;
fitGLM = glmnet(X,Y,'binomial',opts);

% extract parameters for the fit GLMNET structure
THTA = zeros(fitGLM.dim(1)+1,fitGLM.dim(2));
for k=1:fitGLM.dim(2)
    THTA(:,k) = [-1*fitGLM.beta(1:numel(gg.kt),k);-1*fitGLM.a0(k)-log(gg.dt/RefreshRate);-1*fitGLM.beta(numel(gg.kt)+1:end,k)];
end
end

