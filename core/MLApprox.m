function [THTA,ftrain,bmele] = MLApprox(X,y,meleCov,maxNumIters,varargin)
options = struct('L2Lambda',[],'method','lbfgs','precFunc',[]);
options = parseNameValueoptions(options,varargin{:});


% Step 1: Construct the marginalized stimulus filter estimator for an LNP model.  
% -------------------------------------------------------------------------   
if isempty(options.L2Lambda)
    bmele = mele(X,y,meleCov,'poisson');
else
     bmele = mpele(X,y,meleCov,options.L2Lambda,'poisson');
end

%-----------------
% Step 2: Run a small number of PCG iterations down the log-likelihood-----------------------------------------------------------


% PCG optimization parameters 
minFuncOpts = optimset('MaxIter',maxNumIters,'Display','final');minFuncOpts.Method = options.method;minFuncOpts.c2 = 0.9;
if ~isempty(options.precFunc)
    minFuncOpts.precFunc = options.precFunc;
end

if isempty(options.L2Lambda)
    loss = @(t) loglikglm(X,y,t);
else
    loss = @(t) loglikglmL2(X,y,t,options.L2Lambda);
end

if maxNumIters > 0
    [THTA,ftrain] = minFunc(loss,bmele,minFuncOpts);           
elseif maxNumIters == 0 
    if nargout > 1
        ftrain = loss(bmele);
    end
    THTA = bmele;
else
    error('ELE:MLApproxPCG:maxNumItersProblem','maxNumIters must be >= 0');
end


