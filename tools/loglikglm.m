function [lp, G, H] = loglikglm(x,y,thta)
% function [lp, G, H] = loglik_glm(x,y,thta)
% calculate the negative log-likelihood of a canonical GLM.  
% 
% 
% lp- log-likelihood of a canonical GLM
% G - gradient of log-likelihood 
% H - Hessian of log-likelihood
%
% ONLY WORKS FOR POISSON DISTRIBUTED RESPONSES
% ADR-04-14

% sufficient statistic
ss = x'*y;

% projection of input data
eta = x*thta;
feta = exp(eta);

% log-liklihood
lp = -thta'*ss + sum(feta);

if nargout > 1
    % gradient
    G = -ss + x'*feta;
end
if nargout > 2
    H = -x'*spdiags(feta,0,length(feta),length(feta))*x;
end




end