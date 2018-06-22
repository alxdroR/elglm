function [lp, G, H] = loglikglmL2(x,y,thta,k)
% function [lp, G, H] = loglik_glm(x,y,thta,k)
% calculate the negative log-likelihood of a canonical GLM with constant L2
% penalty paramerized by k.
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
lp = -thta'*ss + sum(feta) + 0.5*k*norm(thta)^2;

if nargout > 1
    % gradient
    G = -ss + x'*feta + k*thta;
end
if nargout > 2
    H = x'*spdiags(feta,0,length(feta),length(feta))*x + k*speye(length(thta));
end




end