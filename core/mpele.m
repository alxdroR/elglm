function [thta] = mpele(x,y,C,k,distr)
% MPELE 
% fit a canonical generalized linear model, using the
% L2-norm penalized marginalized expected log-liklihood approximation.
%
%   T = MPELE(X,Y,C,K) uses the marginalized expected log-likelihood 
%   and ridge parameter K to fit a generalized linear model using the Nxp predictor
%   matrix X, Nx1 response vector Y, and the known covariance matrix, C, X.  
%   If K is a scalar the result is T, a p-dmensional column vector.  If K has m elements 
%   T is pxm.  The fit assumes normally distributed responses and uses the canonical link.
%
%   C is either a square symmetric matrix with
%   dimensions equal to the number of predictors OR a function of three
%   arguments
%   C(b,n,k) that quickly determines the vector 'a' that solves the system of linear
%   equations (C*n+keye(p))a = b.  
%
%   T = MPELE(X,Y,C,K,DISTR) can be used to specify poisson distributed responses
%   by setting the string DISTR to 'poisson'. 
%
%   For certain applications it is desirable to use Johnathan Pillow's GLMspiketools package
%   to compute the maximum expected log-likelihood estimator. MELE can do
%   this if one uses the call: 
% 
%    T = mele(X,gg,C,'poisson') 
%   where gg is a structure of parameters specific to the Pillow package
%   (see package help files/testscripts for more details). If the second
%   argument is a structure, mele will assume it has the fields given by
%   the Pillow package. This only works with poisson distributed data

% adr 
% 2013

if nargin < 4
    error('ELE:MPELE:TooFewInputs','TooFewInputs');
end

if nargin < 5 || isempty(distr)
    distr = 'normal';
else
    distr = lower(distr);
end


% number of observations and number of coefficients
[n,p] = size(x);

% distribution information
switch distr
case 'normal'   
    normfac = n; % amount we multiply covariance matrix by in determining mele solution
    case 'poisson'
        if ~isstruct(y)
            if any(y < 0)
                error('ELE:mele:BadPoissonData','BadPoissonData');
            end
            normfac = sum(y);
        else
            normfac = length(y.tsp);
        end
    otherwise
        error('ELE:MPELE:BadDistribution','BadDistribution');
end

m = length(k); 


% compute the sufficient statistic required for both Poisson and normal
% solutions
if ~isstruct(y)
    thta = zeros(p,m);
    ss = x'*y;
else
    nkt = size(y.k,1);
    thta = zeros(p*nkt,m);
    sta0 = simpleSTC(x,y.tsp,nkt);
    ss = normfac*sta0(:);
end


% compute solution with regularization
for i = 1:m
    % add regularization to covariance matrix C
    if isnumeric(C)
        Cr = C*normfac + k(i)*eye(p);
        thta(:,i) = Cr\ss;
    else
        thta(:,i) = C(ss,normfac,k(i));
    end
end




