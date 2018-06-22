function [thta,ss,normfac] = mele(x,y,C,distr,varargin)
% mele Quickly fit a canonical generalized linear model, using the
% maximum expected log-liklihood estimator.
%
%   T = mele(X,Y,C,DISTR) uses the marginalized expected log-likelihood to 
%   fit a generalized linear model using the predictor
%   matrix X, response Y, the known covariance matrix, C, of X.
%   DISTR specifies the distribution of the responses.  The result is T, a vector
%   of coefficient estimates fit using the canonical link associated with DISTR. 
%   DISTR can be 'normal' (default) or 'poisson'
%
%
%   X is a matrix with rows corresponding to observations, and columns to
%   predictor variables.  mele automatically includes a constant term in the
%   model (do not enter a column of ones directly into X).  Y is a vector of
%   response values.  It is assumed that X is drawn from a known distribution with 
%   covariance matrix C.  C is either a square symmetric matrix with
%   dimensions equal to the number of predictors OR a function of a single argument
%   C(b) that quickly determines the vector 'a' that solves the system of linear
%   equations Ca = b.
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

if nargin < 3
    error('ELE:mele:TooFewInputs','TooFewInputs');
end

if nargin < 4 || isempty(distr)
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
        end
    otherwise
        error('ELE:mele:BadDistribution','BadDistribution');
end


% compute the sufficient statistic required for both Poisson and normal
% solutions
if ~isstruct(y)
    ss = x'*y;
    ssn = ss/normfac;
else
    nkt = size(y.k,1);
    sta0 = simpleSTC(x,y.tsp,nkt);
    ssn = sta0(:);
end

% multiply the sufficient statisic by the inverse covariance matrix
if isnumeric(C)
    thta = C\ssn;
else
    thta = C(ssn);
end


