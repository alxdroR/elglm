function [X,y,thtaM] = LNPsimulation
%[X,y] = LNPsimulation. simulate from LNP model
% OUTPUT : 
%        X  - n x p stimulus matrix 
%        y  - "spike" times from the LNP model 
%        thta - model parameters
% 
% adr
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make temporal filter 
tau = 16;            % number of temporal filter bins
tauv = (0:tau-1)';
b1 = 0.03*tau;
b2 = 2*b1;
tthta = flipud(gampdf(tauv,6,b1) - gampdf(tauv,6,b2));

% make spatial filter
nx = 9;  % number of spatial bins
xv = (1:nx)';
xthta = normpdf(xv,nx/2,3);

thtaM = xthta*tthta'; % make space-time filter (theta vector arranged as a matrix)
thtaM = 0.3*thtaM/norm(thtaM);                      % normalize to theta 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample stimulus from a Gaussian distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100*nx*tau; % set number of samples equal to 100*number of parameters

x = randn(N,nx);  %  Gaussian white noise stimulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate poisson parameter and simulate from poisson distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = toeplitzblk(x,[x(1,:) sparse(1,(tau-1)*nx)]);
thta = fliplr(thtaM); 
eta = X*thta(:);
% if X is too large for memory, eta can also be constructed using 
% eta = stimConv(x,thta');
lambda = exp(eta);
y = poissrnd(lambda);



end

