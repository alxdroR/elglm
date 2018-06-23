function [X,tsp,ggsim] = spikeGLMsimulation
%[X,y] = spikeGLMsimulation. Use Johnathan Pillow's GLMspiketools package
%        to simulate responses from a spiking GLM model
% OUTPUT : 
%        X  - n x p stimulus matrix 
%        tsp  - "spike" times from the LNP model 
%        ggsim - model parameters stored as a structure
% 
% adr

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a single-cell simulation structure 
global RefreshRate MAXSIZE;  
RefreshRate = 100; % Stimulus refresh rate (Stim frames per second)
MAXSIZE = .625e8; % Maximum amount to store in memory (1.25e8 ~ 1/2 GB)

DTsim = .01; % Bin size for simulating model & computing likelihood.
nkt = 15;  % Number of time bins in filter;
ggsim1 = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params


% make a spatio-temporal filter 
kt = ggsim1.k;  % Temporal filter

% make a spatial filter;
nkx = 5; 
xxk = [1:1:nkx]';
kx = 1./sqrt(2*pi*4).*exp(-(xxk-nkx/2).^2/5);
Filt = kt*kx'; % Make space-time separable filter
ggsim1.k = Filt./norm(Filt(:))*3; % Insert into simulation struct

ggsim1.dc = 1;

% create spike-history filter using the default 
% basis functions coded in the script for making fitting structures

ggGroundTruth = makeFittingStruct_GLM_adrMod(zeros(nkt,nkx),DTsim);
%[~,ihbas,ihbasis] = makeBasis_PostSpike(ggGroundTruth.ihbasprs,DTsim);
ihGroundTruth = [-10 -5 0 2 -2]';
ggsim1.ih = ggGroundTruth.ihbas*ihGroundTruth;
ggsim1.iht = ggGroundTruth.iht;


% use single-cell template to make coupled GLM structure for 2 cells
ggsim = makeSimStruct_GLMcpl(ggsim1,ggsim1); 

% increase the number of cells simulated -----
numcells = 10; % number of cells

% Add spatial filters and offsets for other cells
ggsim.k(:,:,3:numcells) = repmat(ggsim1.k(:,:),[1 1 numcells-2]);
ggsim.dc(3:numcells) = repmat(ggsim1.dc,[1 numcells-2]);

% create spike-history connectivity matrix
ntpsh = size(ggsim.ih,1);
IHAll = zeros(ntpsh,numcells,numcells);
% all self-history filters are identical
for i=1:numcells
    IHAll(:,i,i) = ggsim.ih(:,1,1);
end

% connections:
% connection basis function 
ggGroundTruth.ihbas2 = exp(-ggGroundTruth.iht*0.5);

% connections to cell 1 
nExcon = 1; % number of excitatory connections 
nIncon = nExcon; % number of inhibitory connections
nzero = numcells-1-nExcon-nIncon; % number of non-coupled filters

% spike-history basis coefficients for connections to cell 1
coef = [0.2 + (0.25-0.2)*rand(1,nExcon),...  % excitatory basis fnc coefficients
    -0.2 + (-0.25 + 0.2)*rand(1,nIncon), ... %  inhibitory basis fnc coefficients
    zeros(1,nzero)]; % non-coupled coefficients

IHAll(:,2:numcells,1) =  ggGroundTruth.ihbas2*coef;


% connections to cell 2
IHAll(:,1,2) = IHAll(:,2,1); % recipricol connection to cell 1;

% all other connections are zero
ggsim.ih = IHAll;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample stimulus from a Gaussian distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set number of samples determined by ratio of number of parameters/number of samples
pnratio = 0.04;
p =  (numcells-1) + ggGroundTruth.ihbasprs.ncols + 1 + numel(ggGroundTruth.kt);
N = round(p/pnratio);

swid = size(ggsim.k,2);% number of stimulus filter timebins and width (stimframes,pixels).
X = randn(N,swid);  %  Gaussian white noise stimulus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate GLM response and store response as spike times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsp = simGLM_adrMod(ggsim, X);  


end

