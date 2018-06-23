function gg = makeFittingStruct_GLM_adrMod(sta,DTsim,glmstruct,cellnumToFit);
% gg = makeFittingStruct_GLM(sta,DTsim,glmstruct,cellnumToFit);
%
% Initialize parameter structure for fitting of GLM model,
% normal parametrization of stim kernel
%
% Inputs:  sta = initial guess at kernel (or use all zeros if unknown)

% Set up structure
gg.k = [];
gg.dc = 0;
gg.ih = [];
gg.iht = [];
gg.ihbas = [];
gg.ihbas2 = [];
gg.ihbasprs = [];
gg.kt = [];
gg.ktbas = [];
gg.kbasprs = [];
gg.tsp = [];
gg.tspi = [];
gg.dt = DTsim;
gg.ih2 = [];
gg.ihbas2 = [];
gg.ihbasprs2 = [];
gg.tsp2 = [];
gg.couplednums = [];

% === Make temporal basis for stimulus filter =======================
[nkt,nkx] = size(sta);
% % ----- Set up temporal basis for stimulus kernel -----------
kbasprs.neye = 5; % Number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % Number of raised-cosine vectors to use  
kbasprs.kpeaks = [0 round(nkt/3)];  % Position of first and last bump (relative to identity bumps)
kbasprs.b = 3; % Offset for nonlinear scaling (larger -> more linear)
if nkt>1
ktbas = makeBasis_StimKernel(kbasprs,nkt);
else
    ktbas=1;
end
gg.ktbas = ktbas;
gg.kbasprs = kbasprs;

% ======================================================================
% Set up basis for post-spike kernel
global RefreshRate
onems = .001*RefreshRate;% 1 ms
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [onems 3];  % Peak location for first and last vectors
ihbasprs.b = .4;  % How nonlinear to make spacings
ihbasprs.absref = onems; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);

gg.iht = iht;
% Orthogonalization seems to put the refractory basis fnc at index 4
%gg.ihbas = [ihbas(:,1:3) ihbas(:,5)];
gg.ihbas = ihbasis;
gg.ihbasprs = ihbasprs;
gg.ih = zeros(size(gg.ihbas,2),1);
%gg.ihfix = ihbasis(:,1)*-10;


% coupling terms
ihbasprs2.ncols = 4;  % Number of basis vectors for post-spike kernel
ihbasprs2.hpeaks = [onems 3];  % Peak location for first and last vectors
ihbasprs2.b = .4;  % How nonlinear to make spacings
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs2,DTsim);
gg.ihbas2 = ihbasis;
gg.ihbasprs2 = ihbasprs2;
gg.ih2 = zeros(size(gg.ihbas2,2),1);


% % ==================================================================
% set up initial K params
gg.kt = inv(gg.ktbas'*gg.ktbas)*gg.ktbas'*sta;
gg.k = gg.ktbas*gg.kt;

% % ==================================================================
% If full param struct passed in, match other params as well
if (nargin >= 3) 
    gg.dc = glmstruct.dc;
    gg.ih = glmstruct.ih;
    gg.iht = glmstruct.iht;

    %---Extract correct ih basis params, if present----
    if isfield(glmstruct, 'ihbasprs')
        if ~isempty(glmstruct.ihbasprs);
            ihbasprs = glmstruct.ihbasprs;
            [iht,ihbas ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);
            % Orthogonalization seems to put the refractory basis fnc at index 4
             ihbas = [ihbas(:,1:3) ihbas(:,5)];
             gg.ihfix = ihbasis(:,1)*-10;
   
            % -- Do some error-checking ----
            if length(iht) ~= length(glmstruct.iht)
                error('mismatch between iht and h-kernel params ihbasprs');
            end
            if size(glmstruct.ih,2)>1 & (nargin < 4)
                error('multi-cell glm struct passed in without cell # to fit');
            end

            %--- Put glmstruct params into gg ----
            gg.iht = glmstruct.iht;
            gg.ihbas = ihbas;
            gg.ihbasprs = ihbasprs;
            if nargin == 3  % single-cell only
                gg.ih = inv(ihbas'*ihbas)*ihbas'*glmstruct.ih;
                gg.dc = glmstruct.dc;
            else % mulitcell-cell
                ncells = size(glmstruct.ih,2);
                ih0 = glmstruct.ih(:,:,cellnumToFit);
                ih1 = ih0(:,cellnumToFit);
                ih2 = ih0(:,setdiff(1:ncells,cellnumToFit));
                gg.ih = inv(ihbas'*ihbas)*ihbas'*[ih1 ih2];
                gg.dc = glmstruct.dc(cellnumToFit);
            end

        end
    end


end

