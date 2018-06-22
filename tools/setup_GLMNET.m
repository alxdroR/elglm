function [X,Y] = setup_GLMNET(gg,Stim,fixs)

% Step 2: Keeping stimulus current fixed, update self-history terms and then coupling terms
% -------------------------------------------------------------------------
% specify pre-conditioner
[~,nc] = size(gg.ih2);% number of sh coefficients and number of coupled cells


% SELF-HISTORY TERM ASSUMING FIXED STIMULUS CURRENT

fix.kt = 0; fix.dc = 0; fix.ih = 0; fix.ih2 = ones(1,nc);
extractFitPrs_GLM_adrMod(gg,fix,Stim);


global MSTM MMntrp OPRS SPNDS SPNDS2
if fixs
    SS = MSTM*gg.kt(:);
    S = MMntrp*SS; % Optimize gain and offset
else
    S = MMntrp*MSTM;
end
clear MSTM MMntrp 

% determine fixed refractory period term
iwin = OPRS.ichunk([1 2]); % abs fine bin win
iwin1 = [iwin(1)+1,iwin(2)];
ilen = diff(iwin);


% CONSTRUCT SELF-HISTORY COVARIATES (CONVOLVE SPIKES WITH BASIS FNCS)
Sh = spikeconv_mex(SPNDS, OPRS.ihbas, iwin1);


% DESIGN MATRIX -- offset is already included in glmnet
X = [S Sh];
clear S Sh
for cellind = 1:nc
    Ih0 = spikeconv_mex(SPNDS2{cellind}, OPRS.ihbas2, iwin1);
    X = [X Ih0];
end

% hack for using GLMNET wrapper to avoid underflow error when p << 1
Y = 2*ones(ilen,1); % no response is labeled 2
Y(SPNDS(SPNDS<=ilen))=1;

