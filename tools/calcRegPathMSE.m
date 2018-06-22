function gtMSE = calcRegPathMSE(thetaMLE,thetaGT,gg0)
% gtMSE = calcRegPathMSE(thetaMLE,thetaGT,gg0)
% calculate the mean-squared error for simulated data for a series of 
% infered parameter values. Function is specialized to relate basis coefficient 
% estimates used in Johnathan Pillow's GLMspiketools package to "true" fully expanded parameters
% in the spiking GLM model
%
% INPUT :
% thetaMLE = p1 x q matrix of infered parameter values, where p1 is total
%            number of parameters used for inference, q is total number of inferences 
% thetaGT = p x 1 vector of ground truth parameter values. p is total number of 
%           parameters in the spiking model (p typically > p1).
%
% gg0  = structure specific to GLMSpiketool package, created using 
%        makeFittingStruct_GLM_adrMod
%
% OUTPUT 
% gtMSE = q x 1 vector of root mean squared error values 
numLambdas = size(thetaMLE,2);
gtMSE = zeros(numLambdas,1);
fix.kt = 0; fix.dc = 0; fix.ih = 0; fix.ih2 = zeros(size(gg0.ih2)); % binary variables 1 to fix component 0 to vary component

[nkt,nkx] = size(gg0.k); numcells = 1 + size(gg0.ih2,2);
for lindex = 1 : numLambdas
    bmle=thetaMLE(:,lindex);
    % populate gg0 sructure with current parametrs
    ggMLE = reinsertFitPrs_GLM_adrMod(gg0,bmle,fix);
    
    % expand basis functions to determine "full" filter parameters
    [~,gge] = expandFitPrs(ggMLE,bmle);
    bmleExpanded = [reshape(flipud(gge.k),nkt*nkx,1);gge.dc;reshape(gge.ih,length(gg0.iht)*numcells,1)];
    
    % root MSE 
    gtMSE(lindex) = norm(bmleExpanded-thetaGT);
end
end

