function [theta,ggmin] = minGTMSE(thetaMatrix,mse,gg0)
% [theta,ggmin] = minGTMSE(thetaMatrix,mse,gg0)
% return the parameters at the value of minimum mean-squared error given a series of 
% infered parameter values and associated mean-squared errors. Function also returns
% parameters as a structure specific to Johnathan Pillow's GLMspiketools package 
%
% INPUT :
% thetaMatrix = p1 x q matrix of infered parameter values, where p1 is total
%            number of parameters used for inference, q is total number of inferences 
% mse = q x 1 vector of root mean squared error values 
% 
% gg0  = structure specific to GLMSpiketool package, created using 
%        makeFittingStruct_GLM_adrMod
%
% OUTPUT 
% theta = p1 x 1 vector of inferred parameter values at the index where mse
% is lowest 
% ggmin = GLMSpiketool structure populated with values from theta
[~,minindex] = min(mse);
theta = thetaMatrix(:,minindex);

% populate gg0 sructure with current parametrs
fix.kt = 0; fix.dc = 0; fix.ih = 0; fix.ih2 = zeros(size(gg0.ih2)); % binary variables 1 to fix component 0 to vary component
ggmin = reinsertFitPrs_GLM_adrMod(gg0,theta,fix);
end

