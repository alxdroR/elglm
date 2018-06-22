function prs = extractFitPrs_GLM_adrMod(gg,fix,Stim)
% prs = extractFitPrs_GLM(gg,Stim,fix);
%
% Use Fix to determine the indices corresponding to fixed and variable paramters
% extract initial starting paramter

% make sure that number of coupled cells matches number of 
% cells specified in fix coupling term
nc = size(gg.ih2,2);
if nc ~= length(fix.ih2)
    if isempty(fix.ih2)
        error(['GLM structure,gg, has instructions to fit with ' num2str(nc) ' coupled cells but fix.ih2 does not specify which cells to hold fixed!?!?']);
    else
        error('# of fixed cells does not match number of cells in GLM structure')
     end
end

setupfitting_GLM_adrMod(gg,Stim,fix);  % Precompute quantities for optimization


if fix.kt
    keepind = fixc2ind(gg,fix,Stim);
else
     keepind = fixc2ind(gg,fix);
end
    
prs = [gg.kt(:); gg.dc; gg.ih(:); gg.ih2(:)];
prs = prs(keepind);
