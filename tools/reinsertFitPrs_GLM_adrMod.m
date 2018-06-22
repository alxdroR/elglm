function gg = reinsertFitPrs_GLM_adrMod(gg,prs,fix)
% gg = reinsertFitPrs_GLM_cond(gg,prs,indEXT);
%
% After fitting, reinsert params into param structure

global indEXT

[nt nx] = size(gg.kt);

if ~isempty(indEXT.kt)
    gg.kt = reshape(prs(indEXT.kt),nt,nx);
    gg.k = gg.ktbas*gg.kt;
end
if ~isempty(indEXT.dc)
    gg.dc = prs(indEXT.dc);
end
if ~isempty(indEXT.ih)
    gg.ih(:,1) = prs(indEXT.ih);
end

cells = 1:size(gg.ih2,2); % coupled cell indices
accells = (1-fix.ih2).*cells; % active cells
acind = find(accells); % active cell indices 
nacells = sum(1-fix.ih2); % number of active cells 
 
if ~isempty(indEXT.ih2)
        gg.ih2(:,acind) = reshape(prs(indEXT.ih2),size(gg.ih2,1),nacells);
end




