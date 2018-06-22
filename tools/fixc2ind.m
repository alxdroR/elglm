function keepind  = fixc2ind(gg,fix,Stim)
% convert the cell, fix, of boolean variables to an array of indices, keepind, 
% indicating which elements to keep from the theta vector containing all 
% possible elemtents.  Extractind yields the indices of the saved parameters 
% in the new theta vector that doesn't contain the paramaters that are fixed.
global OPRS indEXT 

if isfield(gg,'ihfix')
    OPRS.ihfix = [0;gg.ihfix(1:end-1)];
end


nh = OPRS.nh;
nh2 = OPRS.nh2;
nktot =  numel(gg.kt);
ncells = 1+size(gg.ih2,2); % shouldn't be OPRS.ncpl bc this number only reflects number 
                         % of active cells (which is important for memory
                         % mngment reasons).

% Initialized global vector used for appropriate variable selection in log-likelihood
% code
indEXT.kt = []; indEXT.dc = []; indEXT.ih = []; indEXT.ih2 =[];
indEXT.ikprs = []; indEXT.idc = []; indEXT.iihprs = []; 

% Initialize variable indices of the full parameter, without fixing any elements, constructed in extractFitPrs_GLM.m
indKT = 1:nktot;
indDC = nktot+1;
indih = nktot+2:nktot+nh+1;
        
prsoutidx = 0; % counts how many variables are being varied 
% If stimulus filter isn't fixed during optimization
if(~fix.kt)
    prsoutidx = nktot;
else
    indKT = [];
    % fix stimulus current using filter stored in gg
    if any(gg.kt(:))
        indEXT.ikprs = gg.kt(:);
    end
end
indEXT.kt = indKT; % Store indices related to stimulus filter parameters

% If dc offset isn't fixed during optimization
if(~fix.dc)
    prsoutidx = prsoutidx+1;
    indEXT.dc = prsoutidx;
else
    indDC = [];
    % fix dc current
    indEXT.idc = gg.dc;
    indEXT.dc = [];
end

% If spike-history filter isn't fixed during optimization
if(~fix.ih)
    indEXT.ih = prsoutidx+1:prsoutidx+nh;
    prsoutidx = indEXT.ih(end);
    
    % columns we will extract from the spike-convolved basis-fnc matrix
    tis = 1:nh; % self term
    
else
    %if any(gg.ih)
    %    indEXT.iihprs = [];
    %end
    indih = [];
    tis = [];
    % if any(gg.ih(:,1)) 
         indEXT.iihprs = gg.ih(:,1);
    % end
end

indih2 = [];
indEXT.ih2 = [];
tic = [];
stoffset = nktot+1+nh;
indEXT.active = zeros(1,sum(1-fix.ih2));
indEXT.inactive = 0;


cnt =1; cnti = 1;
for q = 1:(ncells-1)
    % If coupling terms are not fixed during optimization
    if(~fix.ih2(q))
        indih2 = [indih2 stoffset+1+(q-1)*nh2:stoffset+q*nh2];
        indEXT.active(cnt) = q;
        indEXT.ih2 = [indEXT.ih2 prsoutidx+1:prsoutidx+nh2];
        % columns we will extract from the spike-convolved basis-fnc matrix
        tic = [tic nh+1+(q-1)*nh2:nh2+q*nh2];
        prsoutidx = indEXT.ih2(end);
        cnt = cnt+1;
    else
        
        if any(gg.ih2(:,q))
             % if the filter isn't a string of zeros add the fixed coupled-history currents
            if ncells > 1
                     indEXT.iihprs = [indEXT.iihprs; gg.ih2(:,q)];
                     %indEXT.inactive = [indEXT.inactive, q]; 
                     indEXT.inactive(cnti) =  q; 
                    cnti = cnti+1;
            end
        end
    end
end


keepind = [indKT indDC indih indih2];

