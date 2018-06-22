function [prse gge] = expandFitPrs(gg,prs)
% [prse gge] = expandFitPrs(gg,prs)
IHB = gg.ihbas; % spike-history basis functions
IHB2 = gg.ihbas2;

KTB = flipud(gg.ktbas);% stimulus history basis functions

[nbf,swid] = size(gg.kt); % proper dimensions for stimulus-related paramteres
nkp = nbf*swid;% number of stimulus-related paramters being fit
nlags = size(KTB,1);  % number of stimulus lags
nCoupled = size(gg.ih2,2); % number of cells

if isfield(gg,'ihfix')
    Ihfix = gg.ihfix;
else
    Ihfix = 0;
end
for i=1:size(prs,2)
     K = KTB*reshape(prs(1:nkp,i),nbf,swid);
     b = prs(nkp+1,i);
     
     iht = reshape(prs(nkp+2:nkp+1+size(IHB,2),i),size(IHB,2),1);
     iht2 = reshape(prs(nkp+2+size(IHB,2):end,i),size(IHB2,2),nCoupled);
     SH = IHB*iht + Ihfix;
     SH2 = IHB2*iht2;
     
     % couldn't the stuff below be replace by 
% prse(:,i) = [K(:);b;SH(:)] ? 
     prse(:,i) = [K(:);b;SH(:);SH2(:)];
%      for k=2:ncells
%        p0 = [p0;SH(:,k)];
%      end

%      prse(:,i) = p0;
gge.ih(:,:,i) = [SH SH2];
end

% stimulus history coefficients
nsdim  = nlags*swid; % number of stimulus-related dimensions
gge.k = reshape(prse(1:nsdim,:),nlags,swid,size(prs,2));
gge.dc = prse(nsdim+1,:);

