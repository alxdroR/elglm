function plotEstimatorsSpikeGLM(ggsim,cell2fit,thetaMLE,thetaMLEApprox,gg0,mlet,melet)

% compute mean-squared error of fits
%% view exact solution using regularization value that minimizes se
thtaM = ggsim2theta(ggsim,cell2fit); % ground truth parameter
mseExact = calcRegPathMSE(thetaMLE,thtaM,gg0);
mseApprox = calcRegPathMSE(thetaMLEApprox,thtaM,gg0);

% determine parameter with minimum squared error
[bmle,ggMLE] = minGTMSE(thetaMLE,mseExact,gg0);
[bapp,ggMLEApprox] = minGTMSE(thetaMLEApprox,mseApprox,gg0);

% covert parameters into a structure relevant for viewing spiking GLM
% functions
[~,gge] = expandFitPrs(ggMLE,bmle);
[~,ggeApprox] = expandFitPrs(ggMLEApprox,bapp);

nkx = size(ggsim.k,2);
xxk = [1:1:nkx]';

nkt = size(ggsim.k,1);
ttk = [-nkt+1:0]';
numcells = size(ggsim.ih,2);


figure; numrows = 2; numcols = 5;

subplot(numrows,numcols,[1 2 numcols+1 numcols+2]); plot(bmle,bapp,'.'); title('parameter values at min \lambda');xlabel(['MLE']);
ylabel({'MLE approximation computed ' [num2str(mlet/melet) ' times faster']});
hold on; plot([min(bmle,bapp),max(bmle,bapp)],[min(bmle,bapp),max(bmle,bapp)],'k--')

subplot(numrows,numcols,[3]);imagesc(xxk,ttk,ggsim.k(:,:,cell2fit)); title({'ground truth(gt)' ['offset = ' num2str(ggsim.dc(cell2fit))]});
box off; ylabel('stimulus kernal time');
subplot(numrows,numcols,[4]);imagesc(xxk,ttk,flipud(gge.k)); title(['MLE-' num2str(gge.dc)]); box off
xlabel('stimulus filter space');
subplot(numrows,numcols,[5]);imagesc(xxk,ttk,flipud(ggeApprox.k)); title(['MELE-' num2str(ggeApprox.dc)]);

subplot(numrows,numcols,numcols+[3 4]);plot(ggsim.iht, squeeze(ggsim.ih(:,cell2fit,cell2fit)),'b--',...
    ggMLE.iht, gge.ih(:,1),'r',...
    ggMLEApprox.iht, ggeApprox.ih(:,1),'b',...
    ggsim.iht, ggsim.iht*0, 'k--');title('self-history terms');legend('gt','MLE','MLE w fixed SF')
xlim([0 10]);ylim([-10 2])


subplot(numrows,numcols,numcols+5);plot(1:(numcells-1), squeeze(ggsim.ih(1,setdiff(1:numcells,cell2fit),cell2fit)),'b--',...
    1:(numcells-1), gge.ih(1,setdiff(1:numcells,cell2fit)),'r',...
    1:(numcells-1), ggeApprox.ih(1,setdiff(1:numcells,cell2fit)),'b',...
    1:(numcells-1), zeros(numcells-1,1), 'k--');title('coupling terms')
xlim([0 7])
end

