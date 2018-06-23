function gg0 = prepareSpikeGLMFit(Stim,tsp,nkt,dt,cell2fit)
% prepareSpikeGLMFit  - Pre-compute necessary quantities for fitting spiking model



% Compute STA and use as initial guess for k
sta0 = simpleSTC(Stim,tsp{cell2fit},nkt);
sta = reshape(sta0,nkt,[]);


%  Initialize params for fitting --------------
gg0 = makeFittingStruct_GLM_adrMod(sta,dt);  % Initialize params for fitting struct
gg0.tspi = 1; % set 1st spike to use for computing likelihood to 1

% PUT RELEVANT SPIKES INTO FITTING STRUCTURE
gg0.tsp = tsp{cell2fit};   % cell 1 spike times (vector)

numcells = length(tsp);
gg0.tsp2 = tsp(setdiff(1:numcells,cell2fit));

gg0.ih = zeros(length(gg0.ih),1);
gg0.ih2 = zeros(1,numcells-1);
gg0.ihbas2 = exp(-gg0.iht*0.5);


end

