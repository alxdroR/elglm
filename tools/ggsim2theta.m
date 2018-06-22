function thtaM = ggsim2theta(ggsim,cell2fit)
% ggsim2theta. Convert parameters stored in a simulation structure-ggsim- for cell at index cell2fit into a vector
% of parameters, thtaM. Simulation structure is from makeSimStruct_GLM

nkt = size(ggsim.k,1);
nkx = size(ggsim.k,2);
numcells = size(ggsim.ih,2);

thtaM = [reshape(squeeze(ggsim.k(:,:,cell2fit)),nkt*nkx,1);ggsim.dc(cell2fit);...
    squeeze(ggsim.ih(:,cell2fit,cell2fit))...
    ;reshape(squeeze(ggsim.ih(:,setdiff(1:numcells,cell2fit),cell2fit)),length(ggsim.iht)*(numcells-1),1)];

end

