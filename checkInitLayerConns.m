%% Get A->C cxs and angles

% Assumes that data comes from outstruct
synloc = outstruct.ntwkconfig.synloc;
conns  = outstruct.ntwkconfig.conns;
sz     = outstruct.layerconfig.sz;

%% Plot histogram retina angle from layer B
nN     = 9;
% nind   = sort(randi(prod(sz.C),nN,1));
nind   = 1:prod(sz.C);
ntheta = cell(nN,1);
sxAB   = cell(nN,1);

% for each layer C neuron of interest, extract B->C connections to it, and
% then work backwards to extract A->B connections, and the angle associated
% with each of those connections.
figure; hold on; 

for ni=1:length(nind)
   ntheta{ni} = [];
   sxAB{ni}   = [];
   cxBC       = conns.BC{ni};
   for ci=1:length(cxBC)
      sxAB{ni}   = [sxAB{ni}; synloc.AB{cxBC(ci)}];
%       ntheta{ni} = [ntheta{ni}; toVec(RF.theta(cxAB))];
%       sxAB{ni}   = [sxAB{ni}; synloc.AB{cxAB}];
   end
   plot(sxAB{ni}(:,1),sxAB{ni}(:,2),'o'); 
   xlim([0 sz.A(1)]);
   ylim([0 sz.A(2)]);
   pause(1);
%    clf(gcf,'reset');
end



