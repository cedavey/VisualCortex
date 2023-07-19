%% test distribution of presyn to postsyn connections
% sum occurences of postsyn relative indices from presyn neurons
plotRelativeInd = true ;
generateConns   = false;
pre='B'; post='C';
try
   outstruct2variables;
end

if plotRelativeInd
   sub = cell2mat(arrayfun(@(i) getRelativePosition(sz.(post),i,sz.(pre)),1:prod(sz.(post)),...
                         'uniformoutput',false)'); 
   ind = subv2ind(sz.(post),sub); 
   a   = accumarray(sub,ones(size(ind))); 
   figure; imagesc(a); axis image; colorbar;
end

if generateConns
   sz.(pre)  = [10 10];
   sz.(post) = [100 100];
   vect = cell(prod(sz.(post)));
   ind  = cell(prod(sz.(post)));
   for i=1:prod(sz.(post))
      % in matlab it's faster to allocate to a cell array than to an object (cuz
      % object needs permissions check etc every time you assign)
      vect{i} = initVisualCortexGaussianConnections(N.([pre post]),r.([pre post]),sz.(pre),sz.(post),i,[]); 
      ind{i}  = subv2ind(sz.(pre),vect{i}); % convert 2D location index to univariate index
   end
else
   % Evaluate an existing network
   if ~exist('conns','var')
      load(outfile);
      outstruct2variables;
   end
   ind  = conns.([pre post]);
   vect = synloc.([pre post]);
end

% sum connections from each presyn neuron
lc = prod(sz.(post)); % length postsyn indices
Nc = min([prod(sz.(pre)) 1000]);  % number of presyn to collate connections across
nc = cell(Nc,1);
for k=1:Nc
   p=arrayfun(@(i) ternaryOp(any(ind{i}==k),i,0),1:lc); 
   p(p==0)=[]; 
   nc{k}=p; 
end
l=cellfun(@length,nc); 
figure; plot(l);
title('Num cxs from each presyn neuron');
xlabel('Presyn neuron ID');
ylabel('Number of connections');

% Colour code B -> A connections
c=cell2mat(ntwkconfig.conns.([pre post])(:)'); 
figure; imagesc(c); colorbar; axis square; 
title('Matrix of presyn connection IDs');





