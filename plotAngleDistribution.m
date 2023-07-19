%% Plot histogram retina angle from layer B
% Plots angles from retina that connect through 2 layers to end in layer C,
% or whatever layer the user specifies. It's assumed to go from A->B, and
% then from B->' '. 
% The network and output data is all assumed to be in outstruct.
nN         = 16;
layer_name = 'L';
randInd    = true;
connlabels = ['AB'; 'BL'; 'LC'];
% connlabels = ['AB'; 'BC'];
fopts      = {'fontsize',20,'fontweight','bold'};

if randInd
   nind   = sort(randi(prod(sz.(layer_name)),nN,1));
else
   maxDisplay = ceil(sqrt(nN));
   N = min([maxDisplay outstruct.layerconfig.sz.(layer_name)(1)]);
   M = min([maxDisplay outstruct.layerconfig.sz.(layer_name)(2)]);
   if outstruct.layerconfig.sz.(layer_name)>maxDisplay
      dN    = outstruct.layerconfig.sz.(layer_name)(1) / maxDisplay;
      dM    = outstruct.layerconfig.sz.(layer_name)(2) / maxDisplay;
      dNM   = round(dN*dM);
      nind  = dNM:dNM:prod(outstruct.layerconfig.sz.(layer_name));
   else
      nind  = 1:prod(outstruct.layerconfig.sz.(layer_name));
   end
   if length(nind)>nN, nind(nN+1:end)=[]; end
end

ntheta = cell(nN,1);

figure; 
nbins = 30;
% for each connection, for each postsyn neuron, extract its presyn conns
nr = ceil(sqrt(nN)); nc = ceil(nN/nr);
for ni=1:length(nind)
   ntheta{ni} = [];
   cx = cell(size(connlabels,1),1);
   % connection indices of presyn neurons in last network layer connection 
   cx{end}{ni} = outstruct.ntwkconfig.conns.(connlabels(end,:)){nind(ni)};
%    for cci=(size(connlabels,1)):-1:1
   for cci=(size(connlabels,1)-1):-1:1
      % get the number of connection in the previous network layer connection
      cx{cci} = cell(outstruct.layerconfig.N.(connlabels(cci,:)),1);
      % for each of the presyn neurons in the current cx, get all of the
      % indices of the previous layer connection
      for ci=1:length(cx{cci+1}{ni})
         cx{cci}{ci} = outstruct.ntwkconfig.conns.(connlabels(cci,:)){cx{cci+1}{ni}(ci)};
%          cx{cci}{ci} = outstruct.ntwkconfig.conns.(connlabels(cci,:)){cx{cci+1}(ci)};
%          cx{cci} = outstruct.ntwkconfig.conns.(connlabels(cci,:)){cx{cci+1}(ci)};
         if cci==1
            ntheta{ni} = [ntheta{ni}; toVec(outstruct.layerconfig.RF.theta(cx{cci}{ci}))];
         end
      end
   end
   subplot(nr,nc,ni)
   hist(ntheta{ni},nbins); 
   xlim([0 180]);
   title(num2str(ind2subv(outstruct.layerconfig.sz.(layer_name),nind(ni))))
end

figure; 
Ntheta = cell2mat(ntheta');
hist(Ntheta(:),cell2mat(outstruct.layerconfig.RF.angleSet));
xlabel('\theta',fopts{:});
xlim([0 180]);


