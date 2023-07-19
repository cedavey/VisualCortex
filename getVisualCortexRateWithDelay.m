%  inrates = getVisualCortexRateWithDelay(rates,conninds,delayinds,psp,dt)
function inrates = getVisualCortexRateWithDelay(rates,conninds,delayinds,psp,dt)
   % Need to apply epsp for every postsyn neuron since there'll
   % be different propagation delay from pre to post for each.
   % Need to calculate matrix index for each synaptic input
   % using the appropriate delay (lag 0 input will be the last
   % entry, so higher delay means smaller index)
   
   % delta psp may be empty, but make delay 1 since we want last sample
   pspdelay  = ternaryOp(isempty(psp),1,length(psp));
   pre_sz    = size(rates);
   delayinds = ternaryOp(isempty(delayinds), zeros(1,pre_sz(2)),delayinds);
   
   % To allow recurrent connections we want inputs to be from last
   % sample, else if there are recurrent connections & we calculate inputs
   % using current time point, then when we calculate inputs going the
   % other way in the recurrent connections, they will get new inputs
   endind   = pre_sz(1) - delayinds;
   startind = endind - pspdelay + 1;
   inrates  = zeros(pspdelay,length(conninds));
   for j=1:length(conninds)
      inrates(:,j) = rates(startind(j):endind(j),conninds(j));
   end
   if ~isempty(psp)
      inrates = toVec(applyEPSP(inrates,psp,dt,false));
   end

end