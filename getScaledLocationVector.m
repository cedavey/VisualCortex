%     ind_rel = getScaledLocationVector(ind, toSz, fromSz)
% 
% Get the relative position of the postsynaptic neuron, with the 2 layers 
% centred relative to one another.
%  [ind_rel] = getGroupsLocationArray(sz,index,fromSz)
% Inputs:
%  sz      - size of postsynaptic layer that the neuron resides in
%  index   - index of current neuron into the layer (scalar)
%  fromSz  - size of presynaptic layer if relative position required
% Outputs:
%  ind_rel - gets the location in presynaptic coordinates assuming 
%            that the 2 layers are centred around each other, & 
%            stretched so that the end points meet
%
%  See also: getGroupsLocationArray
function ind_rel = getScaledLocationVector(ind, fromSz, toSz)
   if size(toSz,2) ~= size(ind,2) || size(toSz,2)~=size(fromSz,2)
      fprintf('\ngetScaledLocationVector - size mismatch between layers or indices, exiting...\n\n');
      ind_rel = [];
      return;
   end
   % If the 'from' size is different to the 'to' size, then everything has
   % to be scaled, because the 'from' neuron that is directly 'above' the
   % 'to' neuron needs to be centred so that the 'from' layer is centred
   % directly above the 'to' layer. So calculate the scale of 'from' layer
   % to 'to' layer, & then shift so that middle points are aligned
   numInd     = size(ind,1);
   scale      = repmat((toSz)./(fromSz),[numInd 1]); % (fromSz-1)./(sz-1);
   centreFrom = repmat((fromSz-1)/2 + 1, [numInd 1]); % (fromSz-1)/2 + 1;
   centreTo   = repmat((toSz  -1)/2 + 1, [numInd 1]);
   % Round indices down if coming from larger slice, else round up (i.e.
   % need to choose whether 2.5 etc need to go up or down for even dist of
   % neuron locations across the laminar)
   ind_scaled = zeros(size(ind));
   if toSz(1)>fromSz(1)
      ind_scaled(ind> centreFrom) = ind(ind> centreFrom)*1.001;
      ind_scaled(ind<=centreFrom) = ind(ind<=centreFrom)*0.999;
   else
      ind_scaled(ind> centreFrom) = ind(ind> centreFrom)*0.999;
      ind_scaled(ind<=centreFrom) = ind(ind<=centreFrom)*1.001;
   end

   ind_rel    = round((ind_scaled - centreFrom).*scale + centreTo);
   ind_rel(ind_rel==0) = 1;

   if any(ind_rel < 1)
      fprintf('\ngetScaledLocationVector: relative index became less than max size %d \n',...
          toSz(1));
      fprintf('\tindex = [%s]\t', num2str(ind)); fprintf('\trelative index = [%s]\n\n', num2str(ind_rel));
   elseif any(ind_rel > toSz(1))
      fprintf('\ngetScaledLocationVector: relative index became greater than max size %d \n',...
          toSz(1));
      fprintf('\tindex = [%s]\t', num2str(ind)); fprintf('\trelative index = [%s]\n\n', num2str(ind_rel));
      error();
   end

end