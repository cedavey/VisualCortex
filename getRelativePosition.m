%   [locationVector, locationRelative] = getRelativePosition(sz,index,fromSz)
% Get a location vector for the
% appropriate dimension - e.g. [i j k] for 3D, or [i j] for 2D - if an
% index position is provided. If the size of the presynaptic layer is
% provided in 'fromSz', then the relative position of the postsynaptic
% neuron is also provided, with the 2 layers centred relative to one
% another.
%
% Inputs:
%  sz               - size of postsynaptic layer that the neuron resides in
%  index            - index of current neuron into the layer (scalar)
%  fromSz           - size of presynaptic layer if relative position required
% Outputs:
%  locationVector   - multi-dimensional index giving location 
%  locationRelative - gets the location in presynaptic coordinates assuming 
%                     that the 2 layers are centred around each other, & 
%                     stretched so that the end points meet
function varargout = getRelativePosition(preSz,index,postSz)
% Note: test over a range of inputs using:
%  sub = cell2mat(arrayfun(@(i) getRelativePosition(sz,i,fromSz),1:prod(sz),'uniformoutput',false)'); 
%  ind = subv2ind(sz,sub); 
%    a = accumarray(sub,ones(size(ind))); 
%  figure; imagesc(a); axis image; colorbar;
  % If the 'from' size is different to the 'to' size, then everything has
  % to be scaled, because the 'from' neuron that is directly 'above' the
  % 'to' neuron needs to be centred so that the 'from' layer is centred
  % directly above the 'to' layer. So calculate the scale of 'from' layer
  % to 'to' layer, & then shift so that middle points are aligned
  if nargin==3 && ~isempty(index)
    % get neuron's position as matrix subscripts
    if isscalar(index)
   	 locVector   = ind2subv(postSz,index); % neuron ID is an indices
    else
       locVector   = index; % index now input as multi-dim
    end

                % getScaledLocationVector(ind, fromSz, toSz)
    locRelative = getScaledLocationVector(locVector, postSz, preSz);
    if nargout==0 || nargout==1
       varargout{1}   = locRelative;
    else
       varargout{1}   = locVector;
       varargout{2}   = locRelative;
    end
    
  elseif nargin==2 && ~isempty(index)
    locVector      = cell(1,length(postSz));
    [locVector{:}] = ind2sub(postSz,index); % multi dimensional index
    locVector      = cell2mat(locVector);
	 varargout{1}   = locVector;		
  end
end
