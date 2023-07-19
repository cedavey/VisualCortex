% [preIDs,postRel] = initVisualCortexGaussianConnections(N,radius,preSize,postSize,postID,Nsigma,radialconfig)
%
% This function determines which presynaptic neurons from the 
% are connected to each postsynaptic neuron. It returns an
% exact number of connections every time if no std dev in number of 
% connections is provided, else it returns a number of connections with an 
% expected value of N and stddev of Nsigma. Radius of connection density
% is determined by a Gaussian with provided variance.
% Inputs:
%  n        - number of connections to return exactly
%  radius   - variance of gaussian connection density
%  fromSize - size of presynaptic layer
%  toSize   - size of postsynaptic layer
%  postID   - index of current postsynaptic neuron
%  Nsigma   - std dev in number of connections
%  radialConfig - parameters associated with a radial configuration, in
%                 which connections are drawn from a non-circular
%                 distribution, with radius pontentially changing as a
%                 function of distance from the postsyn neuron
%
% Outputs:
%  preIDs   - ID of postsynaptic cell (i.e. index in matrix)
%  postRel  - relative position of postsyn cell when layers are stretched
%             and centred over each other
function [preIDs,postRel] = initVisualCortexGaussianConnections( N, radius, preSize, postSize, postID, Nsigma, radialconfig )
   if ~exist('Nsigma','var') || isempty(Nsigma)
      Nsigma = 0; 
   elseif isstruct( Nsigma ) && ~exist('radialconfig','var')
      % Nsigma was a late variable add on, so older code may have
      % radialconfig as an input in its place
      radialconfig = Nsigma;
   end
   if ~exist('radialconfig','var') || isempty( radialconfig )
      % if no radial configuration options provided, set all to false,
      % except periodic, which is set to true
      [radialDist,rayleighDist,radialSquRoot,uniqueLoc] = deal(false);
      periodic      = true;
      rayleighDist  = true;
   else
      periodic      = radialconfig.periodic;     % periodic boundaries
      radialDist    = radialconfig.radialDist;   % prob of conn changes with radial distance from postsyn neuron
      rayleighDist  = radialconfig.rayleighDist; % use rayleigh instead of normal distribution
      radialSquRoot = radialconfig.radialSquRoot;% radius is std dev instead of variance
      uniqueLoc     = radialconfig.uniqueLoc;    % no repeated conns btwn presyn and postsyn neurons
   end
   % allow scalar layer size inputs for square layers
   if isscalar( preSize ),  preSize = [ preSize  preSize]; end
   if isscalar(postSize),  postSize = [postSize postSize]; end
   debug = false;
   
   % In the fixed_pre connection type, the total number of connections is
   % fixed, but the exact neurons being connected is still random.
   % To be compatible with Linsker's layers, I'm going to connect from 
   % near the neuron's position, & emanate outwards using a Gaussian dist
    
   % Get size arrays in required format (no singleton dimensions etc), &
   % get location vector for post neuron relative to pre & post layer sizes
   % (e.g. if pre and post neuron layers are not the same size, centre the
   % postsyn layer above the presyn layer, and calculate postsyn neuron 
   % positions relative to presyn layer, so we can draw synapse connections
   % using our spatial Gaussian distribution
   [postInd,postRel] = getRelativePosition( preSize, postID, postSize );
   indRel = subv2ind(preSize,postRel); % get index of relative postsynaptic subscripts
   preIDs = [];
   if N > prod(preSize)
       fprintf('\n\n initialiseSpikeSim: connect error: n (%d) is too large for presyn group size\n',...
               N);
     return
   end
   radius = radius/sqrt(2); % conform to linsker's definition of radius
   
   % if variance of synaptic density is 0 conns are divergent instead of convergent, 
   % so single presyn neuron connects to many postsyn neurons, so just return 
   % relative location of postsynaptic neuron, in presynaptic layer coordinates
   if radius==0
      preIDs = repmat(postRel,[N 1]);
      return;
   end
   % if radius changes with distance from centre of laminar, calc new radius
   if radialDist
      radius = getRadialVar( radius, postSize, postInd, radialSquRoot );
   end
   
   % Put this in a while loop for the case that uniqueLoc is true, because
   % multiple connections to a postsyn neuron must be removed, and new
   % connections generated
   iteration = 0;
   
   % make exact number of connections random - assume gaussian
   N = randn(1)*Nsigma + N;
   N = ternaryOp( N<=0.51, 1, round(N) );
   while size(preIDs,1)<N && iteration<10 % cap time spent on this
     iteration = iteration + 1;

     % If we're generating with periodic boundaries then you can use an 
     % infinite Gaussian. If we're not using periodic boundaries then 
     % use a truncated multivariate gaussian.  
     if periodic
       % rician dist --> multivar gaussian coordinates
       if rayleighDist 
         % Generate random indices of presyn neurons, using a 
         % multivar Gaussian distribution around the current neuron. 
         if isscalar(radius)
            Sigma = eye(2)*radius^2;
         else
            Sigma = diag(radius.^2);
         end
         ind = mvnrnd([0 0], Sigma, N); % sigma is variance
         ind = bsxfun(@plus, ind, postRel); % shift to be centred at postsyn
         
       % gaussian in each dimension 
       else
         % Generate random indices of presyn neurons, using a Gaussian
         % distribution around the current neuron for the radius, and a
         % uniform distribution for the angle, then convert to cartesian

         % Randomly generate polar coordinates, with radius according
         % to Linsker's distribution & angle evenly distributed on [0, 2pi]
         r     = abs(normrnd(0,radius,N,1)); % sigma is std dev
         theta = rand(size(r))*2*pi;
         % 0 indexing until after mod, so subtract 1 from postRel
         ind   = [r.*cos(theta)  r.*sin(theta)] + repmat(postRel,[length(r) 1]);
       end
       
       % Apply periodic boundary conditions
       ind = periodicCellPosition(round(ind),preSize);

     else
        % truncated multivar Gauss bounded by hyperplanes AX<=B
        % x = rmvnrnd(mu, sig, N, A, B)
        % since we always have < to make AX > B you need to multiply by -1
        % so that the following gives -x <= -1 so that x >= 1 & x <= preSize
        ind  = rmvnrnd(postRel,eye(2)*radius, N, [-eye(2); eye(2)], [-1; -1; preSize(:)]);
        % Round indices to nearest integer, also converting to 1-indexing
        ind = round(ind); % generated indices on [0 fromSize-1], with postRel-1, so add 1 here
     end
     
      % If single synaptic connection allowed between neurons, remove duplicates
     if uniqueLoc
        ind = unique(ind,'rows'); % extract unique connections
     end
     
     % Add new positions to the presyn neuron matrix
     preIDs  = [preIDs; ind];
   end
   
   if size(preIDs,1)<N
      fprintf('Too many iterations to generate fixed Gaussian connections, exiting with %d instead of %d\n',...
              size(preIDs,1),N);
   end
   
	preIDs = preIDs(1:N,:);
   % since unique returns sorted array, choose randomly wh preIDs to keep
	preIDs = preIDs(randperm(size(preIDs,1),N),:);
   preIDs = fliplr(sortrows(fliplr(preIDs)));
end

