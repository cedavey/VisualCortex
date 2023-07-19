%    img = genCentredWeightImage(sz_pos, sz_pre, postID, weight_img)
% OR
%    img = genCentredWeightImage(sz_pos,sz_pre,postID,weights,conninds)
% Convert weights image to centred weight image. Can input an image, or a
% list of weights & connection indices
% Inputs:
%  sz_pre   - size of presynaptic layer
%  postID   - index of postsynaptic neuron
%  weights  - a 2D weight image OR a vector of weights
%  conninds - if a vector of weights provided, connection indices into
%             presynaptic layer must be provided to enable img constrution
%  defVal   - default value of pixels with no connection (defaults to NaN)
% Outputs:
%  img      - 2D centred weight image
%  col      - colour coded distance btwn presyn locations & postsyn neuron
function [img,col] = genCentredWeightImage(sz_pos, sz_pre, postID, weights, conninds, defVal)
   if nargin==0, help genCentredWeightImage; return; end
   
   if ~exist('defVal','var') || isempty(defVal), defVal = NaN; end
   img      = zeros(sz_pre)*defVal;  % make a 2D image of final weights
   cent_pre = fix(sz_pre/2);         % centre of presyn layer
   cent_pos = fix(sz_pos/2);         % centre of postsyn layer
   [~,loc,loc_rel]  = getGroupsLocationArray(sz_pos,postID,sz_pre);
   
   if nargin<5 || isempty(conninds)
      if isDim(weights,2)
         img = flipud(weights);
         conninds = 1:prod(sz_pre);
      else
         error('genCentredWeightImage: if conn inds not provided, weights must be a 2D img');
      end
   else
      if isempty(weights)
         weights = ones(size(conninds));
      end
   end
   
   dist_pos = loc_rel - cent_pre;

   % Get connection location, centre, mod for periodic boundaries, &
   % convert to indices
   nC = length(conninds);
   conn_loc = ind2subv(sz_pre,conninds);
   conn_loc = conn_loc - repmat(dist_pos,[nC 1]);
   conn_loc = mod(conn_loc,repmat(sz_pre-1,[size(conn_loc,1) 1]))+1;
   conninds = subv2ind(sz_pre,conn_loc);

   img(conninds) = weights; % populate image
   
   % Plot timeseries of weights, & shade the plot lines to encode distance 
   % of presynaptic neuron to postsynaptic neuron
%    d        = @(j) sqrt(sum(min(abs([cent_pre-j; sz_pre+1-j])).^2,2));
%    dist_pre = sqrt(sum((conninds - repmat(cent_pre,[length(conninds) 1])).^2,2));
   d        = @(j) sqrt(sum(j-repmat(cent_pre,[nC 1]),2).^2); % centred so no periodic boundaries
   dist_pre = d(conn_loc);
      ming  = 0.05;
      maxg  = 0.95;
   if length(unique(dist_pre))==1
      scale = unique(dist_pre)*maxg;
   else
      scale = (maxg-ming)/(max(dist_pre)-min(dist_pre)); % scale onto [0.2 0.9] according to distance
   end
   col = repmat((dist_pre-min(dist_pre))*scale,[1 3]) + ming; % RGB;
   col = mat2cell(col,ones(size(col,1),1),3);
end





