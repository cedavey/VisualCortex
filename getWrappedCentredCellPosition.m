% wrappos = getWrappedCentredCellPosition(sz, centre, pos)
% 
% Define a new layer centre & give cell position with respect to the new 
% centre, which can now be negative. So if layer size is [10 10], and we 
% have cell in position [10, 0], and we want the new centre to be [1, 1], 
% give the cell position as [-1, 0], since the distance of the cell to the
% new centre is 1, but in the wrapped direction, so negative
function wrappos = getWrappedCentredCellPosition(sz, origcentre, origpos)
   pos       = origpos - 1; % so start indices at 0 instead of 1
%    sz        = sz  - 1; % since indices go from 1 to max
   centre    = origcentre - 1; % start indices from 0
   pos       = pos - centre; 
   wrappos   = pos;
   
   % for each x, y coord, get actual distance & wrapped distance to centre
   [periodicDist, i] = min(abs([pos + sz ; pos - sz]));
   actualDist   = abs(pos);
   % set x & y separately coz easier
   wrappedDist = [pos + sz ; pos - sz];
   if actualDist(1) > periodicDist(1)
      wrappos(1) = wrappedDist(i(1),1); % get out periodic ind that minimised dist
   end
   if actualDist(2) > periodicDist(2)
      wrappos(2) = wrappedDist(i(2),2); 
   end

   wrappos = wrappos + 1; % so indices start from 1 again
end
