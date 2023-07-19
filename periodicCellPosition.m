% pos = periodicCellPosition(pos,sz)
% Get cell position for periodic boundaries
function pos = periodicCellPosition(pos,sz)
   modDim   = @(x,szx) mod(x,szx);
   pos(:,1) = modDim(pos(:,1), sz(1));
   pos(pos(:,1)==0,1) = sz(1);
   pos(:,2) = modDim(pos(:,2), sz(2));
   pos(pos(:,2)==0,2) = sz(2);
end
