%    combos = allCombos( vlengths )
% OR
%    [y, x,...] = allCombos(nX, nY)
% Create a matrix of all combinations of input vectors. For each vector,
% just the length is supplied, and the vector is assumed to be 1:length.
% Inputs:
%   vlengths - vector containing length of each dimension's vector
% Outputs:
%   combos - matrix of all combinations
%
% See also: combvec    fullfact      combnk
%     e.g. to list all pairs with no repeated elements
%          x = fliplr(fullfact([N N]))
%          x(~diff(x')',:)=[]
%     e.g. the same but don't care about the order of elmts within each pair
%          x = combnk(1:N,2);
function varargout = allCombos(vlengths,varargin)
   % I keep inputting sizes in different input parameters - merge
   if ~isempty(varargin)
      vlengths = [toVec(vlengths) toVec(varargin{:})];
   end
    vlengths = vlengths(:)';
    nv = length(vlengths);  % number of vectors
    % each vector has the number of times it's pattern repeats, & the
    % number of times each number in the pattern is repeated
    vr = [1 cumprod(vlengths)];
    nc = vr(end); % num combos is last vector index 
    vr = vr(1:end-1); % don't need total combos agains
    nr = [1 cumprod(vlengths(end:-1:1))]; 
    nr=nr(end-1:-1:1); % reverse & ditch last value 
    
    combos = zeros(nc,nv); % final vector of combinations
    
    for i=1:nv
        v = 1:vlengths(i);
        n = nr(i);
        r = vr(i);
        p = repmat(v(:),[1 n])'; % pattern to repeat
        c = repmat(p(:),[r 1]);
        
        combos(:,i)=c;
	 end
	 % to make it compatible with sub2ind functions, flip it left -> right
	 combos = fliplr(combos);
    varargout = cell(1,nargout);
    if nargout==1
       varargout{1} = combos;
    elseif nargout>1
       for ii=1:size(combos,2)
          varargout{ii} = combos(:,ii);
       end
    end
end
