%       [rad, theta] = periodicCellDistance(pos1,pos2,sz,scalar,periodic)
%
% Gets norm2 distance between 2 cells assuming periodic boundaries
% Inputs:
%  pos1   - single cell position or matrix of 2D positions
%  pos2   - single cell position or matrix of 2D positions or, if empty, a
%           copy of pos1 is made & pairwise distance is calculated between
%           pos1 cells
%  sz     - size of cell layer to wrap around
%  scalar - return result as norm 2 scalar distance (default) or, if false,
%           return results in vector format (each distance a row vector)
%           TO DO: vector option not yet implemented
%  periodic - true for periodic boundaries (defaults to true)
%
% Outputs:
%  dist   - minimum distance btwn 2 positions, accounting for periodicity
%           if requested
%  angle  - angle between 2 positions
%           
% See also:          periodicCellPosition       periodicDistance
function [ dist, angle ] = periodicCellDistance( pos1, pos2, sz, scalar, periodic )
   if nargin==0
      help periodicCellDistance
      return;
   end
   dist=[]; angle=[];
   % scalar = false --> return 2D Cartesian distances (i.e. not abs val)
   if ~exist('scalar','var')   || isempty(scalar),     scalar = true; end
   if ~exist('periodic','var') || isempty(periodic), periodic = true; end
   if ~exist('sz','var') || isempty(sz)
      cprintf('Error','\nperiodicCellDistance input error - need layer size\n\n');
      return;
   end
   if isscalar(sz),    sz = [sz sz]; end
   if isempty(pos2), pos2 = pos1;    end
   if isempty(pos1), pos1 = pos2;    end
   % if location is 1D (i.e. scalar) convert to 2D with 0 position in 2nd dim
   if size(pos1,2)==1, pos1 = [pos1 ones(size(pos1))]; end
   if size(pos2,2)==1, pos2 = [pos2 ones(size(pos2))]; end
   
   if periodic
      distfn = @(d,l) min(abs([mod(d,l) l-mod(d,l)]),[],2);
   else
      distfn = @(d,l) d;
   end
   
   % If pos1 or pos2 has only a single 2D position
   if size(pos1,1)==1 || size(pos2,1)==1
      if scalar
         dist = sqrt( distfn(pos1(:,1)-pos2(:,1), sz(1)).^2 + distfn(pos1(:,2)-pos2(:,2), sz(2)).^2 );
      else
         dist = [distfn(pos1(:,1)-pos2(:,1),sz(1)),  distfn(pos1(:,2)-pos2(:,2),sz(2))];
      end
      angle = atan2(distfn(pos1(:,2)-pos2(:,2),sz(1)), distfn(pos1(:,1)-pos2(:,1),sz(2)));
      
   % If both pos1 && pos2 have multiple 2D positions
   else
      if scalar
         dist  = cell2mat(arrayfun(@(i) sqrt(distfn(pos1(:,1)-pos2(i,1),sz(1)).^2 ...
                                           + distfn(pos1(:,2)-pos2(i,2),sz(2)).^2),...
                          1:length(pos2), 'uniformoutput',0));
         angle = cell2mat(arrayfun(@(i) atan2(distfn(pos1(:,2)-pos2(i,2),sz(1)), ...
                                              distfn(pos1(:,1)-pos2(i,1),sz(2))),...
                          1:length(pos2), 'uniformoutput',0));
      else
         dist  = cell2mat(arrayfun(@(i) [distfn(pos1(:,1)-pos2(i,1),sz(1)),  ...
                                         distfn(pos1(:,2)-pos2(i,2),sz(2))],...
                                         (1:length(pos2))', 'uniformoutput',0));
         angle = cell2mat(arrayfun(@(i) atan2(distfn(pos1(:,2)-pos2(i,2),sz(1)), ...
                                              distfn(pos1(:,1)-pos2(i,1),sz(2))),...
                                         (1:length(pos2))', 'uniformoutput',0));      
      end
   end

end



% function [dist,angle] = periodicDistance(pos1,pos2,sz,scalar)
%    % periodic distance with layer of size l in specific direction
%    distfn = @(d,l) min(abs([mod(d,l) l-mod(d,l)]),[],2); 
%    
%    % If pos1 or pos2 has only a single 2D position
%    if size(pos1,1)==1 || size(pos2,1)==1
%       if scalar
%          dist = sqrt(distfn(pos1(:,1)-pos2(:,1),sz(1)).^2 + distfn(pos1(:,2)-pos2(:,2),sz(2)).^2);
%       else
%          dist = [distfn(pos1(:,1)-pos2(:,1),sz(1)),  distfn(pos1(:,2)-pos2(:,2),sz(2))];
%       end
%       angle = atan2(distfn(pos1(:,2)-pos2(:,2),sz(1)), distfn(pos1(:,1)-pos2(:,1),sz(2)));
%       
%    % If both pos1 && pos2 have multiple 2D positions
%    else
%       if scalar
%          dist  = cell2mat(arrayfun(@(i) sqrt(distfn(pos1(:,1)-pos2(i,1),sz(1)).^2 ...
%                                            + distfn(pos1(:,2)-pos2(i,2),sz(2)).^2),...
%                           1:length(pos2), 'uniformoutput',0));
%          angle = cell2mat(arrayfun(@(i) atan2(distfn(pos1(:,2)-pos2(i,2),sz(1)), ...
%                                               distfn(pos1(:,1)-pos2(i,1),sz(2))),...
%                           1:length(pos2), 'uniformoutput',0));
%       else
%          dist  = cell2mat(arrayfun(@(i) [distfn(pos1(:,1)-pos2(i,1),sz(1)),  ...
%                                          distfn(pos1(:,2)-pos2(i,2),sz(2))],...
%                                          (1:length(pos2))', 'uniformoutput',0));
%          angle = cell2mat(arrayfun(@(i) atan2(distfn(pos1(:,2)-pos2(i,2),sz(1)), ...
%                                               distfn(pos1(:,1)-pos2(i,1),sz(2))),...
%                                          (1:length(pos2))', 'uniformoutput',0));      
%       end
%    end
% end
% 
% % Calculate the effective 2D distance btwn 2 neurons, NOT using
% % periodic boundaries. 
% function [dist,angle] = nonPeriodicDistance(pos1,pos2,sz,scalar)
%    distfn = @(d,l) d; 
%    
%    % If pos1 or pos2 has only a single 2D position
%    if size(pos1,1)==1 || size(pos2,1)==1
%       if scalar
%          dist = sqrt(distfn(pos1(:,1)-pos2(:,1),sz(1)).^2 + distfn(pos1(:,2)-pos2(:,2),sz(2)).^2);
%       else
%          dist = [distfn(pos1(:,1)-pos2(:,1),sz(1)),  distfn(pos1(:,2)-pos2(:,2),sz(2))];
%       end
%       angle = atan2(distfn(pos1(:,2)-pos2(:,2),sz(1)), distfn(pos1(:,1)-pos2(:,1),sz(2)));
%       
%    % If both pos1 && pos2 have multiple 2D positions
%    else
%       if scalar
%          dist  = cell2mat(arrayfun(@(i) sqrt(distfn(pos1(:,1)-pos2(i,1),sz(1)).^2 ...
%                                            + distfn(pos1(:,2)-pos2(i,2),sz(2)).^2),...
%                           1:length(pos2), 'uniformoutput',0));
%          angle = cell2mat(arrayfun(@(i) atan2(distfn(pos1(:,2)-pos2(i,2),sz(1)), ...
%                                               distfn(pos1(:,1)-pos2(i,1),sz(2))),...
%                           1:length(pos2), 'uniformoutput',0));
%       else
%          dist  = cell2mat(arrayfun(@(i) [distfn(pos1(:,1)-pos2(i,1),sz(1)),  ...
%                                          distfn(pos1(:,2)-pos2(i,2),sz(2))],...
%                                          (1:length(pos2))', 'uniformoutput',0));
%          angle = cell2mat(arrayfun(@(i) atan2(distfn(pos1(:,2)-pos2(i,2),sz(1)), ...
%                                               distfn(pos1(:,1)-pos2(i,1),sz(2))),...
%                                          (1:length(pos2))', 'uniformoutput',0));      
%       end
%    end
% 
% end














