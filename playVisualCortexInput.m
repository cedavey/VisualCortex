% movie = playVisualCortexInput(sz,T,dt,intype,params,lambda)
% Play visual cortex input over a given time. 
% Inputs:
%  sz     - size of input (2D)
%  T      - length of time to play for
%  dt     - time resolution
%  intype - time of visual cortex input to generate
%     gaussian 
%     uniform
%     gratings
%  params - params of input type
%     gaussian: requires std dev
%     uniform:  requires single bound for [-bound bound]
%     gratings: requires angle of grating, temporal frequency, and 
%               spatial frequency
% E.g. 
%  playVisualCortexInput([10 10],1,0.1,'gratings',{70 1, 10},1);
% 
% See also:    generateVisualCortexInput
function M = playVisualCortexInput(varargin)
   nargs = length(varargin);
   optargs = {[10 10],10, 0.1,'gratings',[45 10], 10};
   optargs(1:nargs) = varargin(:);
   [sz,T,dt,intype,params,lambda] = optargs{:};
   
   % Grating input requires addition time input
   if any(strcmpi(intype,{'grate','gratings','grating'}))
      grating_angle = params{1};
      temporal_freq = params{2};
      spatial_freq  = params{3};
   end
   figure; colorbar;
   in = zeros(sz(1),sz(2),round(T/dt));
   for ti=1:round((T/dt)),
      if any(strcmpi(intype,{'grate','gratings','grating'}))
         params = {ti*dt, grating_angle, temporal_freq, spatial_freq};
      end
      in(:,:,ti) = generateVisualCortexInput(lambda,sz,intype,params); 
      imagesc((in(:,:,ti))); 
      M(:,ti)=getframe; 
   end; 
end