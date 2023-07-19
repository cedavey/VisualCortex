%% Use visual_cortex network to process input
% Network defined by:
%  - sz (layer sizes)
%  - layernames (layers included in the network)
%  - conns_loc (location of synapses connecting the network)
%  - weights (weights of synapses)
%  - PSP (psp structure)
%  - delays (delays of synapses)
%  - lambda (background layer rate)
%  - RF (receptive field structure)
%

% Optional 
%  - delayinds (delays in indices - specific to choice of dt)
%  - synloc (connection indices in subscripts)

% TO DO:
% - calculate phase as a function of temporal_freq
% - make number of timesteps as a function of number of phases
% - give input & calculate output
% - get maximal angular response
% - make script that gets tuning curve, evolves plasticity, & recalculates

d_angle         = 90/sz(1);  % assumes square layer
input_type     = 'grating';
grating_angles = 0:d_angle:180;
temporal_freq  = 2;
spatial_freq   = 2;

in = zeros(sz(1),sz(2),round(T/dt));
for ti=1:round((T/dt)),
   if any(strcmpi(intype,{'grate','gratings','grating'}))
      params = {ti*dt, grating_angle, temporal_freq, spatial_freq};
   end
   in(:,:,ti) = generateVisualCortexInput(lambda,sz,intype,params); 
   imagesc((in(:,:,ti))); 
   M(:,ti)=getframe; 
end; 



