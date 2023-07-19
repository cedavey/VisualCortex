%  rates = applyEPSP(inrates,psp,dt)
% Apply an EPSP filter to rates generated in this timestep for a vector of
% neurons, to generate the actual current rates of the neurons
% Inputs:
%  inrates  - matrix of input rates, inrates = timestep x neuron
%  psp      - postsynaptic potential kernel, assumed to be generated using
%             the same timestep at which rates are calculated
% recursive - true -> apply to all samples except beginning ones that are
%             used to initialise (i.e. first length(psp) samples)
%             false -> calculate only end sample for all neurons
function rates = applyEPSP(inrates,psp,dt,recursive)
   % Note: it's more than 10 times faster doing my way then it is using
   % conv when you have a matrix rather than a vector
   if isempty(psp)
      rates = inrates;
   end
   nt = length(psp);
   if nargin<4 || isempty(recursive) || ~recursive
      nt    = min([size(inrates,1) nt]); % cull included timesteps if at beginning of sim
      rates = psp(1:nt)'*inrates(end:-1:(end-nt+1),:)*dt;
   else
%       % do for all timepoints, even if we don't have enough samples
%       for ti=1:size(inrates,1)
%          num = min([ti nt]);
%          rates(ti,:) = psp(1:num)'*inrates(ti:-1:(ti-num+1),:)*dt;
%       end
      % start only after we have enough samples for whole psp
      for ti=nt:size(inrates,1)
         rates(ti,:) = psp'*inrates(ti:-1:(ti-nt+1),:)*dt;
      end
   end
end