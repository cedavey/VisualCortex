% outrates = updateVisualCortexRates(Ra,Rb,curr_rates,inrates,weights)
% Calculate postsynaptic layer output rate for a single postsynaptic neuron
%  as a function of input rates & synaptic connection strengths.
% Inputs:
%  Ra       - background noise rate
%  Rb       - scale of input
%  rates    - current rates of all neurons in the input layer to allow
%             summing of multiple input groups feeding in to layer
%  weights  - connection strength from each input to this neuron
%  conns    - indices of presynaptic neurons
%  delays   - axonal delays (in number of timesteps) for each connection
%             defaults to 0 for all connections if not provided. If it is
%             provided, there must be an entry for every connection
%  psp      - to convolve inputs with a postsynaptic potential kernel,
%             provide the kernel (generated with same timestep as rates)
%  dt       - for epsp you can also input dt (defaults to 1/sum(psp), which
%             assumes that your psp kernel intergrates to 1
% Outputs:
%  outrates - calculated output rates as F_j = Ra + Rb*sum(weight*F_i(t-delay))
function outrates = updateVisualCortexRates(Rb,inrates,weights,curr_rates,rectify)
   if ~exist('rectify','var'), rectify = false; end
   % In case of multiple layers feeding into neuron, sum input & current
   % value of neuron's rate (initialised to 0 at start of timestep)   
   if exist('curr_rates','var')
      outrates = curr_rates + Rb*(weights'*toVec(inrates));
   else
      outrates = Rb*(weights'*toVec(inrates));
   end
   if rectify, outrates(outrates<0) = 0; end
end






