% wdot = updateVisualCortexWeights(k1,k2,N,weights,Q,preind)
% Calculate change in weights as a function of input covariance, using
% Linsker's learning equation. This is for all synapses to a single
% postsynaptic neuron
% Inputs:
%  k1      - homeostatic constant
%  k2      - first order weight change
%  N       - expected number of synapses going to this postsynaptic neuron
%  weights - strengths of all synapses going to this post syn neuron
%  Q       - covariance of inputs (estimated or analytic)
%  conns   - list of all connections to postsynaptic neuron
%  n       - proportion of excitatory neurons (determines weight bounds)
%  dt      - length of timestep
% Outputs:
%  wdot    - change in weight (scalar)
function w = updateVisualCortexWeights(k1, k2, N, weights, Q, conns, dt)
   % wdot_ij = k1 + 1/N*sum((Q_ij + k2)*c_ij
   w = zeros(size(conns));
   for j=1:length(conns)
      % Q only includes covariance of connection neurons
      wdot(j) = k1 + 1/N * sum(weights * k2) ...
                   + 1/N * Q(j,:)*weights(:);
      w(j)    = weights(j) + wdot(j)*dt;
   end          
   % Update weights at the end of all layer inputs being calculated, as you
   % may remove the impact of a layer with negative covariance
%    w(w> n   ) = n;
%    w(w<(n-1)) = n-1;
end

