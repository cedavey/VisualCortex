% wdot = updateJointInputWeights(k1,k2,N,weights,Q,preind)
% Calculate changes in weights induced in the plastic pre -> post connections
% by a second presyn input to the postsyn layer, where the 2 presynaptic
% layers are correlated in some way.
% This function calculates all changes in presyn connections to a single
% postsynaptic neuron in the postsynaptic layer
% Inputs:
%  k1        - homeostatic constant
%  k2        - first order weight change, calculated btwn the pre layers
%  N         - expected number of synapses going to this postsynaptic neuron
%  weights   - strengths of all synapses from presyn layer to this postsyn neuron
%  weights_other - strength of synapses from co-presyn layer to this post
%  Q         - covariance of inputs (estimated or analytic)
%  conns     - list of all plastic connections to postsynaptic neuron
%  conns_co  - list of all connections between pre & other presyn layer, both
%              of which connect to the postsynaptic layer
%  conn_other- connections between the other presynaptic layer & postsyn layer
%  n         - proportion of excitatory neurons (determines weight bounds)
%  dt        - length of timestep
% Outputs:
%  wdot      - change in weight (scalar)
function w = updateJointInputWeights(k2,N,weights,weights_other,Q_co,conns,conns_co,conns_other,n,dt)
%    weights.(cxlabel){i} =                             ...
%          updateJointInputWeights(colayer.k2,          ...
%                              N.(other_input),         ...
%                              weights.(other_input){i},...
%                              Q.(colabel){i},          ...
%                              conns.(cxlabel){i},      ...
%                              conns.(other_input){i},  ...
%                              n.(cxlabel), dt);

   % conns - connection indices from other input layer -> postsyn layer
   % extract these covariance values from presyn layer -> joint input layer
   w    = zeros(size(conns)); 
   % initialise weight change to change induced by k2 term
   wdot = zeros(size(conns)); % 
%    wdot = ones(size(conns))*(1/N * sum(weights_other * k2)); % zeros(size(conns)); 
   for j=1:length(conns) % pre -> post connections
      % Get cov btwn pre & co_pre, scaled by weight btwn co_pre & post
      for l=1:length(conns_other)
         % See if neuron in the co-presynaptic layer has connection to
         % presynaptic neuron (i.e. are they correlated)
%          if any(j==conns_co{l})
            % presyn layer neurons connected --> get cov btwn presyn
            % layers, & wgt from co-presyn layer to postsyn layer
%             jl_ind = find(j==conns_co{l},1,'first'); % if multiple cov they'll be identical
%             q_jl   = Q_co(jl_ind,l);
            q_jl   = Q_co(conns(j),l);
            
            % is this co-presyn neuron connected to postsyn neuron?
%             if any(conns_other==l)
%                il_ind  = conns_other==l; % if multiple add impact of all
%                w_il    = weights_other(il_ind);
%                wdot(j) = wdot(j) + 1/N * sum(w_il*q_jl);
               
               wdot(j) = wdot(j) + 1/N * weights_other(l)*(q_jl + k2);
% w_(end+1) = 1/N * weights_other(l)*(q_jl + k2);
%             end
%          end
      end
      w(j) = weights(j) + wdot(j)*dt;
%       Q_ = Q(j,conns_other);
%       wdot(j) = 1/N * sum(weights * k2) ...
%                 + 1/N * Q*weights(:);
   end
%    w(w> n   ) = n;
%    w(w<(n-1)) = n-1;
end

