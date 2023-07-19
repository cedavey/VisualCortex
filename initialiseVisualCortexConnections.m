% weights = initialiseLinskerConnections(conns, wgtdist, wgtparams)
% Given connection indices and a distribution, initialise weight value for
% each connection using the given distribution.
% Inputs:
%  conns     - vector of univariate connection indices
%  wgtdist   - name of initialising distribution 
%  wgtparams - parameter values for initialising distribution
function weights = initialiseVisualCortexConnections(conns, wgtdist, wgtparams)
   N = size(conns,1);
   switch wgtdist
      case 'constant'
         c = wgtparams(1);
         weights = ones(N,1)*c;
      case 'gaussian'
         m = wgtparams(1);
         s = wgparams(2);
         % 2D gaussian distribution of weights
         % TO DO: perhaps make this periodic?? Does that make sense for weights?
         weights = rmvnrnd(m,s,N,[-eye(2); eye(2)], [0; 0; N; N]);
      case 'uniform'
         l = wgtparams(1);
         u = wgtparams(2);
         weights = rand(N,1)*(u-l) + l;
      otherwise
         error('initialiseVisualCortexConnections - unknown weight distribution (%s)',wgtdist);
   end
end