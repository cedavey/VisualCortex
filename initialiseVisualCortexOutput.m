%
% [outrates,outweights,time,outQ,outsig,outmu] = ...
%          initialiseVisualCortexOutput(layerconfig,interval,intervalQ,rates,weights,co_inputs);
%
% or
% 
% [outrates,outweights,time,outQ,outsig,outmu] = ...
%          initialiseVisualCortexOutput(layerconfig,interval,intervalQ);
%
% From linsker configuration & output intervals, allocate memory for output
% rates & weights. Requires weights because number of connections to each
% neuron may vary depending on connection mechanism used (i.e. is mean
% number of connections equal to N, or does each postsyaptic neuron have
% exactly N connections to it?). If rate & weight structures aren't
% provided, it's assumed that each neuron has exactly N connections to it.
% Note: require rates to be non-empty if delay is non-trivial.
%
% If endpoints is set to true then allocate memory for first & last
% timestep regardless of whether they're included from interval or not.
% Qinterval is the interval between covariance matrix Q saving - these
% require a fuckload of memory so typically are saved less frequently.
%
% See also: visual_cortex,  processVisualCortexInput
function [outrates,outweights,time,outQ,outsig,outmu,co_outputs] = ...
            initialiseVisualCortexOutput(layerconfig,interval,intervalQ,varargin)
   tic;
   nargs = length(varargin);
   optargs = {[],[],[]};
   optargs(1:nargs) = varargin(:);
   [rates,weights,co_inputs] = optargs{:};
         
   [T, dt, estCov, layernames,layerconns, sz, N, endpoints, plastic] = ...
             struct2v(layerconfig,'T','dt','estCov','layernames','layerconns',...
                                  'sz', 'N' ,'endpoints','plastic');
   endpoints    = ternaryOp(isempty(endpoints),false,endpoints);
   estCov       = ternaryOp(isempty(estCov),   false,estCov);
   intervalQ    = ternaryOp(isempty(intervalQ),interval,intervalQ);
   interval     = ternaryOp(ischar(interval)  && strcmpi(interval, 'all'),1,interval);
   intervalQ    = ternaryOp(ischar(intervalQ) && strcmpi(intervalQ,'all'),1,intervalQ);
   [time,nT,nQ] = getLinskerTimeVec(T,dt,interval,endpoints,intervalQ);   
   outweights   = []; co_outputs = [];
   
   % Allocate memory for rate output 
   for li=1:length(layernames)
      layer   = layernames{li};
      if ~exist('rates','var') || isempty(rates)
         outrates.(layer) = zeros(max([nT 2]),prod(sz.(layer))); % min size
      else
         outrates.(layer) = repmat(rates.(layer)(end,:),[nT 1]);
      end
   end
   % Allocate memory for weights & covariance related terms
   for ci=1:size(layerconns,1)
      pre     = layerconns{ci,1};
      post    = layerconns{ci,2};
      cxlabel = [pre post];
      M       = prod(sz.(post)); % length(rates.(post)(end,:));
      
      if ~isempty(plastic) && plastic.(cxlabel)
         % If we're estimating covariance, allocate memory for cov & mean
         for i=1:M
            if isempty(weights)
               nC = N.(cxlabel); % assume all postsyn have exactly N presyn cxs
               outweights.(cxlabel){i} = zeros( nT, nC );
            else
               nC = length( weights.(cxlabel){i} );
               outweights.(cxlabel){i} = zeros(nT,nC);
               if endpoints
                  outweights.(cxlabel){i}(1,:) = weights.(cxlabel){i};
               end
            end
            if estCov
%                for j=1:M
                  outQ.(cxlabel){i}  = zeros(nQ,nC,nC);
                  outsig.(cxlabel){i}= zeros(nQ,nC);
                  outmu.(cxlabel){i} = zeros(nQ,nC);
%                end
            else
               outQ = []; outsig = []; outmu = [];
            end
         end
         if ~isempty(co_inputs)
            for cc=1:length(co_inputs.(cxlabel))
               co_layer      = co_inputs.(cxlabel){cc};
               co_out.Q      = zeros([nQ size(co_layer.Q)]);
               co_out.mu_in  = zeros([nQ size(co_layer.mu_in)]);
               co_out.mu_out = zeros([nQ size(co_layer.mu_out)]);
               co_out.sig_in = zeros([nQ size(co_layer.sig_in)]);
               co_out.sig_out= zeros([nQ size(co_layer.sig_out)]);
               co_outputs.(cxlabel){cc} = co_out;
            end
         end
      end
   end
   
   cprintf( 'Strings', printTime(toc,'\tInitialising visual cortex output structures tooks') );

end
