% [ntwkconfig, layerconfig] = prepareVisualCortexNetwork(layerconfig)
function [ntwkconfig, layerconfig] = prepareVisualCortexNetwork(layerconfig)
   ntwkconfig = [];
   [layernames, layerconns, inlayer,         ...
    retina, retina_on, retina_off,           ...
    RF, sz, N, lambda, plastic, wgtparams,   ...
    noisedist, noiseparams, r,               ...
    kb, k1, k2, Rb,estCov,continueSim] = struct2v(layerconfig,    ...
                       'layernames', 'layerconns',  'inlayer',    ...
                       'retina',     'retina_on',   'retina_off', ...
                       'RF',         'sz',          'N',          ...
                       'lambda',     'plastic',     'wgtparams',  ...
                       'noisedist',  'noiseparams', 'r',          ...
                       'kb',         'k1',          'k2', 'Rb',   ...
                       'estCov',     'continueSim');

   haveRetina = ternaryOp(isempty(RF), false, true); 
	if ~estCov
      initTime = 0;              % nothing to initialise if not estimating cov
      layerconfig.initTime = 0;  % update for when saved to file
   end
   % if inputting natural images then get names of images to use now
   % use user input param of number of timesteps to display each img for,
   % with total number of timesteps in simulation, to determine num imgs
   if strcmpi(noisedist, 'natural')
      noiseparams = naturalImageParams( layerconfig, noiseparams );
      if isempty( noiseparams )
         return;
      end
      layerconfig.noiseparams = noiseparams;    % udpate layerconfig struct with changes
   end
   
   %% Layer specific rate update parameters   
%    f02.LC      = Ra.C - lambda.B*N.LC*Rb.LC*k1.LC/k2.LC;  
   eta_fn = @(N,kb,Rb,f02) N*kb*Rb*f02;
   % Init Ra to 0 for each layer, & then update based on cx info
   for li=1:length(layernames)
      Ra.(layernames{li}) = 0;
   end
   if haveRetina
      % first layer (retina) needs to have Ra set by the mean value of 
      % the RF kernel multiplied by input rate
      muKernel   = sum(RF.retina_kernel(:)); 
      if ~isempty(retina)
         Ra.(retina)     = lambda.(retina)     - muKernel*lambda.(inlayer); 
      end
      if ~isempty(retina_on)
         Ra.(retina_on)  = lambda.(retina_on)  - muKernel*lambda.(inlayer); 
      end
      if ~isempty(retina_off)
         Ra.(retina_off) = lambda.(retina_off) + muKernel*lambda.(inlayer); 
      end
   end
   for ci=1:length(layerconns)
      pre     = layerconns{ci,1};
      post    = layerconns{ci,2};
      cxlabel = [pre post];
      Ra_curr = 0;
%          if exist('Ra','var') && isfield(Ra,post)
%             Ra_curr = Ra.(post);
%          else
%             Ra_curr = 0;
%          end
%          if Ra.(post)==0 % only allow Ra to be updated once if already >0
      if ~plastic.(cxlabel)
         % If connections aren't plastic then mean input rate will depend
         % on number of presyn neurons & weight initialisation
         w = mean(wgtparams.(cxlabel)); % works for constant & uniform, but not Gaussian
%             Ra.(post)  = lambda.(post) - N.(cxlabel)*Rb.(cxlabel)*lambda.(pre)*w; % Background noise
         Ra_curr    = lambda.(post) - N.(cxlabel)*Rb.(cxlabel)*lambda.(pre)*w; % Background noise
      else
         % For plastic connections the mean input rate will eventually
         % depend on weight fixed point, which is -k1/k2 for small covariance
         if k2.(cxlabel)==0
%                Ra.(post)  = lambda.(post);
            Ra_curr    = lambda.(post);
         else
%                Ra.(post)  = lambda.(post) + N.(cxlabel)*Rb.(cxlabel)*k1.(cxlabel)/k2.(cxlabel)*lambda.(pre);
            Ra_curr    = lambda.(post) + N.(cxlabel)*Rb.(cxlabel)*k1.(cxlabel)/k2.(cxlabel)*lambda.(pre);
         end
%             Ra.(post)     = ternaryOp(isinf(Ra.(post)),lambda.(post), Ra.(post));
         Ra_curr       = ternaryOp(isinf(Ra_curr), lambda.(post), Ra_curr);
         f02.(cxlabel) = lambda.(pre)*N.(cxlabel)*Rb.(cxlabel); % layer rate
         eta.(cxlabel) = eta_fn(N.(cxlabel),kb.(cxlabel),Rb.(cxlabel),f02.(cxlabel));  % layer learning rate
      end
      Ra.(post)   = Ra.(post) + Ra_curr;
%          end
   end
   layerconfig.Ra = Ra; layerconfig.f02 = f02; layerconfig.eta = eta;

   %% Initialise
	if haveRetina
      if any(strcmpi(RF.biasdist, {'levick', 'uniform'}))
         % Levick's paper
         Lbins = (0.05:0.05:0.5) - 0.025; % the bin positions are edges not centres
         Lbars = [23 54 55 59 28 20 6 3 1 1];
         Lbars = Lbars / sum(Lbars); 

         maxbias = RF.sd_cent_ver * 3;
         minbias = RF.sd_cent_ver + 1;
         Nbias   = length(Lbins); 
         dbias   = (maxbias - minbias) / (Nbias-1);
         biasvec = minbias : dbias : maxbias;
         rf      = zeros( RF.kernel_size, RF.kernel_size, Nbias );
         for i=1:Nbias
            RF.sd_cent_hor = biasvec(i);
            rf(:,:,i)      = retina_kernel(RF,false);
         end

         % Generate samples of bias from requested distribution
         if haveRetina
            if ~isempty(retina)
               NR = prod(layerconfig.sz.(layerconfig.retina)); % number of retinal cells
            elseif ~isempty(retina_on)
               NR = prod(layerconfig.sz.(layerconfig.retina_on)); % number of retinal cells
            elseif ~isempty(retina_off)
               NR = prod(layerconfig.sz.(layerconfig.retina_off)); % number of retinal cells
            end
            switch RF.biasdist
               case 'levick'
                  Lpdf  = Lbars; % convert to probability 
                  Lcdf  = cumsum(Lpdf);
                  u     = rand(NR, 1); % uniform random number for each cell
                  a     = arrayfun(@(r) find(r<=Lcdf,1,'first'), u); % compare random num to cum prob
                  bias  = Lbins(a);
                  % if we assume that the spread of levick's bias between 0 & 0.5
                  % matches the spread of our bias from minbias to maxbias, then we can
                  % just grab our kernels in rf at the desired frequency,
                  % so make 'a' a list of indices of specified distribution
                  bias  = a;
               case 'uniform'
                  bias  = ceil(rand(NR, 1)*10);
            end

            % Generate RF based on bias (all at 0 degrees - rotate next)
            RF.retina_kernel = cell(NR,1); % remove existing kernel & replace it with cell of kernels
            for ri=1:NR
               RF.retina_kernel{ri,1} = rf(:,:,bias(ri));
            end
            RF.bias = bias; % record distribution
         end
      end
      
      % Now rotate kernels 
      if ~isempty(retina)
         [ RF.kernels, RF.theta ] = initVisualCortexOrientation(RF.retina_kernel, sz.(retina), RF.angleSet, RF.angleProb);
      elseif ~isempty(retina_on)
         [ RF.kernels, RF.theta ] = initVisualCortexOrientation(RF.retina_kernel, sz.(retina_on), RF.angleSet, RF.angleProb);
      elseif ~isempty(retina_off)
         [ RF.kernels, RF.theta ] = initVisualCortexOrientation(RF.retina_kernel, sz.(retina_off), RF.angleSet, RF.angleProb);
      end
      layerconfig.RF = RF;
      % if DEBUG, plot a few example kernels to get a sense of
      % extent of bias & its distribution
      % plot example kernels of particular bias & orientation
      %          figure; 
      %          for i=1:min([NR 30])
      %             subplot(5,6,i), 
      %                imagesc(RF.kernels{i}); 
      %                axis image; 
      %                colorbar; 
      %          end

	end

   ntwkconfig = initialiseVisualCortexNetwork(layerconfig, inlayer, retina, retina_on, retina_off);


end