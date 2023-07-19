% Make a layer of retinal kernels with different bias & calculate the
% tuning curves for it
distype     = 'levick'; % 'levick', 'uniform'
sfreqs      = [0.2:0.1:0.5];   % spatial frequency/ies to get tuning curves for
genkernels  = true ; % generate retinal kernels according to a bias distribution
plotkernels = true ; % plot kernels for full range of bias parameter (once per param value)
plotbias    = false; % plot distribution of bias params we're giving to kernels
plotresults = true ; % plot estimated bias of retinal cells
NR          = 10;    % num retina cells in each dimension - total number of cells is NR^2

if genkernels
   % Levick's paper
   Lbins = (0.05:0.05:0.5) - 0.025; % the bin positions are edges not centres
   Lbars = [23 54 55 59 28 20 6 3 1 1];
   Lbars = Lbars / sum(Lbars); 
   Lmean = sum( Lbars .* Lbins ); 
   Lvar  = sum( (Lbins - Lmean).^2 .* Lbars );
   % from isabelle's tracy mc-traceface
   % 23.03554365	54.20632318	55.23393126	59.25873011	28.1735872	20.20962143	6.251280407	3.168456162	1.284510102

   ksz = 10;
   RF.kernel_size     = ksz;
   RF.sd_surr_hor     = 7; 
   RF.sd_surr_ver     = 7;
   RF.sd_cent_ver     = 5; 
   RF.coeff_cent      = 6.00e-4;
   RF.coeff_surr      = 1.30e-4;
   RF.ratio_surr_cent = 2;
   RF.normfn          = 'power';
   RF.input_var       = 5;
   [RF.angleSet, RF.angleProb] = catRetinaOS(10);
   RF.angleSet        = num2cell(RF.angleSet);         % angle options, in cell to allow 'radial'

   maxbias = 14;
   minbias = RF.sd_cent_ver + 1;
   Nbias   = length(Lbins); 
   dbias   = (maxbias - minbias) / (Nbias-1);
   biasvec = minbias : dbias : maxbias;
   rf   = zeros(ksz,ksz,Nbias);
   for i=1:Nbias
      RF.sd_cent_hor = biasvec(i);
      rf(:,:,i)      = retina_kernel(RF,false);
   end

   if plotkernels
      figure; nr = 2; nc = 5;
      cmap = redgreencmap(128);
      rf   = zeros(ksz,ksz,Nbias);
      for i=1:Nbias
         subplot(nr,nc,i),
         RF.sd_cent_hor = 5+i-1;
         rf(:,:,i)      = retina_kernel(RF,false); % receptive field kernel
         imagesc(rf(:,:,i), [-0.4 0.4]); axis image; 
         colormap(cmap(end:-1:1,:));
         colorbar;
      end
      % Plot cut through the middle
      % figure; for i=1:10, subplot(nr,nc,i), plot(rf(:,:,i)); end
      % figure; for i=1:10, subplot(nr,nc,i), plot(rf(:,:,i)'); end
   end

   % Allocate distribution of bias uniformly, or like Levick or ???
   switch distype
      case 'levick'
         Lpdf  = Lbars; % convert to probability 
         Lcdf  = cumsum(Lpdf);
         u     = rand(NR^2, 1); % uniform random number for each cell
         a     = arrayfun(@(r) find(r<=Lcdf,1,'first'), u); % compare random num to cum prob
         bias  = Lbins(a);
         % if we assume that the spread of levick's bias between 0 & 0.5
         % matches the spread of our bias from minbias to maxbias, then we can
         % just grab our kernels in rf at the desired frequency
         bias  = a;
      case 'uniform'
         bias  = ceil(rand(NR^2, 1)*10);
   end

   for ri=1:NR^2
      RF.retina_kernels{ri,1} = rf(:,:,bias(ri));
   end
   RF.bias = bias;
   [RF.kernels,RF.theta] = initVisualCortexOrientation(RF.retina_kernels,[NR NR],RF.angleSet,RF.angleProb);
   
   if plotbias
      figure; 
         % plot distribution we've created of sd_cent_hor bias parameter
         [y, x] = hist(bias);
         subplot(1,2,1), bar(x, y/sum(y));
         subplot(1,2,2), bar(Lbins, Lbars);
      
      % plot example kernels at specific bias & orientation
      figure; 
      for i=1:min([NR^2 30])
         subplot(5,6,i), 
            imagesc(RF.kernels{i}); 
            axis image; 
            colorbar; 
      end
   end

   %% Configure
   runSim             = true ;                % run simulation or just plot figures
   evolveSep          = false;                % evolve plastic layers separately
   estCov             = true;                 % estimate cov or use analytical result
   initCovWithAnal    = true;                 % initialise empirical cov with analytical cov
   calcCorr           = false;                % plot covariance matrix at the end
   endpoints          = true;                 % record start & end rates/weights
   layerdist          = 0;                    % distance between layers (to calc. delay)

   initTime           =  0;                   % num timesteps used to stabilise cov est
   interval           = 10;                   % num samples of rates/weights to output
   intervalQ          = 10;                   % num samples of cov/mu to include
   normQ              = true;                 % normalise covariance after estimation

   dt                 = 1e-1;                 % time resolution
   T                  = 2e-1;                 % num seconds in the simulation


   % EPSP configuration
   psptype            = 'delta';              % postsynaptic response after spike
   pspparams          = 0.15;                 % params for PSP
   %  psptype            = 'double';      
   %  pspparams   = [0.225 0.25]; % rise time, decay time
   pspdelay           = 1;                    % length of psp in time samples
   psp                = [];

   %% LAYER PARAMETERS
   layernames         = {'R'};                % neuron layers 
   layerconns         = {};                   % layer cxs
   inlayer            = 'I'; % name of input layer for readability
   retina             = 'R'; % name of retina layer for readability

   sz.I               = [100 100];                              % input image
   sz.R               = [ NR  NR];                              % retina
   spatialJitter.R    = 0;                                      % jitter the position of retinal cells

   % Mean layer rate 
   lambda.I           = 15;
   lambda.R           = 15;                                     % mean layer rate
   Ra.R               = 0;

   noisedist          = 'poisson';                              % noise distribution 
   noiseparams        = lambda.I;                               % Poisson lambda is mean & var of noise

   %% LEARNING PARAMETERS
   % Variability of Q (covariance matrix)
   b                  = 1;                                      
   beta               = 10;                                     % Controls cov variance
   delta              = 1;                                      % Grid resolution
   N                  = [];
   r                  = [];
   conndist           = [];
   wgtdist            = [];
   wgtparams          = [];
   plastic            = [];
   deldist            = [];
   delparams          = [];
   ka                 = [];
   kb                 = [];
   n                  = [];
   k1                 = [];
   k2                 = [];
   Rb                 = [];

   layerconfig        = v2struct(T,dt,initTime,interval,intervalQ,normQ,...
                                layernames, layerconns, inlayer, retina,...
                                sz,spatialJitter,N,r,conndist,lambda,...
                                noisedist,noiseparams,RF,...
                                wgtdist, wgtparams,b,beta,delta,...
                                psptype,pspparams,pspdelay,psp,...
                                plastic,deldist,delparams,layerdist,...
                                ka,kb,n,k1,k2,Rb,Ra,...
                                endpoints,runSim,evolveSep,...
                                estCov,initCovWithAnal,calcCorr);

   outstruct.outQ        = [];
   outstruct.outmu       = [];
   outstruct.outrates    = [];
   outstruct.outtime     = [];
   outstruct.outweights  = [];
   outstruct.ntwkconfig  = [];
   outstruct.config_file = 'generateRetinaKernelBias';
   outstruct.outfile     = '';
   outstruct.layerconfig = layerconfig;

   %    ntwkconfig  = v2struct(rates, weights, cell_loc, conns, synloc, ...
   %                           totaldelay, maxdelay, time, delayinds, ...
   %                           initTime, interval, intervalQ, ...
   %                           Q, sig, mu, co_inputs...
   %                           );
   ntwkconfig.totaldelay  = 2;
   ntwkconfig.intervalQ   = 0;
   ntwkconfig.synloc      = 0;
   ntwkconfig.rates       = [];
   ntwkconfig.conns       = [];
   ntwkconfig.weights     = [];
   ntwkconfig.maxdelay    = [];
   ntwkconfig.time        = [];
   ntwkconfig.delayinds   = [];
   ntwkconfig.maxdelay    = [];
   ntwkconfig.initTime    = [];
   ntwkconfig.interval    = [];
   ntwkconfig.Q           = [];
   ntwkconfig.sig         = [];
   ntwkconfig.mu          = [];
   ntwkconfig.co_inputs   = [];
   ntwkconfig.cell_loc.R  = zeros(NR^2, 2);
   outstruct.ntwkconfig   = ntwkconfig;

   % tuning = processVisualCortexInput(outstruct,1,[],...
   %             {'rectify',false,'max_phases',16,'num_freqs',16,'max_freq',0.3}); 
   tuning = processVisualCortexInput(outstruct,1,[],...
               {'rectify',false, 'max_phases',16, 'num_angles',10, 'spatial_freqs', sfreqs}); 

   curves = distributionOSbias(tuning, 'R', 'levick', 'max', 'start', 1, 0);
   [nr, nc, next] = getNumSubplots(gcf); 
   % get num row & cols of figure, & just hope the last one is free!
   if ~isempty(next)
      subplot(nr,nc,next)
         bar(Lbins,Lbars/sum(Lbars)); title('Levick');
   else
      subplot(nr,nc,nr*nc),
         hold on; 
         bh = bar(Lbins,Lbars/sum(Lbars)); % title('Levick'); 
         set(bh, 'facecolor', 'r');
         legend(bh, 'levick');
   end
end

if plotresults
   % results cover all spatial freqs + cell's optimal spatial frequency
   freqs = tuning.spatial_freq; Nf = length(freqs);
   mu    = zeros(Nf+1,1);
   v     = zeros(Nf+1,1);
   for fi=1:(Nf+1)
      mu(fi,1) = mean(curves.start.R.bias(:,fi));
      v(fi,1)  = var(curves.start.R.bias(:,fi));
   end
   figure; 
   ah = subplot(1,2,1);
      % plot mean of bias at each spatial freq, & prepend bias calculated for 
      % each cell's preferred freq
      a1 = plot([0; freqs(:)], [mu(end); mu(1:end-1)], 'b.', 'linewidth',4, 'markersize',20);
      lab = get(ah, 'xticklabel');
      set(ah, 'xticklabel', ['Pref freq'; lab(2:end)]);
      hold on; 
      a2 = plot(0, Lmean, 'rx');
      xlabel('Orientation bias - mean');
      ylabel('Proportion of retinal cells');
      legend( a2, 'Levick' );
   ah = subplot(1,2,2);
      a1 = plot([0; freqs(:)], [ v(end);  v(1:end-1)], 'b.', 'linewidth',4, 'markersize',20);
      lab = get(ah, 'xticklabel');
      set(ah, 'xticklabel', ['Pref freq'; lab(2:end)]);
      hold on; 
      a2 = plot(0, Lvar, 'rx');
      legend( a2, 'Levick' );
      xlabel('Orientation bias - var');
      ylabel('Proportion of retinal cells');
end
         
         
         
         
         
         
         
         