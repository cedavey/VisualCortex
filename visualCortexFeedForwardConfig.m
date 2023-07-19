%   layerconfig  = visualCortexNoInhibConfig(fname,verbose)
%
% Input a filename to save the parameter configuration to. This file will
% be of type .mat.
% Set verbose to true to print warnings and error messages when parsing the
% network configuration (defaults to true).
function [layerconfig, fname]  = visualCortexFeedForwardConfig(fname,verbose)

%% Parameter configuration without lateral inhibitory connections

   if nargin<2,verbose= false; end
   if nargin<1 || isempty(fname)
     fname = 'visual_cortex_paramConfig'; 
     cprintf('*comments',sprintf('\nSaving visual cortex parameter configuration to %s\n',fname));
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
   % if epsp is delta make it empty so we can identify it easily in calcs.
   if strcmpi(psptype,'delta')
      psp             = [];
      pspdelay        = 1;
   else
      psp             = EPSP(psptype,pspparams,dt,pspdelay*dt);
   end
%    PSP = {'psptype',psptype,'pspparams',pspparams,'pspdelay',pspdelay,'psp',psp};

   %% LAYER PARAMETERS
   layernames         = {'A', 'B', 'C'};                   % neuron layers 
   layerconns         = {'A','B'; ...
                         'B','C'; ...
                        };     % layer cxs
   inlayer            = 'I'; % name of input layer for readability
   retina             = 'A'; % name of retina layer for readability

   layerSz            = @(r) ceil([r r]*4);
   sz.I               = [100 100];                              % input image
   sz.A               = [ 10  10];                              % retina
   sz.B               = [ 30  30];                              % lgn
   sz.C               = [  5   5];                              % cortex
   spatialJitter.A    = 0;                                      % jitter the position of retinal cells
   spatialJitter.B    = 0;
   spatialJitter.C    = 0;

   % Mean layer rate 
   lambda.I           = 15;
   lambda.A           = 15;                                     % mean layer rate
   lambda.B           = 15;
   lambda.C           =  0;
   
   % Receptive field parameters
   RF.kernel_size     = 10;
   RF.sd_surr_hor     = 7; 
   RF.sd_surr_ver     = 7;
   RF.sd_cent_hor     = 8;
   RF.sd_cent_ver     = 5; 
   RF.coeff_cent      = 6.00e-4;
   RF.coeff_surr      = 1.30e-4;
   RF.ratio_surr_cent = 2;
   RF.normfn          = 'power';
   RF.retina_kernel   = retina_kernel(RF,false);                 % receptive field kernel
   RF.input_var       = 5;
   [RF.angleSet, RF.angleProb] = catRetinaOS(10);
   RF.angleSet        = num2cell(RF.angleSet);                  % angle options, in cell to allow 'radial'
%    RF.angleSet        = {0,90,'radial'};                      % angle options
%    RF.angleProb       = [1/3 1/3 1/3];                        % prob of each angle
%    dangle             = 22.5;
%    angleSet           = 0:dangle:(180-dangle);
%    RF.angleSet        = mat2cell(angleSet,1,ones(length(angleSet),1)); % angle options, in cell to allow 'radial'
%    RF.angleProb       = ones(length(RF.angleSet),1)/length(RF.angleSet); % prob of each angle
   noisedist          = 'poisson';                              % noise distribution 
   noiseparams        = lambda.I;                               % Poisson lambda is mean & var of noise
   
   %% LEARNING PARAMETERS
   % Variability of Q (covariance matrix)
   b                  = 1;                                      
   beta               = 10;                                     % Controls cov variance
   delta              = 1;                                      % Grid resolution

   %% Connection specific weight update parameters
   %% Connection specific weight update parameters
   pre = 'A'; post = 'B'; cxlabel = [pre post];
   plastic.(cxlabel)  = false;                                  % plasticity on conns
   % Delay distributions
   deldist.(cxlabel)  = 'constant';                             % spatial dist of axonal delay
   delparams.(cxlabel)= 0;                                      % parameters for delay dist
   % M_pre presyn cells, each with N conns to each of M_post postsyn cells
   % --> N*M_pre / M_post connection per postsyn cell
   N.BA               = round(prod(sz.B)/prod(sz.A));
   N.(cxlabel)        = max([1 round(prod(sz.A)*N.BA/prod(sz.B))]);
   r.(cxlabel)        = 1.5;                                    % radius of Gaussian connection density
   Rb.(cxlabel)       = b/N.(cxlabel);                          % Rate input scale
% Rb.(cxlabel) = 1;
   conndist.(cxlabel) = 'fixed_norm';                           % spatial dist of connections
   wgtdist.(cxlabel)  = 'constant';                             % initialising distribution of cx weights
   wgtparams.(cxlabel)=  1;                                     % initialising weight value range
   n.(cxlabel)        =  1;                                     % E/I proportions

   pre = 'B'; post = 'C'; cxlabel = [pre post];
   plastic.(cxlabel)  = true;                                   % plasticity on conns
   deldist.(cxlabel)  = 'constant';                             % spatial dist of axonal delay
   delparams.(cxlabel)= 0;                                      % parameters for delay dist
   N.(cxlabel)        = 150;
   r.(cxlabel)        = 5;                                      % radius wrt presyn layer  
   Rb.(cxlabel)       = b/N.(cxlabel);                          % Rate input scale
% Rb.(cxlabel) = 1;
   conndist.(cxlabel) = 'fixed_norm';
   wgtdist.(cxlabel)  = 'constant';
   wgtparams.(cxlabel)=  0.5;
%    a.(cxlabel)        = 3;                                    % Radius ratio btwn layers
   n.(cxlabel)        = 1;                                      % E/I proportions
   % plasticity parameters 
   ka.(cxlabel)       = 5e-4;                                   % Weight background
   kb.(cxlabel)       = 5e-2;                                   % Weight cov scale
   k1.(cxlabel)       = 1.0;                                    % Homeostatic equilibrium
   k2.(cxlabel)       =-5.0;                                    % First order wgt dependence
   
   layerconfig        = v2struct(T,dt,initTime,interval,intervalQ,normQ,...
                                layernames, layerconns, inlayer, retina,...
                                sz,spatialJitter,N,r,conndist,lambda,...
                                noisedist,noiseparams,RF,...
                                wgtdist, wgtparams,b,beta,delta,...
                                psptype,pspparams,pspdelay,psp,...
                                plastic,deldist,delparams,layerdist,...
                                ka,kb,n,k1,k2,Rb,...
                                endpoints,runSim,evolveSep,...
                                estCov,initCovWithAnal,calcCorr);
                       
   layerconfig = parseVisualCortexConfigFile(layerconfig,verbose);

	save(fname,'layerconfig');
   
end
