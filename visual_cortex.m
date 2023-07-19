%
% [outstruct, outfile] = visual_cortex(config_file, doplot, fprefix, configFileHandle)
%
%     to evolve network weights according to specification file (e.g.
%     @visualCortexWithInhibConfig )
%
% [outstruct, outfile, tuning] = visual_cortex(config_file, doplot, fprefix, configFileHandle, genTuningCurves)
%
%     to evolve network weights & then generate tuning curves
%
%                                visual_cortex(outstruct, 1); % to plot previous results
%
%     to simply replot evolved network 
%
% [outstruct, outfile] = visual_cortex(outstruct, doplot, fprefix) 
%
%     to continue/extend simulation, requires modifying the simulation 
%     duration first, e.g. 
%        outstruct.layerconfig.T = 1000;
%        outstruct.layerconfig.continueSim = true;
%
% Simulate the visual cortex according to configuration parameters in 
% config_file, if they are pre-saved, or extract the config parameters from 
% the config function provided in configFunctionHandle (defaults to 
% @visualCortexNoInhibConfig). 
%
% To continue a simulation that hasn't converged yet type:
% [outstruct, outfile] = visual_cortex(outstruct, true, outfile)
% Inputs:
%   config_file      - Name of .mat file containing configuration params.
%                      If configFunctionHandle is provided too, then config
%                      parameters are read from the function, & saved to
%                      config_file. If no function handle is provided then
%                      try to load config_file. If this fails, exit.
%   fprefix          - file prefix for simulation output .mat file, which 
%                      will be named according to primary config values
%   configFileHandle - name of network configuration fn if not using .mat
%   genTuningCurves  - generates tuning curves & plots for plastic layers
%
% You can also use the function to plot previously generated data by 
%        visual_cortex(outstruct, doplot)
function varargout = visual_cortex(config_file, varargin)
% TO DO:
% - implement each layer in frequency rather than spatially
% - f1 function gives weird increase at higher frequencies?
% - ask re Sagar's orthogonal dvlpmt thing - need recurrent input for this?
% - I have changed the code to only allow a single contribution to Ra (line 173)
% - could f1 ever work given it's amplitude? Can't subtract inhib from tuning curve 
% - 

% - change to using natural input by making noisedist and noiseparams based
%   on natural input, i.e. outstruct.layerconfig.noisedist = 'natural', and
%   outstruct.layerconfig.noiseparams = 50 (number of timesteps to show
%   each image for). Then if you're continuing a simulation you need to
%   increase T, i.e. outstruct.layerconfig.continueSim = true,
%   outstruct.layerconfig.T = 500. All together this becomes (from MEMORY!):
%        outstruct.layerconfig.noisedist   = 'natural';
%        outstruct.layerconfig.noiseparams = 50;
%        outstruct.layerconfig.continueSim = true;
%        outstruct.layerconfig.T = 500;
% - average sum(weight*input) from inhib compared to excitatory conns?
% - can we use Rb to scale the excitatory minus inhibitory blip?
% - normalise f1 by f0 (complex versus simple cells)
% - check what is saved in the output structures saved to file
% - add different connection distributions - I thought I had already??!
% - colour code final weight squares wrt retinal cell synapse emerged from
% - does strength of orientation preference vary with kernel angle
% - change layer declaration so retina is just another layer, but have a
%   continuous/DT setting that states how to process 
% - what sort of temporal filtering should be in the input to the retina? 
% - determine the most appropriate rate update equation
% - make N.B a function of ratio from pre to post & desired connections
%   from pre to post
     outstruct = []; outfile = []; varargout = cell( nargout, 1 );
    [layerconfig, doplot, fprefix, genTuningCurves] = parseVisualCortexInputs(config_file, varargin);
    if isempty(layerconfig)
       fprintf('visual_cortex: error parsing inputs, exiting...\n\n');
       return;
    end
    % only prepare network if runnimg sim (not if just plotting results)
    if layerconfig.runSimInputs && ~layerconfig.continueSim
       [ntwkconfig, layerconfig] = prepareVisualCortexNetwork(layerconfig);
       if isempty(ntwkconfig)
          return;
       end
       
    else % if layerconfig.continueSim
       outstruct = config_file;
      % If continuing a previous sim, update the ntwk time vector to new
      % end time, which is retrieved from layerconfig.T
      ntwkconfig = config_file.ntwkconfig;
      
      % if continuing simulation to input natural images after noise, we
      % need to make sure that we're converted using parameter
      % configuration to the parameter format that the functions expect
      if strcmpi( layerconfig.noisedist, 'natural' )
         noiseparams = naturalImageParams( layerconfig, layerconfig.noiseparams );
         if isempty( noiseparams )
            return;
         end
         layerconfig.noiseparams = noiseparams;    % udpate layerconfig struct with changes
      end

    end
   
   % Extract necessary params from the config file
   [T,          dt,             calcCorr,                   ...
    initTime,   runSim,         runSimInputs,   continueSim,...
    layernames, layerconns,     inlayer,                    ...
    retina,     retina_on,      retina_off,     RF,         ...
    sz,         N,              lambda,                     ...
    noisedist,  noiseparams,                                ...
    plastic,    wgtparams,      normQ,                      ...
    psptype,    pspdelay,       psp,                        ...
    deldist,    delparams,      layerdist,    velocity      ...
    kb, n, k1, k2, Ra, Rb,                                  ...
    endpoints,  evolveSep,      thresholdWgts,              ...
    estCov,     initCovWithAnal] = struct2v(layerconfig,...
                       'T',         'dt',         'calcCorr',                  ...
                       'initTime',  'runSim',     'runSimInputs','continueSim',...
                       'layernames','layerconns', 'inlayer',                   ...
                       'retina',    'retina_on',  'retina_off',  'RF',         ...
                       'sz',        'N',          'lambda',                    ...
                       'noisedist', 'noiseparams',                             ...
                       'plastic',   'wgtparams',  'normQ',                     ...
                       'psptype',   'pspdelay',   'psp',                       ...
                       'deldist',   'delparams',  'layerdist',   'velocity',   ...
                       'kb', 'n', 'k1', 'k2', 'Ra', 'Rb',                      ...
                       'endpoints', 'evolveSep',  'thresholdWgts',             ...
                       'estCov',    'initCovWithAnal');
                    
   %% Run simulation
   if runSim && runSimInputs % run sim from param config & inputs are valid
      [time, initTime, interval, intervalQ,...
       rates, weights, cell_loc, synloc, conns,...
       delayinds, totaldelay, maxdelay,...
       co_inputs Q, sig, mu] ...
                 = struct2v(ntwkconfig,'time',      'initTime',   'interval', 'intervalQ',...
                                       'rates',     'weights',    ...
                                       'cell_loc',  'synloc',     'conns',    ...
                                       'delayinds', 'totaldelay', 'maxdelay', ...
                                       'co_inputs', 'Q',          'sig',      'mu');

      nlayers = length(layerconfig.layernames);  % number of layers
      nconns  = size(layerconns,1);  % number of layer cxs
      
      % If we're evolving plastic layers separately we need multiple
      % iterations of the simulation, where a plastic layer is only included
      % when all layers before it have evolved to maturity
      numit   = ternaryOp(evolveSep, sum(cell2mat(struct2cell(plastic))),1);
      evolve  = cell2mat(struct2cell(plastic)); % indices of plastic layers to evolve

      %% Extending existing simulation
      % If continuing simulation, copy last output to current rates etc
      if continueSim
         [outQ, outmu, outsig, outrates, outtime,outweights,co_outputs] ...
                    = struct2v(outstruct,'outQ', 'outmu', 'outsig', 'outrates',...
                                         'outtime','outweights','co_outputs');
         ntwkconfig.time = 0:dt:layerconfig.T; % extend time to continue simulation
         time = ntwkconfig.time;
         % Rotate rates & output if requested
         for li=1:nlayers
            post = layernames{li};
            % Rotate rates (keeping recent history for axonal delays)
            rates.(post)(end,:) = outrates.(post)(end,:);
         end
         % get last recorded weights
         for ci=1:nconns
            pre     = layerconns{ci,1};        % presynaptic layer name
            post    = layerconns{ci,2};        % postsynaptic layer name
            cxlabel = [pre post];              % connection label/name
            if plastic.(cxlabel)
               for i=1:prod(sz.(post))
                  weights.([pre post]){i} = outweights.([pre post]){i}(end,:)';
                  if estCov
                     Q.(cxlabel){i}   = squeeze(outQ.(cxlabel){i}(end,:,:));
                     sig.(cxlabel){i} = toVec(outsig.(cxlabel){i}(end,:,:));
                     mu.(cxlabel){i}  = toVec(outmu.(cxlabel){i}(end,:,:));
                  end
               end
            end
            for cc=1:length(co_inputs.(cxlabel))
               co_layer = co_outputs.(cxlabel){cc};
               [Q_co,sig_in,sig_out,mu_in,mu_out] = struct2v(co_layer,...
                     'Q','sig_in','sig_out','mu_in','mu_out');
                co_inputs.(cxlabel){cc}.Q = squeeze(Q_co(end,:,:));
                co_inputs.(cxlabel){cc}.sig_in  = squeeze(sig_in(end,:,:));
                co_inputs.(cxlabel){cc}.sig_out = squeeze(sig_out(end,:,:))';
                co_inputs.(cxlabel){cc}.mu_in   = squeeze(mu_in(end,:,:));
                co_inputs.(cxlabel){cc}.mu_out  = squeeze(mu_out(end,:,:))';              
            end
         end
         % If iterating each layer separately & continuing sim, can only continue last iteration
         startit = numit; 
         
      %% New simulation
      % If starting new sim, initialise output structures
      else
         % Initialise output structures
         [outrates,outweights,outtime,outQ,outsig,outmu,co_outputs] = ...
                      initialiseVisualCortexOutput(layerconfig,interval,intervalQ,rates,weights,co_inputs);
         if ~estCov
            outQ = Q;
         end
         startit = 1; % start iteration number is for when continueSim is re-implemented
      end

      %% Simulate
      tic;

      %% Iterate learning over layers
      % For each iteration - if allowing each layer to evolve separately
      for it=startit:numit
         % Get  plastic fields for this iteration, & include only layers up to
         % and including the plastic layer
         if evolveSep && ~continueSim
            pInd   = find(evolve,1,'first');
            layers = layerconns(1:pInd);
         % If evolving simultaneously then use all layers straight away
         else
            layers = layerconns; 
         end
         if ~continueSim
            starttime = (-ceil(initTime));
            outi = double(endpoints); % reset rates/weights output index
            outq = double(endpoints); % reset cov/mu output index
            outr = double(endpoints); % starting output rates index depends on if we're incl. endpoints
            
         else
            % get last timestep & output indices from old last output, 
            % scaled by increase in total time
            ti   = outtime(end)/dt; % last timestep of previous simulation
            outi = floor(layerconfig.interval  * outtime(end) / T);
            outq = floor(layerconfig.intervalQ * outtime(end) / T);
            outr = floor(layerconfig.interval  * outtime(end) / T);
            % update interval btwn output samples - user gives num samples
            % required, so calc time btwn samples
            interval  = layerconfig.interval; intervalQ = layerconfig.intervalQ;
            interval  = round(ternaryOp(interval>0,round(T/dt)/interval,Inf));  % interval between rate/wgt outputs
            intervalQ = ternaryOp(intervalQ>0,round(T/dt)/intervalQ,Inf);% interval between cov/mu outputs
            starttime = ti; % redo last timestep just in case continuing from error
            % udpate outtime - don't do it before we get current ti from it
            outtime   = getLinskerTimeVec(T,dt,interval,endpoints,intervalQ);   
         end

         %% Simulate learning - for each timestep
         prevprog = []; % progress to user
         for ti=starttime:length(time)
%             if ti==0 && initTime>0, printTime(toc,'Initialising covariance matrices took '); end
            % determine if this is an output timestep & update progress to user
            if ti>=0 && (mod(ti,interval)==0 || (endpoints && (ti==length(time) || ti==0)))
               prog = sprintf( '\t%.2f percent complete\n', (ti-starttime)/length(time)*100 );
               % refreshdisp( prog, prevprog );
               cprintf( 'keywords*', prog );
               prevprog = prog;
               outr     = outr + 1;
               output   = true;
            else
               output   = false;
            end
            % If we're in the simulation (ti>0) and this is a cov output step
            % If no cov output required, intervalQ will be Inf
            if ti>=0 && (mod(ti,intervalQ)==0 || (~isinf(intervalQ) && endpoints && (ti==length(time) || ti==0)))
               outq    = outq + 1;
               outputQ = true && estCov;
            else
               outputQ = false;
            end

            %% Update layer rates
            % When updating retina from input, use delayed time samples, so
            % that v(t) = v(t-1) + sum(RF*input(t-1))
            
            % if input type is natural image, provide timestep as a
            % parameter so we can tell whether to switch images
            if strcmpi(noisedist, 'natural')
               noiseparams = updateNaturalImageParams(noiseparams, ti); 
               if isempty(noiseparams{end})
                  prog = sprintf('saved everything in tmp coz no image so exiting\n');
                  cprintf('Keywords*', prog);
                  return;
               end
            end
            rates.(inlayer)(end,:)       = toVec( generateVisualCortexInput( lambda.(inlayer), sz.(inlayer), noisedist, noiseparams) );
            if ~isempty(retina)
               rates.(retina)(end,:)     = toVec( generateRetinaInput(sz.(retina), sz.(inlayer), [0 0], ...
                                                   cell_loc.(retina), Ra.(retina), rates.(inlayer), RF) );
            end
            if ~isempty(retina_on)
               rates.(retina_on)(end,:)  = toVec( generateRetinaInput(sz.(retina_on), sz.(inlayer), [0 0], ...
                                                   cell_loc.(retina_on), Ra.(retina_on), rates.(inlayer), RF, false ) );
            end
            if ~isempty(retina_off)
               rates.(retina_off)(end,:) = toVec( generateRetinaInput(sz.(retina_off), sz.(inlayer), [0 0], ...
                                                   cell_loc.(retina_off), Ra.(retina_off), rates.(inlayer), RF, true) );
            end
            for ci=1:nconns
               pre     = layerconns{ci,1};        % presynaptic layer name
               post    = layerconns{ci,2};        % postsynaptic layer name
               cxlabel = [pre post];              % connection label/name
               for i=1:prod(sz.(post))
                  % To allow recurrent connections we want inputs to be from last time
                  % sample, else if there are recurrent connections & we calculate inputs
                  % using current time point, then when we calculate inputs going the
                  % other direction in the recurrent connections, they might get half
                  % built inputs. So treat last time step as rates that we're still
                  % building, calculating each group's input to the rate one at a time,
                  % and the previous timestep is the current, fully built, rate.
                  inrates = getVisualCortexRateWithDelay(rates.(pre)(1:end-1,:),...
                                       conns.(cxlabel){i}, delayinds.(cxlabel){i},...
                                       psp, dt);

                  % Update rates: F_j = Ra + Rb*sum(c_ij*F_i)
                  rates.(post)(end,i) =                                    ...
                              updateVisualCortexRates(Rb.(cxlabel),        ...
                                                      inrates,             ...
                                                      weights.(cxlabel){i},...
                                                      rates.(post)(end,i)); % need for multiple inputs into layer
                  if isnan(rates.(post)(end,i))
                     cprintf('*Errors',sprintf('Timestep %d:\tLayer %s rates, neuron %d is nan\n',...
                                              ti, post, i));
                  end
               end
            end
            
            %% Update cov estimates & layer plasticity
            % Gotta do this after rate updates because weight update may
            % need updated rates for later layers in connection list
            for ci=1:nconns
               pre     = layerconns{ci,1};        % presynaptic layer name
               post    = layerconns{ci,2};        % postsynaptic layer name
               cxlabel = [pre post];              % connection label/name

               % If this layer is plastic & provides joint input into 
               % another layer, we need to calc. covariance between 
               % the 2 input layers for wgt udpate. This update is not 
               % per postsyn neuron, so take outta for each postsyn loop
               for cc=1:length(co_inputs.(cxlabel))
                  co_layer = co_inputs.(cxlabel){cc};
                  [Q_co, sig_in, sig_out, mu_in, mu_out, r_co, colabel] = struct2v(co_layer, ...
                        'Q','sig_in','sig_out','mu_in','mu_out','r','colabel');
%                   [Q_co,mu_co,sig_co,r_co,conns_co,colabel,jt_input,delay_co] = struct2v(co_layer,...
%                         'Q_co','mu_co','sig_co','r_co','conns_co','colabel','co_pre','delayinds');
                  co_pre  = colabel(1); co_post = colabel(2);
                  for i=1:prod(sz.(co_post))
                     inrates = getVisualCortexRateWithDelay(rates.(co_pre)(1:end-1,:),...
                                          conns.(colabel){i}, delayinds.(colabel){i},...
                                          psp, dt);
                     outrate = rates.(co_post)(end-1,i)';

                     if estCov && ti>=1
                         cind = conns.(colabel){i};
                        [Q_co(cind,i), sig_in(cind,i), sig_out(1,i), mu_in(cind,i), mu_out(1,i)] = ...
                              calculateBetweenLayerCov(Q_co(cind,i),             ...
                                                    sig_in(cind,i), sig_out(1,i),...
                                                    mu_in(cind,i),  mu_out(1,i), ...
                                                    inrates,     outrate,     ...
                                                    ti, dt, normQ);
                     end
                  end
                  % Update co_layer info - v2struct has some sort of error, so 
                  % requires a variable rather than a string for struct name input
                  nameOfStruct2Update = 'co_layer';
                  co_layer = v2struct(Q_co,sig_in,sig_out,mu_in,mu_out,nameOfStruct2Update,...
                                      {'fieldNames','Q','sig_in','sig_out','mu_in','mu_out'});
                  co_inputs.(cxlabel){cc} = co_layer;
                  if outputQ
                     co_out = co_outputs.(cxlabel){cc};
                     co_out.Q(outq,:,:)      = Q_co;
                     co_out.mu_in(outq,:,:)  = mu_in;
                     co_out.mu_out(outq,:,:) = mu_out;
                     co_out.sig_in(outq,:,:) = sig_in;
                     co_out.sig_out(outq,:,:)= sig_out;
                     co_outputs.(cxlabel){cc} = co_out;
                  end
               end
               % Update plastic connections btwn pre & post
               for i=1:prod(sz.(post))
                  inrates  = getVisualCortexRateWithDelay(rates.(pre)(1:end-1,:),...
                                       conns.(cxlabel){i},delayinds.(cxlabel){i},...
                                       psp,dt);

                  % All other layers calculate rates from Linsker's update equation
                  if estCov && plastic.(cxlabel) && (ti-starttime > totaldelay)
                     [Q.(cxlabel){i},sig.(cxlabel){i},mu.(cxlabel){i}] =   ...
                                  calculateVisualCortexCov(Q.(cxlabel){i}, ...
                                                      sig.(cxlabel){i},    ...
                                                      mu.(cxlabel){i},     ...
                                                      inrates,             ...
                                                      ti-starttime+1, normQ);


                  end
                  if outputQ && plastic.(cxlabel)
                     outQ.(cxlabel){i}(outq,:,:) = Q.(cxlabel){i};
                     outsig.(cxlabel){i}(outq,:) = sig.(cxlabel){i};
                     outmu.(cxlabel){i}(outq,:)  = mu.(cxlabel){i};
                  end

                  %% Update weights - plasticity
                  if ti>0 && plastic.(cxlabel) && ti>1
                     % cdot_ij = k1 + 1/N*sum((Q_ij + k2)*c_ij
                     weights.(cxlabel){i} =                                ...
                            updateVisualCortexWeights(k1.(cxlabel),        ...
                                                      k2.(cxlabel),        ...
                                                      N.(cxlabel),         ...
                                                      weights.(cxlabel){i},...
                                                      Q.(cxlabel){i},      ...
                                                      conns.(cxlabel){i},  ...
                                                      dt);
                     if any(isnan(weights.(cxlabel){i}) | isinf(weights.(cxlabel){i}))
                        cprintf('*Errors',sprintf('Timestep %d:\tConn %s weight %d is nan\n',...
                                                 ti, cxlabel, i));
                     end

                     for cc=1:length(co_inputs.(cxlabel))
                     % Update weights for plastic connections based on
                     % weights of other input to postsynaptic layer,
                     % scaled by covariance btwn co-inputs to postsyn layer.
                        co_layer    = co_inputs.(cxlabel){cc};
                        colabel     = co_layer.colabel;
                        other_input = co_layer.other_input;
                        weights.(cxlabel){i} =                             ...
                              updateJointInputWeights(co_layer.k2,         ...
                                                  N.(other_input),         ...
                                                  weights.(cxlabel){i},    ...
                                                  weights.(other_input){i},...
                                                  co_layer.Q,              ...
                                                  conns.(cxlabel){i},      ...
                                                  conns.(colabel),         ...
                                                  conns.(other_input){i},  ...
                                                  n.(cxlabel), dt);
                        if any(isnan(weights.(cxlabel){i}) | isinf(weights.(cxlabel){i}))
                           cprintf('*Errors',sprintf('Timestep %d:\tConn %s colayer %s weight %d is nan\n',...
                                                    ti, cxlabel, colabel, i));
                        end

                     end
                  end
               end
               % Threshold connections after all weight updates have been
               % applied from all input connections, otherwise the inhib
               % may never get a look in!
               for i=1:prod(sz.(post))
                  if ti>0 && thresholdWgts && plastic.(cxlabel)
                     % Upper & lower bound weights using n from Linsker
                     weights.(cxlabel){i}(weights.(cxlabel){i}> n.(cxlabel)   ) = n.(cxlabel);
                     weights.(cxlabel){i}(weights.(cxlabel){i}<(n.(cxlabel)-1)) = n.(cxlabel)-1;
                     if all( weights.(cxlabel){i} == 0 )
                        prog = sprintf('All %s(%d) weights are 0, pausing for inspection', cxlabel, i);
                        cprintf( 'SystemCommands*', prog);
                        prog = 'extra line to pause on';
                     end
                  end
               end
            end % for each post-synaptic neuron

            % Rotate rates & output if requested
            for li=1:nlayers
               lname = layernames{li};
               % Copy rates to output structure
               if output
                  outrates.(lname)(outr,:) = rates.(lname)(end,:);
               end
               % Rotate rates (keeping recent history for axonal delays)
               rates.(lname)(1:end-1,:) = rates.(lname)(2:end,:);
               rates.(lname)(end,:)     = Ra.(lname); % lambda.(post); % initialise next timestep
            end
            % Now copy weight outputs to output structs, if it's an output 
            if output
               for ci=1:nconns
                  pre     = layerconns{ci,1};
                  post    = layerconns{ci,2};
                  cxlabel = [pre post];
   %                outrates.input(outi,:) = rates.input(end,:);
                  if plastic.(cxlabel)
                     for i=1:prod(sz.(post))
                        outweights.(cxlabel){i}(outr,:) = weights.(cxlabel){i};
                     end
                  end
               end
            end
         end % end of for each timestep
         outr = min([outr length(outtime)]); % may have incremented once too often on final timestep

         % Now that this iteration is complete, make the last plastic layer
         % non-plastic, & remove it from the evolve vect
         if evolveSep
            evolve(pInd) = false; % remove layer we just did from evolve list
            pInd   = find(evolve,1,'first'); % find next closest layer to evolve
         % if evolving simultaneously then use all layers straight away
         else
            layers = layerconns; 
         end
      end % for each simulation iteration (to evolve layers separately)
      printTime(toc,'Simulating Linsker layers took ');
      if doplot && ti<length(time)
        try
          plot_visual_cortex; % update weight figures as simulation evolves
        catch ME
          cprintf('*comment','\n\tError plotting final visual cortex figures (%s)\n',...
                  ME.message);
        end
      end

      % Prep output - don't need current weights, Q etc, cuz have output versions
      % (actually, we do if we want to be able to continue the simulation)
%       ntwkconfig = v2struct(time,interval,intervalQ,...
%                             cell_loc,conns,weights,delayinds,...
%                             initTime,totaldelay,maxdelay,synloc,co_inputs);

      % Create filename & save pertinent variables
      optionsStr = ternaryOp(strcmpi(psptype,'delta'),[],['_' psptype '_PSP' num2str(pspdelay)],'');
      optionsStr = [optionsStr ternaryOp(evolveSep,'_evolveSep','')];
      connsStr   = strjoin_CD(cell2mat(layerconns),'_');
      outfile    = sprintf('%s_%s_sz%d%s',fprefix,...
                           connsStr,sz.(inlayer)(1),...
                           ternaryOp(isempty(optionsStr),'',['_' optionsStr]));
      outstruct  = v2struct(outQ,outmu,outsig,outrates,outtime,outweights,co_outputs);
      layerconfig.continueSim = false;  % if we were continuing set it to false
      outstruct.layerconfig   = layerconfig;
      outstruct.ntwkconfig    = ntwkconfig;
      outstruct.config_file   = config_file;
      outstruct.outfile       = outfile;

      cprintf('*comments',sprintf('\nSaving output to %s\n\n',outfile));
      save(outfile,'outstruct','-v7.3');

      if nargout>=1 
         varargout{1} = outstruct; 
      end
      if nargout>=2
         varargout{2} = outfile; 
      end
      
      if genTuningCurves
         % plot postsynaptic plastic layers
         plasticLayers = getPlasticLayers( outstruct );                  
         optargs = { 'genMovies', false, 'iterations', { 'start', 'perc30', 'end' }, ...
                     'num_freqs',  15, 'max_freq', 0.3, 'min_phases', 7, ...
                     'num_angles', 20, ...
                     'responsefn', { 'f1', 'max' }, 'dc_offset', 200, 'amplitude', 100 };

         [tuning, curves] = processVisualCortexInput( outstruct, doplot, plasticLayers, optargs);
         
         if nargout>=3
            varargout{3} = tuning; 
         end
         if nargout>=4
            varargout{4} = curves; 
         end
      end
      

   %% Don't run simulation - just plot results
   else

      [outQ, outmu, outsig, outrates, outtime,outweights,co_outputs] ...
                 = struct2v(outstruct,'outQ', 'outmu', 'outsig', 'outrates',...
                                      'outtime','outweights','co_outputs');
      outr = length(outtime); % rate output index can skip straight to end
      % Need to extract info from output struct & param config structs
      [time,interval,intervalQ,initTime,...
       cell_loc,synloc,conns,...
       delayinds,totaldelay,maxdelay] ...
                 = struct2v(ntwkconfig,'time','interval','intervalQ','initTime',...
                                       'cell_loc','synloc','conns',...
                                       'delayinds','totaldelay','maxdelay');
                                    
   end % end running simulation
   
   % Final plot - goes straight to this if runSim==false
   % If evolving weights in layer connections separately, only 
   % plot weights for final evolution
   if doplot
     try
       plot_visual_cortex; % update weight figures as simulation evolves
     catch ME
        errStr = sprintf('\n\tError plotting final visual cortex results (%s, line %d in %s)\n',...
                          ME.message, ME.stack(1).line, ME.stack(1).name);
       cprintf('*Keywords', errStr);
       return;
     end
   end
end
   















