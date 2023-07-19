%      [tuning, curves] = processVisualCortexInput(outfile)
%      [tuning, curves] = processVisualCortexInput(outstruct)
% OR 
%      [tuning, curves] = processVisualCortexInput(outfile, doplot, plotlayers, optargs)
%
% Use processVisualCortexInput to process spatial gratings for a network
% simulated in visual_cortex. 
% You MUST run visual_cortex prior to this function, and save the output in
% outfile, which is then used by this function, or assign the output to
% outstruct, which can then be input instead.
%
% Inputs:
%  outfile    - name of output file containing output struct from
%               visual_cortex plasticity simulation
%  doplot     - true or false
%  plotlayers - cell array of layers to plot the tuning curves of 
%               e.g. {'B', 'C', 'L'}
%  optargs    - cell array of parameter name and value pairs
%     DEBUG          - boolean (defaults to false)
%     amplitude      - amplitude of input sinusoidal grating (float)
%     dc_offset      - DC offset of input sinusoidal grating (float)
%     genMovies      - boolean (defaults to true)
%     grating_angles - vector of angles to plot in degrees (don't set
%                      num_angles if using this)
%     input_type     - 'grating' (sinusoidal grating), 'wgn' (white Gaussian 
%                       noise), 'wpn' (white Poisson noise)
%     iterations     - cell array or string with time in simulation to plot 
%                      tuning curves for (current options are 'start', 'end', 
%                      or perc<dd> to use the network a certain percentage 
%                      of the way through the simulation, e.g. perc40)
%     max_freq       - maximum spatial frequency to include for gratings
%     max_phases     - if using variable length time for each spatial frequency 
%                      (higher spatial freqs need fewer time points to
%                      capture them), then can set a max here to save time
%                      on lower spatial frequencies
%     min_phases     - if using variable length time for each spatial frequency 
%                      (higher spatial freqs need fewer time points to
%                      capture them), then can set a min here to help
%                      visualise (becomes difficult with only 2 phases!)
%     num_angles     - can set number of orientations to include, so that 
%                      orientation vector evenly samples on [0, 180]
%                      (don't set grating_angles if using this)
%     num_freqs      - can set number of frequencies to include, so that 
%                      frequency vector evenly samples from min to max freq
%                      (don't set spatial_freqs if using this)
%     num_phases     - can set number of phases to include, so that 
%                      number of phases (i.e. time) is constant across all
%                      spatial frequencies (defaults to 1/min_freq)
%     phase_type     - 'variable' if diff spatial freqs have diff number of
%                      sinusoidal phases, or 'constant' to give all spatial 
%                      freqs the same (high freqs have fewer number of
%                      samples in a period, so 'variable' sets using this)
%     rectify        - boolean (defaults to true)
%     responsefn     - response function used to generate tuning curves 
%                      in the call to generateVisualCortexTuningCurves {'f1', 'max'}
%     spatial_freqs  - vector of spatial frequencies to simulate (don't set
%                      num_freqs if using this)
%     stddev         - if input type is WGN you can set the std dev of noise
%
% Output:
%  tuning structure that gives responses to input grating of different
%  spatial frequencies and orientations, which are determined by sampling
%  rate and number of samples, as well as input size. Responses are
%  calculated at requested iterations, drawn from prior to plasticity 
%  (iteration = 'start'), and after plasticity has been run (iteration = 
%  'end'), and a percentage of the way through (e.g. iteration = 'perc40')
%          tuning.iterations.start/perc/end.weights % weights for iteration
%          tuning.iterations.start/perc/end.output  % output rates for each layer
%          tuning.input                        % input grating as movie & matrix
%          tuning.temporal_freq                % temporal phase vector;
%          tuning.spatial_freq                 % spatial frequency vector;
%          tuning.grating_angles               % grating angles vector;
%          tuning.input_phase                  % input phase = time * temp_phase;
%          tuning.time                         % time vector;
% 
% You can also request to plot the tuning curves for specified layers
% 
% See also: generateVisualCortexTuningCurves
function [tuning, curves] = processVisualCortexInput( outfile, doplot, plotlayers, varargin )
%% TO DO
   % Ratio of horizontal to vertical from 0 to 2 (1, 1.5, 2)
   % Save retinal convolution computations to save time 
   % Wavers law - double spatial frequencies each round
   % Rate should be ~100 or so when active
     
   % Are phases appropriate for all spatial frequencies? Currently using
   % the same phases for all --> complete 1 cycle --> always uses index 2
   % of the fft of response, regardless of spatial frequency
   
   tstart = tic;
   curves = []; % return empty curves if we don't make it to the end
   if ischar(outfile)
      load(outfile);
   elseif isstruct(outfile)
      outstruct = outfile; %
   end      

   try
      outstruct2variables;
   catch ME
      str = sprintf('\nprocessVisualCortexInput: cannot open file or struct to get network and layer configs (%s), exiting...\n\n', ME.message);
      cprintf('Errors', '%s',str);
      return;
   end
   % Network defined by:
   %  - sz (layer sizes)
   %  - layernames (layers included in the network)
   %  - conns_loc (location of synapses connecting the network)
   %  - weights (weights of synapses)
   %  - PSP (psp structure)
   %  - delays (delays of synapses)
   %  - lambda (background layer rate)
   %  - RF (receptive field structure)

   % Optional network config parameters
   %  - delayinds (delays in indices - specific to choice of dt)
   %  - synloc (connection indices in subscripts)

   % Notes:
   % Need to run for longer than temporal_freq so that the minimum
   % frequency is not the one we're testing for, else all the orientations
   % with a shorter period will have power put in the temp_freq frequency
   % band, since that's the lowest freq. Need a higher resolution. 
   
   % Note: interpret each input cell as being a timestep (i.e. 1 timepoint,
   % or space point), and each retinal cell as being a single sample.
   % Therefore, dt = number of input cells in each retinal cell. Everything
   % follows from this interpretation. The total number of samples is the 
   % number of cells in the retinal layer, and the maximum frequency is half 
   % of the number of input cells in each retinal cell (look at a line 
   % across, not at the total area). Total number of timesteps is the
   % number of input cells in a line. 
   
   % Note: humans have 120 samples per degree, and the lens filters out
   % spatial variations finer than 60 cycles/degree (i.e. half the max)
   
   input_type     = 'grating'; % 'wpn' % 'wgn' % 
   haveRetina     = ternaryOp(isempty(retina),false,true);
   sz_input       = sz.(inlayer);

   optargs        = parseInputs(varargin, sz, inlayer, lambda);
   [amplitude,  dc_offset,  DEBUG,      genMovies,  grating_angles, ...
    input_type, iterations, max_phases, min_phases, num_angles,    num_freqs, ...
    num_phases, phase_type, rectify,    responsefn, spatial_freqs, stddev] = struct2v(optargs, ...
             'amplitude',      'dc_offset',     'DEBUG',      'genMovies',  ...
             'grating_angles', 'input_type',    'iterations', ... 
             'max_phases',     'min_phases',    'num_angles', 'num_freqs', ...
             'num_phases',     'phase_type',    'rectify',   ...
             'responsefn',     'spatial_freqs', 'stddev');
        
   if haveRetina
      sz_retina= sz.(retina);
      sz.kernel= ones(1,2)*(RF.kernel_size(1));
      ratioIn2Ret = sz_input./sz_retina; % ratio of num inputs to retinal cells
      margin   = ceil(ratioIn2Ret);   % number of input pixels per retinal cell
   else
      margin   = [0 0];
   end
   
   % num phases = number of samples in lowest spatial frequency
   T              = 1;
   temporal_freq  = 1;
   Nfreqs         = length(spatial_freqs);
   temporal_phase = 2*pi*temporal_freq; % want to cover 0 - 2pi --> max_time*freq = 2pi
   visual_input   = cell( length(grating_angles), Nfreqs, 1 );
   visual_output  = cell( length(grating_angles), Nfreqs, 1 );
   visual_weights = cell( length(grating_angles), Nfreqs, 1 );
   interval       = 1; % record output from every iteration

   nlayers        = length(layerconfig.layernames);
   nconns         = size(outstruct.layerconfig.layerconns,1);
   totaldelay     = totaldelay + nconns; % give rates time to propagate down thru connections
   layerconfig.endpoints = false; 
   layerconfig.interval  = interval; 

   fprintf('\tSpatial frequencies (%d):\t', Nfreqs); 
   fprintf('%.2g ', spatial_freqs); fprintf('\n');
   fprintf('\tGrating angles (%d):\t',length(grating_angles)); 
   fprintf('%d ', round(grating_angles)); fprintf('\n');
   fprintf('\tPhases (%d):\t',num_phases);
   if ~isempty(num_phases), fprintf('%.2g ', 0:(T/num_phases):(T-T/num_phases)); end
   fprintf('\n\n');
 
   if DEBUG
      figdebug = figure;
      theta_debug = zeros(num_angles,1);
      for ai=1:length(grating_angles)
         angle = grating_angles(ai);
         theta_debug(ai) = find(outstruct.layerconfig.RF.theta == angle, 1, 'first');
         [theta_sorted,index_debug] = sort(theta_debug);
      end
   end
   if genMovies, fh = figure; end
   
   % for each iteration requested from 'start', 'perc<dd>', and 'end'
   for iti=1:length(iterations)
      for si=1:Nfreqs
         tstart_sfreq = tic;

         sfreq = spatial_freqs(si);
         cprintf('Comments',sprintf('Recording response at %s, freq %g\n',...
                                    iterations{iti},sfreq));
         % Number of phases (i.e. number of time points) depends on
         % frequency
         if strcmpi(phase_type,'variable')
            if sfreq==0
               num_phases = 1;
            else
               num_phases  = min([max_phases 1/sfreq]);
               num_phases  = max([min_phases num_phases]);
            end
         end
         dt             = T/num_phases;
         time           = (0:dt:T)';
         input_phase    = time*temporal_phase;
         layerconfig.T  = T; 
         layerconfig.dt = dt;
         % Initialise output structures for this frequency/number of phases
         outrates       = initialiseVisualCortexOutput(layerconfig,interval,intervalQ);

         % Get weights to use for tuning curves
         for ci=1:length(layerconns)
            pre     = layerconns{ci,1};
            post    = layerconns{ci,2};
            cxlabel = [pre post];
            if plastic.(cxlabel)
               % Switch the name of the iteration to a time index (e.g.
               % start iteration means we want the tuning curves at the
               % start, before plasticity, so get time index of 1.
               switch iterations{iti}
                  case {'start'} % get initial weights
                     tpt = 1;
                  case {'end'}   % get weights at the end of the simulation
                     tpt = size(outweights.(cxlabel){1}, 1);
                     % for some reason some of the final weights were NaN,
                     % so I had to search for the last non NaN & make it's
                     % time index the last time point collected. This is
                     % obviously stupid, & I seem to recall fixing it?!
%                      for i=1:prod(sz.(post))
%                         tpt = size(outweights.(cxlabel){i}, 1);
%                         while any(isnan(outweights.(cxlabel){i}(tpt,:)))
%                            tpt = tpt - 1;
%                         end
%                         weights.(cxlabel){i} = outweights.(cxlabel){i}(tpt,:)';
%                      end
                  otherwise
                     % to get weights a percentage of the way through the
                     % simulation, set the iteration to 'perc40' for example
                     if strfind( iterations{iti}, 'perc' )
                        tmp  = regexp( iterations{iti}, '\d+', 'match' );
                        perc = str2double( tmp{1} );
                        if perc <= 1 || 100 <= perc
                           str  = sprintf( 'Iteration percentage (%d) must be between 0 and 100 \n',...
                                            perc );
                           cprintf( 'Errors*', str );
                           return;
                        end
                        % get index of the requested percentage of time
                        tpt  = round( length(outtime) *perc/100 ); 
                        
                     else
                        str  = sprintf('Unrecognised tuning iteration (%s)',...
                                        iterations{iti});
                        cprintf( 'Errors*', str );
                        return;
                     end
               end
               for i=1:prod(sz.(post))
                  weights.(cxlabel){i} = outweights.(cxlabel){i}(tpt,:)';
               end
               weights.(cxlabel){i}(isnan(weights.(cxlabel){i})) = 0;
            end
         end

         %% For each orientation, generate & record input/output
         for ai=1:length(grating_angles)
            % Generate oriented input
            angle = grating_angles(ai);
            for li=1:length(layernames)
               outrates.(layernames{li}) = zeros(size(outrates.(layernames{li})));
            end

            telapsed_input = 0;
            %% For each phase, generate oriented grating & input to ntwk
            for tphase=1:length(time)
            % For gratings, update phase & re-make input grating if its
            % first iteration, else grab previously made grating
               switch input_type
                  case {'grate','gratings','grating'}
                     inputparams = {tphase*dt, angle, temporal_phase, ...
                                    sfreq, amplitude, dc_offset, rectify};
                  case {'poisson','poiss','exponential','exp','wpn'}
                     inputparams = {};
                  case {'gauss','gaussian','norm','normal','wgn'}
                     inputparams = sig;
               end
               in = generateVisualCortexInput(dc_offset, sz.(inlayer), input_type, inputparams, margin);
               if iti==1 % first iteration
                  if genMovies 
                     figure(fh);
                     marginx   = ceil(margin(1)/2); marginy = ceil(margin(2)/2);
%                      imagesc(flipud((in(:,:,tphase)))); 
                     pcolor( in(marginx+1:(end-marginx), (marginy+1):end-marginy) ); 
%                      pcolor(((in(:,:,tphase)))); 
                     M(tphase) = getframe;
                  end
%                else % grab previously made gratings if not first time through
%                   in(:,:,tphase) = visual_input{ai,si}.input(:,:,tphase);
               end
               
               % initialise rates for each layer - initialise samples to mean
               % layer rate, except final timepoint, which should be set to Ra
               % for postsynaptic layers, since input from previous layer will
               % + Ra gives mean layer rate for timepoint (previous timepoints
               % will be written over as we advance through processing input)
               rates.(inlayer) = zeros(1,prod(sz.(inlayer)+margin*2));
               for li=1:nlayers
                  % Initialise old timesteps to mean rate, & current
                  % timesteps to Ra (background rate)
                  lname = layernames{li};
                  if ~strcmpi(lname, inlayer)
                     rates.(lname) = ones( totaldelay, prod(sz.(lname)) ) * lambda.(lname);
                     rates.(lname)(end,:)= ones( 1, prod(sz.(lname)) ) * Ra.(lname);
                  else
                     rates.(lname) = ones( totaldelay, prod(sz.(lname) + margin*2) ) * lambda.(lname);
                  end
               end
               for ci=1:size(layerconns,1)
                  post = layerconns{ci,2};
                  rates.(post)(end,:)  = ones( 1, prod(sz.(post)) ) * Ra.(post);% record inhibitory input to layer C 
               end

               % Show static, oriented input for entire trial - trial length is
               % equal to the max delay, so we have 1 full cycle
               rates.(inlayer)(1,:) = in(:); % constant input
               % Prepare first timestep
               if haveRetina
                  rates.(retina)(end-1,:)= toVec(generateRetinaInput(sz.(retina), sz.(inlayer), margin, cell_loc.(retina),...
                                                                     Ra.(retina), rates.(inlayer), RF, 0)); %DEBUG && tphase==1));
               end
               telapsed_input = toc + telapsed_input;
               
               %% For each timestep process static input
               for ti=2:(totaldelay)
                  % If we have a retina, process input via retina first
                  if haveRetina
                     rates.(retina)(end,:) = toVec(generateRetinaInput(sz.(retina), sz.(inlayer), margin, cell_loc.(retina),...
                                                                   Ra.(retina), rates.(inlayer)(end,:), RF));
                  else
                     rates.(inlayer)(end,:)= in(:); % constant input
                  end
                  % Propagate input down each connection
                  for ci=1:size(layerconns,1)
                     pre    = layerconns{ci,1};
                     post   = layerconns{ci,2};
                     cxlabel= [pre post];
%                      pre_sz = size(rates.(pre));
%                      del    = maxdelay.(cxlabel);
                     for i = 1:prod(sz.(post))
                        % To allow recurrent connections we want inputs to be from last time
                        % sample, else if there are recurrent connections & we calculate inputs
                        % using current time point, then when we calculate inputs going the
                        % other direction in the recurrent connections, they might get half
                        % built inputs. So treat last time step as rates that we're still
                        % building, calculating each group's input to the rate one at a time,
                        % and the previous timestep is the current, fully built, rate.
                        inrates = getVisualCortexRateWithDelay(rates.(pre)(1:end-1,:),...
                                                            conns.(cxlabel){i},      ...
                                                            delayinds.(cxlabel){i},  ...
                                                            psp,dt);

                        % Update rates: F_j = Ra + Rb*sum(c_ij*F_i)
                        rates.(post)(end,i) =                                        ...
                                     updateVisualCortexRates(Rb.(cxlabel),           ...
                                                            inrates,                 ...
                                                            weights.(cxlabel){i},    ...
                                                            rates.(post)(end,i)); % need for multiple inputs into layer

                     end % for each post-synaptic neuron
                  end % end for each connection
                  % Rotate rates (keeping recent history for axonal delays)
                  for li=1:nlayers
                     lname = layernames{li};
                     if rectify
                        rates.(lname)(end,rates.(lname)(end,:)<0) = 0;
                     end
                     % Rotate rates (keeping recent history for axonal delays)
                     rates.(lname)(1:end-1,:) = rates.(lname)(2:end,:);
                     rates.(lname)(end,:)     = Ra.(lname); % lambda.(post); % initialise next timestep
                  end
               end % for each timestep
               % Get last rate for each neuron as output for this phase
               for li=1:nlayers
                  lname = layernames{li};
                  if layernames{li}==inlayer
                     % if input layer we need to extract the section within
                     % the margin, which is the actual rates of the layer
                     lrate = reshape(rates.(lname)(end,:),sz.(lname)+margin*2); 
                     outrates.(lname)(tphase,:) = toVec(lrate(margin+1:(end-margin),margin+1:(end-margin)));
                  else
                     outrates.(lname)(tphase,:) = rates.(lname)(end-1,:);
                  end
               end
               
            end % for each phase 
            sprintf('Time taken to generate inputs with %d phases is %.2f\n',...
               length(time),telapsed_input);
            % Record details for this orientation (across all phases)
            if iti==1 % same for each iteration
               if genMovies
                  visual_input{ai,si}.movie = M;
                  clear M;
               else
%                   visual_input{ai,si}.input = in;
               end
               visual_input{ai,si}.angle    = angle;
               visual_input{ai,si}.sfreq    = sfreq;
               visual_input{ai,si}.time     = time;
               visual_input{ai,si}.input_phase = input_phase;               
            end
            visual_weights{ai,si} = weights;
            visual_output{ai,si}  = outrates;
            
            if DEBUG
               rates_debug = outrates.(retina)(:,theta_sorted); 
               response    = zeros(num_angles,1);
               for ang_debug=1:num_angles
                  switch responsefn
                     case 'f1'
                        response(ang_debug) = max(rates_debug(:,ang_debug)) - min(rates_debug(:,ang_debug)); 
                     case 'max'
                        response(ang_debug) = max(rates_debug(:,ang_debug)); 
                  end
               end
               figure(figdebug); hold on; 
               plot(theta_sorted, response)
            end
         end % for each orientation
         
         % Record response details for this iteration (start/end of plasticity)
         tuning.iterations.(iterations{iti}).weights = visual_weights;
         tuning.iterations.(iterations{iti}).output  = visual_output;
         tuning.amplitude     = amplitude; 
         tuning.dc_offset     = dc_offset;
         tuning.grating_angles= grating_angles;
         tuning.input         = visual_input;
         tuning.input_type    = input_type;
         tuning.isplastic     = outstruct.layerconfig.plastic;
         tuning.layernames    = outstruct.layerconfig.layernames;
         tuning.num_freqs     = num_freqs;
         tuning.num_angles    = num_angles;
         tuning.phase_type    = phase_type;
         tuning.spatial_freq  = spatial_freqs;
         tuning.stddev        = stddev;
         tuning.temporal_freq = temporal_phase;
%          tuning.input_phase   = input_phase;
%          tuning.time          = time;

         tstop_sfreq = toc( tstart_sfreq );
         printTime( tstop_sfreq, sprintf( '\tTuning curves for freq %g generated in', sfreq ) );
      end % end for each spatial frequency
      % close figure if spatial grating movies were plotted for first iteration
      if iti==1
         if genMovies, close(fh); end
      end
   end
   
   if DEBUG
      legend(str2legend('\theta=',grating_angles));
   end
   
   tstop = toc(tstart);
   printTime(tstop,'Tuning curves generated in ');
   save('tuningTmp', 'tuning', '-v7.3'); % for variables > 2 GB
   
   % not plotting -> nothing left to do so can exit
   if ~exist('doplot','var') || ~isempty(doplot)
      doplot = true; 
   end
   
   % plotting resulting tuning curves
   if doplot
      if ~exist('plotlayers','var') || isempty(plotlayers)
         plotlayers = getPlasticLayers(tuning);
         if isempty(plotlayers), plotlayers = layernames; end
      end
      try
         curves  = generateVisualCortexTuningCurves( tuning, plotlayers, responsefn, iterations, true, 'orient', false );
         %% Plot sharpness of tuning curves
      catch ME
         cprintf('*errors',sprintf('\nError plotting tuning curves (line %d: %s), returning tuning curve data...\n\n',...
                  ME.stack(1).line, ME.message));
         return;
      end
   end
end

%% Default values
function optargs = parseInputs(varargin, sz, inlayer, lambda)
   % Note: Number of spatial freqs depends on total number of samples, in
   %       the same way that df depends on T. 
   DEBUG          = false;
   % Get info from which parameters can be calculated
   dt             = sqrt(2);       % slowest time evolution is on 45 degree angle
   dc_offset      = lambda.(inlayer);
   amplitude      = dc_offset/2;  % size of input stimulus (I guess in rate (1/s)?)
   freq_max       = 1./dt/sqrt(2); % adjust for orientation of 45 degrees
   sz_input       = sz.(inlayer);
   delta_freq     = 2/min(sz_input);
   freq_min       = delta_freq;
   % kernel size determines number of angles if we have a retina, else use
   % layer size, which will determine max number of possible angles
   if isfield(sz,'kernel')
      sz_kernel   = sz.kernel;
   else
      sz_kernel   = sz.(inlayer);
   end
   
   % Set default parameteter values 
   spatial_freqs  = freq_min:delta_freq:freq_max; % set 1st since others depend on it
   num_freqs      = length(spatial_freqs); % have to populate since user can set it
   num_angles     = max(sz_kernel)*2;
   d_angle        = 180/num_angles;
   
   genMovies      = true;
   grating_angles = 0:d_angle:(180-d_angle);
   input_type     = 'grating';
   iterations     = {'start','end'};
   max_phases     = 1./freq_min;
   min_phases     = 1./freq_max; 
   num_phases     = [];
   num_angles     = length(grating_angles);
   phase_type     = 'variable';
   rectify        = true;
   responsefn     = 'max'; % response fn options: 'f1', 'max'
   stddev         = 2; % std dev if input type is wgn 
   
   if d_angle<1
      cprintf('Warning',sprintf('\tdelta angle is less than 1 (%.2g) for kernel size %d \n',...
                                 d_angle, sz_kernel(1)));
   end
   
   % Old way of calculating & generating spatial frequencies 
%    % Assume that each of the input layer squares is 1 spatial sample -->
%    % sz.input./sz_firstlayer gives dt, or ds in this case since its spatial.
%    % Then half of this is the Nyquist frequency.
%    freq_max       = ratioIn2Ret/2;     % max freq = num pixels per layer cell / 2
%    % fmax = 1/dt, and Tmax = 1/df --> df = 1/Tmax = 1/(sz_retina * ratioIn2retina)
%    delta_freq     = min(1./(sz_input))*4;
%    freq_min       = delta_freq;
%    % max unique angles in retina layer - assumes square
%    num_angles     = min(sz_input);
%    d_angle        = 180/num_angles*4;
%    input_type     = 'grating';
%    grating_angles = 0:d_angle:(180-d_angle);
   optargs = v2struct(amplitude,      dc_offset,  DEBUG,      genMovies, ...
                      grating_angles, input_type, iterations, ...
                      max_phases,     min_phases, num_angles, num_freqs, ...
                      num_phases,     phase_type, rectify,    responsefn, ...
                      spatial_freqs,  stddev);

   if isempty(varargin)
      return;
   elseif ~isint(length(varargin{1})/2)
      error('processVisualCortexInput:inputError - number of optional arguments must be even...exiting...\n');
   end
   varname = @(var) inputname(1); 
   % test inputs which gives the length of a vector - e.g num_freqs
   isvec_sz    = @(var) (isfloat(var) && isint(var) && var>0);
   isvec       = @(var) (isfloat(var) && isvector(var));
   isposscalar = @(var) isfloat(var) && isscalar(var) && var>0;
   isbool      = @(var) (islogical(var));
   isfnname    = @(var) ischar(var);

   for vi=1:2:length(varargin{1})
      value = varargin{1}{vi+1};
      switch varargin{1}{vi}
         case 'amplitude'
            if ~isposscalar(value)
               fprintf('processVisualCortexInput:inputError - amplitude must be a positive scalar, using default value %g \n', amplitude);
            else
               amplitude = value;
            end
            
         case 'dc_offset'
            if ~isposscalar(value)
               fprintf('processVisualCortexInput:inputError - dc_offset must be a positive scalar, using default value %g \n', dc_offset);
            else
               dc_offset = value;
            end
            
         case 'DEBUG'
            if ~(isbool(value) || ischar(value) || isnumeric(value))
               fprintf('processVisualCortexInput:inputError - DEBUG must be a boolean, using default value %d \n', DEBUG);
            else
               if ischar(value) 
                  if any(strcmpi(value,{'true','t'})), DEBUG = true;
                  elseif strcmpi(value,{'false','f'}), DEBUG = false;
                     fprintf('processVisualCortexInput:inputError - DEBUG must be boolean, using default value %d \n', DEBUG);               
                  end
               elseif isbool(value)
                  DEBUG = value;
               elseif isnumeric(value)
                  DEBUG = ternaryOp( value > 0, true, false );
               end
            end
            
         case 'genMovies'
            if ~isbool(value)
               if ischar(value) 
                  if strcmpi(value,'false')
                     genMovies = false;
                  else
                     genMovies = true;
                  end
               elseif isnumeric(value) && isscalar(value)
                  if value==0
                     genMovies = false;
                  else
                     genMovies = true;
                  end
               else
                  fprintf('processVisualCortexInput:inputError - genMovies must be a boolean, using default value %d \n', genMovies);
               end
            else

               genMovies = value;
            end
            
         case {'grating_angles','gratings','angles'}
            if ~isvec(value)
               fprintf('processVisualCortexInput:inputError - %s must be a numerical vector, using default value \n', grating_angles);
            else
               grating_angles = value;
               num_angles = length(grating_angles);
               d_angle    = mean(diff(grating_angles)); % just populate but it's redundant
            end
            
         case 'input_type'
            if ischar(value) 
               if any( strcmpi(value, {'grating','gratings','gauss','wgn','poisson','wpn'}) )
                  input_type = value;
               else
                  error('processVisualCortexInput:inputError - input type must be ''grating'' or ''wgn'' or ''wpn'' \n');
               end
            else
               fprintf('processVisualCortexInput:inputError - genMovies must be a string, using default value %d \n', input_type);
            end
            
         case {'iterations','iteration','iter'}
            if ~iscell(value) && ischar(value), value = {value}; end
            if ~(iscell(value) && all( cellfun( @(it) (ischar(it) && any( contains(it,{'start','end','perc'}) ) ), value) ) )
               error('processVisualCortexInput:inputError - iterations must be a cell array of strings containing ''start'' or ''end'' or ''perc'' \n');
            end
            iterations = value;
            
         case 'max_freq'
            if ~isscalar(value)
               fprintf('processVisualCortexInput:inputError - max_freq must be a scalar number, using default value %d \n',max_freq);
            else
               spatial_freqs(spatial_freqs > value) = [];
               num_freqs = length(spatial_freqs);
            end
            
         case 'max_phases'
            if ~isvec_sz(value)
               fprintf('processVisualCortexInput:inputError - max_phases must be a positive integer, using default value %d \n',max_phases);
            else
               max_phases = value;
            end
            
         case 'min_phases'
            if ~isvec_sz(value)
               fprintf('processVisualCortexInput:inputError - min_phases must be a positive integer, using default value %d \n',min_phases);
            else
               min_phases = value;
            end
            
         case 'num_freqs'
            if ~isvec_sz(value)
               fprintf('processVisualCortexInput:inputError - num_freqs must be a positive integer, using default value %d \n',num_freqs);
            else
               num_freqs     = value;
               delta_freq    = freq_max / num_freqs;
               freq_min      = delta_freq;
               spatial_freqs = freq_min:delta_freq:freq_max;
            end
            
         case 'num_phases'
            if ~isvec_sz(value)
               fprintf('processVisualCortexInput:inputError - num_phases must be a positive integer, using default value %d \n',num_phases);
            else
               num_phases = value;
               phase_type = 'constant';
            end
            
         case 'num_angles'
            if ~isvec_sz(value)
               fprintf('processVisualCortexInput:inputError - num_angles must be a positive integer, using default value %d \n',num_angles);
            else
               num_angles     = value;
               d_angle        = 180/num_angles;
               grating_angles = 0:d_angle:(180-d_angle);
            end
            
         case 'rectify'
            if ~isbool(value)
               error('processVisualCortexInput:inputError - rectify must be a boolean, using default value %d \n',rectify);
            end
            rectify = value;
            
         case 'responsefn'
            if ~( isfnname( value ) || iscell( value ) )
               error('processVisualCortexInput:inputError - responsefn must be string ''max'' or ''f1'', using default value %s',responsefn);
            end
            responsefn = value;
            
         case {'spatial_freq','spatial_freqs','freq','freqs','sfreq','sfreqs'}
            if ~isvec(value)
               fprintf('processVisualCortexInput:inputError - spatial_freq must be a numerical vector, using default value ') 
               fprintf(' %g' ,spatial_freqs); fprintf('\n \n');
            else
               spatial_freqs = value;
               num_freqs     = length(spatial_freqs);
               delta_freq    = mean(diff(spatial_freqs)); % just populate but it's redundant
            end

         case 'stddev'
            if ~isposscalar(value)
               fprintf('processVisualCortexInput:inputError - stddev must be a positive float, using default value %d \n',num_angles);
            else
               stddev        = value;
            end
            
      end
   end
   
   % if input type is noise, we don't need grating angles & frequencies
   if any( strcmpi( input_type, {'wgn', 'wpn'} ) )
      grating_angles = 1;
      max_phases = 1;
      min_phases = 1;
      num_angles = 1; 
      num_freqs  = 1;
      num_phases = 1;
      phase_type = 'constant';
      spatial_freqs = 0;
   end
   
   optargs = v2struct(amplitude,  dc_offset,  DEBUG,      dc_offset,  genMovies,     grating_angles, ...
                      input_type, iterations, max_phases, min_phases, num_angles,    num_freqs, ...
                      num_phases, phase_type, rectify,    responsefn, spatial_freqs, stddev);
   
end











