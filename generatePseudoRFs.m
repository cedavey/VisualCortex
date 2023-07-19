%% Calculate pseudo-RFs, or structural RFs
%
%  Sum input kernels weighted by B->C weight to generate pseduo, or
%  structural, receptive fields of layer C cells (it's approx. having an 
%  input grating that is DC)
%
%      rfs = generatePseudoRFs(outfile)
%      rfs = generatePseudoRFs(outstruct)
% OR 
%      rfs = generatePseudoRFs(outfile, doplot, iterations, plotlayers, savefigs)
%
% Use generatePseudoRFs to process constant input for a network simulated 
% in visual_cortex, to generate pseudo, or structural, receptive fields for
% cells deeper in the network hierarchy
% You MUST run visual_cortex prior to this function, and save the output in
% outfile, which is then used by this function, or assign the output to
% outstruct, which can then be input instead.
%
% Inputs:
%  outfile    - name of output file containing output struct from
%               visual_cortex plasticity simulation
%  doplot     - true or false
%  iterations - which iteration (start/end) to generate structure RFs for
%               ( defaults to {'start', 'end'} )
%  plotlayers - cell array of layers to plot the tuning curves of 
%               e.g. {'B', 'C', 'L'} (defaults to all layers)
%  savefigs   - saves figs in all of eps/pdf/fig formats with prefix pseudoRFs
%
% Output:
%  tuning structure that gives responses to input grating of different
%  spatial frequencies and orientations, which are determined by sampling
%  rate and number of samples, as well as input size. Responses are
%  calculated prior to plasticity (iteration = 'start'), and after
%  plasticity has been run (iteration = 'end')
%          tuning.iterations.start/end.weights % weights for iteration
%          tuning.iterations.start/end.output  % output rates for each layer
%          tuning.time                         % time vector;
% 
% You can also request to plot the pseudo RFs for specified layers
% 
% See also: processVisualCortexInput   visual_cortex   
%
function rfs = generatePseudoRFs( outfile, varargin ) % doplot, iterations, plotlayers, savefigs
%% TO DO
   % Ratio of horizontal to vertical from 0 to 2 (1, 1.5, 2)
   % Save retinal convolution computations to save time 
   % Wavers law - double spatial frequencies each round
   % Rate should be ~100 or so when active
   % You can convert to temporal RFs by outputting RF at each timestep
     
   % Are phases appropriate for all spatial frequencies? Currently using
   % the same phases for all --> complete 1 cycle --> always uses index 2
   % of the fft of response, regardless of spatial frequency
   optargs = {true, 'end', [], false};
   nargs   = length(varargin); 
   optargs(1:nargs) = varargin(:);
   [doplot, iterations, plotlayers, savefigs] = optargs{:};
   
   DEBUG = true;
   
   if ischar(outfile)
      load(outfile);
   elseif isstruct(outfile)
      outstruct = outfile; %
   end      

   try
      outstruct2variables;
   catch ME
      str = sprintf('\ngeneratePseudoRFs: cannot open file or struct to get network and layer configs (%s), exiting...\n\n', ME.message);
      cprintf('Errors', '%s',str);
      return;
   end
   if isempty(iterations)
      iterations = {'start','end'};
   elseif ischar(iterations)
      iterations = {iterations};
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

   % Only use this function if you have a retina!!
   haveRetina  = ternaryOp(isempty(retina),false,true);
   if ~haveRetina, return; end
   
   sz_input    = sz.(inlayer);
   sz_retina   = sz.(retina);
   sz.kernel   = RF.kernel_size^2;
   ratioIn2Ret = sz_input./sz_retina; % ratio of num inputs to retinal cells
   margin      = ceil(ratioIn2Ret);  % number of input pixels per retinal cell
   dc_offset   = lambda.(inlayer);
   nL          = length(layernames);
   
   for li=1:nL
      layer = layernames{li};
      rates.(layer) = zeros(prod(sz.kernel), prod(sz.(layer)));
   end

   for iti=1:length(iterations)

      % Get weights to use for RFs
      for ci=1:length(layerconns)
         pre     = layerconns{ci,1};
         post    = layerconns{ci,2};
         cxlabel = [pre post];
         if plastic.(cxlabel)
            switch iterations{iti}
               case {'start'}       % get initial weights
                  for i=1:prod(sz.(post))
                     weights.(cxlabel){i} = outweights.(cxlabel){i}(1,:)';
                  end
               case {'end','final'} % get final weights
                  for i=1:prod(sz.(post))
                     weights.(cxlabel){i} = outweights.(cxlabel){i}(end,:)';
                  end
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
                     for i=1:prod(sz.(post))
                        weights.(cxlabel){i} = outweights.(cxlabel){i}(tpt,:)';
                     end
                  end
            end
         end
      end
               
      % initialise rates for each layer - initialise samples to mean
      % layer rate, except final timepoint, which should be set to Ra
      % for postsynaptic layers, since input from previous layer will
      % + Ra gives mean layer rate for timepoint (previous timepoints
      % will be written over as we advance through processing input)
      rates.(inlayer) = zeros(1, prod(sz.(inlayer)+margin*2));

      % Show static, oriented input for entire trial - trial length is
      % equal to the max delay, so we have 1 full cycle
      rates.(inlayer)(1,:) = dc_offset; % constant input
      rates.(inlayer)(1,:) = 1; % constant input
      % Prepare first timestep
      [~, retina_input] = generateRetinaInput(sz.(retina), sz.(inlayer), margin, cell_loc.(retina),...
                                              Ra.(retina), rates.(inlayer), RF, false);
      rates.(retina) = cell2mat(cellfun(@(c) c(:), retina_input(:)','uniformoutput',false));

      %% For each timestep process static input
      % Process input via retina first, then propagate input down each connection
      includeFiringRate = true;
      for ci=1:size(layerconns,1)
         pre    = layerconns{ci,1};
         post   = layerconns{ci,2};
         cxlabel= [pre post];
         for i = 1:prod(sz.(post))
            inrate_indices = conns.(cxlabel){i};
            R = rates.(pre)(:,inrate_indices);
            W = weights.(cxlabel){i};
            b = Rb.(cxlabel);
            b = 1;
            
            % DEBUG
            if DEBUG 
               shape   = @(x) reshape(x,sz.kernel);
               if strcmpi(cxlabel, 'BC')
                  index   = find(W>0.3); 
                  lgn     = outstruct.ntwkconfig.conns.BC{1}(index); 
                  rgc     = cell2mat(outstruct.ntwkconfig.conns.AB(lgn)); 
                  theta{i}= outstruct.layerconfig.RF.theta(rgc); 
                  kernels = outstruct.layerconfig.RF.kernels(rgc); 
                  rf{i}   = zeros(10,10); 
                  for j=1:length(kernels)
                     rf{i} = rf{i} + kernels{j}; 
                  end
               end
            end
            
            if ~includeFiringRate
               rates.(post)(:,i) = rates.(post)(:,i) + R * W * b;
            else
               l = lambda.(pre); % mean layer firing rate
               a = Ra.(post); % spontaneous firing rate
               % scale R so that its max for each presyn neuron equals the
               % presyn layer's mean firing rate
               R = bsxfun(@times, R, l./max(R));
               R = bsxfun(@rdivide, R, max(R)/l);

%                if strcmpi(cxlabel,'BL'), b = b/3; end
               rates.(post)(:,i) = rates.(post)(:,i) + R * W * b;
               if DEBUG
                  rfs.DEBUG.( iterations{iti} ).(cxlabel)(:,i) = R * W * b;
               end
            end
            
         end % for each post-synaptic neuron
      end % end for each connection

      % Record response details for this iteration (start/end of plasticity)
      rates = rmfield(rates,inlayer);
      rfs.(iterations{iti})  = rates;
      if DEBUG
         rfs.DEBUG.(iterations{iti}).rgc = rf; 
         rfs.DEBUG.(iterations{iti}).theta = theta;
      end
%       rfs.(iterations{iti}).weights = weights;
   end % for each iterations

   
   if isempty(doplot),         doplot = true;        end
   if isempty(savefigs),     savefigs = false;       end
   if isempty(plotlayers), plotlayers = layernames;  end
   if ischar(plotlayers),  plotlayers = {plotlayers};end
   if doplot
      try
         plast = getPlasticLayers( outstruct );
         plotPseudoRFs(rfs, [RF.kernel_size RF.kernel_size], sz, plotlayers, iterations, plast, savefigs);
      catch ME
         cprintf('*errors',sprintf('\nError plotting RFs (%s, line %d), returning RF data...\n\n',...
                  ME.message, ME.stack(1).line));
         return;
      end
   end
end


function plotPseudoRFs(rfs, sz_kernel, sz, plotlayers, iterations, plastic, savefigs)
   warning off
   fh   = zeros(length(plotlayers)*length(iterations),1);
   cmap = redgreencmap(256);
  
   % get min & max across iterations
   for pp=1:length(plotlayers)
       lims.min.(plotlayers{pp}) = Inf; lims.max.(plotlayers{pp}) = -Inf; 
       for ii=1:length(iterations)
          m = min(rfs.(iterations{ii}).(plotlayers{pp})(:));
          if m < lims.min.(plotlayers{pp}), lims.min.(plotlayers{pp}) = m; end 
          m = max(rfs.(iterations{ii}).(plotlayers{pp})(:));
          if m > lims.max.(plotlayers{pp}), lims.max.(plotlayers{pp}) = m; end 
       end
   end
   % Plot layer in figure
   ff = 0;
   for ii=1:length(iterations)
      % do a separate figure for each iteration
      it_name = iterations{ii};

      % do a separate figure for each layer
      for li=1:length(plotlayers)
         layer = plotlayers{li};
         if ~isfield( rfs.(iterations{1}), layer )
            str = sprintf('generatePseudoRFs:inputError: %s layer is not known, ignoring\n',layer);
            cprintf('Errors',str);
         else
            if ii==1 || ( ii>1 && strcmp( layer, plastic ) )
               ff = ff+1; % figure index
               fh(ff) = figure;
               set(fh(ff),'name',[it_name ' pseudo RFs - ' layer]);
               set(fh(ff),'visible','off'); 

               maxDisplay = 5;
               N = min([maxDisplay sz.(layer)(1)]);
               M = min([maxDisplay sz.(layer)(2)]);
               if sz.(layer)>maxDisplay
                  dN    = sz.(layer)(1) / maxDisplay;
                  dM    = sz.(layer)(2) / maxDisplay;
                  ninds  = (1:(N*M))*floor(dN*dM);
      %                ninds = dNM:dNM:prod(sz.(layer_name));
               else
                  ninds = 1:prod(sz.(layer));
               end
      % xx=allCombos(4,4);  + 25;
      % ninds = subv2ind(sz.(layer_name),xx);
               for ni=1:length(ninds) % for each neuron in the layer
                  % select neuron, extract its f1 score
                  loc = ind2subv(sz.(layer),ninds(ni));  % get neuron subscript
      %             loc   = bsxfun(@minus,loc,sz.(layer_name)/2); % centre cell location 
                  y   = rfs.(it_name).(layer)(:,ninds(ni));

                  % y = [angle, spatial_freq, neuron]
                  % y = curves.(it_name).(layer_name).f1(:,:,ninds(ni));
                  figure(fh(ff));
                  subplot(N,M,ni),
                     low   = lims.min.(layer);
                     up    = lims.max.(layer);
                     bound = ternaryOp( abs(low) < abs(lims.max.(layer) ), abs(lims.max.(layer)), abs(low) );
                     imagesc(reshape(y, sz_kernel),[-bound bound]);
   %                   imagesc(reshape(y, sz_kernel));
                     cmap  = redgreencmap(512);
                     colormap(cmap(end:-1:1,:));
         %             colorbar;
                     axis image;
                     str = sprintf('%s: [%d,%d]',layer,loc(1),loc(2));
                     title(str);
                     set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1]);
               end
               colorbar;

               set(fh(ff),'visible','on');
               figname = sprintf('PseudoRF_layer%s_%s',upper(layer),it_name);
               set(fh(ff), 'name', figname);
               if savefigs
                  set(fh(ff),'units','normalized','outerposition',[0 0 1 1])
         %                saveas(figs(ff), figname, 'jpg');
         %                saveas(figs(ff), figname, 'fig');
                  saveFigure(fh(ff), figname, true, {'pdf',fig'});
               end
            end % end if layer is plastic
         end % end if field in rfs
      end % end for each layer
   end % end for each iteration
   
end








