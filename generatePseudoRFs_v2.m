%% Calculate pseudo-RFs, or structural RFs
%
%     rfs = generatePseudoRFs_v2(outfile, doplot, iterations, plotlayers, savefigs)
%
%  Sum input kernels weighted by B->C weight to generate pseduo, or
%  structural, receptive fields of layer C cells (it's approx. having an 
%  input grating that is DC).
%
%  In this version we keep track of where in the input field an input comes
%  from, so we take a weighted sum of each input pixel. The idea is to keep
%  track of the spatial extent of each RF relative to the input field.
%
%      rfs = generatePseudoRFs_v2(outfile)
%      rfs = generatePseudoRFs_v2(outstruct)
% OR 
%      rfs = generatePseudoRFs_v2_v2(outfile, doplot, iterations, plotlayers, savefigs)
%
% Use generatePseudoRFs_v2 to process constant input for a network simulated 
% in visual_cortex, to generate pseudo, or structural, receptive fields for
% cells deeper in the network hierarchy
%
% You MUST run visual_cortex prior to this function, and save the output in
% outfile, which is then used by this function, or assign the output to
% outstruct, which can then be input instead.
%
% Inputs:
%  outfile    - name of output file containing output struct from
%               visual_cortex plasticity simulation, or outstruct itself
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
%          tuning.iterations.start/end.output  % output fields for each layer
%          tuning.time                         % time vector;
% 
% You can also request to plot the pseudo RFs for specified layers
% 
% See also: processVisualCortexInput   visual_cortex   
%
function rfs = generatePseudoRFs_v2( outfile, varargin ) % doplot, iterations, plotlayers, savefigs
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
   vis_input   = sz.(inlayer)+margin*2; % size of visual input, including margin
   dc_offset   = lambda.(inlayer);
   nL          = length(layernames);
   
   for iti=1:length(iterations)
      % Each cell has a representation of the input layer as a spatially  
      % sensitive representation of the summed & weighted inputs - gotta 
      % reset here so removal of margin doesn't muck up next iteration
      for li=1:nL
         layer = layernames{li};
         fields.(layer) = zeros( [ sz.(inlayer)+margin*2 prod(sz.(layer)) ] );
      end

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
      
      % To calculate pseudo RF we're not inputting a sinusoidal grating,
      % but rather inputting each retina cell's kernel, so no phase is reqd
               
      % initialise rates for each layer - initialise samples to mean
      % layer rate, except final timepoint, which should be set to Ra
      % for postsynaptic layers, since input from previous layer will
      % + Ra gives mean layer rate for timepoint (previous timepoints
      % will be written over as we advance through processing input)
      rates.(inlayer)  = zeros(1, prod(sz.(inlayer)+margin*2));

      % Show static, oriented input for entire trial - trial length is
      % equal to the max delay, so we have 1 full cycle
      in = generateVisualCortexInput( dc_offset, sz.(inlayer), 'wgn', dc_offset, margin );
      rates.(inlayer)(1,:) = in(:);
      
%       rates.(inlayer)(1,:) = dc_offset; % constant input
%       rates.(inlayer)(1,:) = 1; % constant input
      % Prepare first timestep
      [~, retinput, visfield] = generateRetinaInput( sz.(retina), sz.(inlayer),    margin, cell_loc.(retina), ...
                                                     Ra.(retina), rates.(inlayer), RF,     false );
      % this is a copy of the input convolved with the retinal cell kernels, so not cell specific
      fields.(retina) = visfield; 

      %% For each timestep process static input
      
      % Process input via retina first, then propagate input down each connection
      includeFiringRate = true;
      
      for ci=1:size(layerconns,1)
         pre     = layerconns{ci,1};
         post    = layerconns{ci,2};
         cxlabel = [pre post];
         
         for i = 1:prod(sz.(post))
            inrate_indices = conns.(cxlabel){i};
            F   = fields.(pre)(:,:,inrate_indices);
            % make field 2D, apply matrix multiplication, then reshape back
            Fsz = size( F );
            F   = reshape( F, [ Fsz(1)*Fsz(2) Fsz(3) ] ); 
            W   = weights.(cxlabel){i};
            b   = Rb.(cxlabel);
%             b = 1;
            
            % DEBUG
            if DEBUG 
               shape = @(x) reshape(x,sz.kernel);
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
               fields.(post)(:,i) = fields.(post)(:,i) + F * W * b;
            else
               l = lambda.(pre); % mean layer firing rate
               a = Ra.(post); % spontaneous firing rate
               % scale R so that its max for each presyn neuron equals the
               % presyn layer's mean firing rate
%                F = bsxfun(@times, F, l./max(F));
%                F = bsxfun(@rdivide, F, max(F)/l);

%                if strcmpi(cxlabel,'BL'), b = b/3; end
               fields.(post)(:,:,i) = fields.(post)(:,:,i) + reshape( F * W * b, vis_input );
               if DEBUG
                  rfs.DEBUG.( iterations{iti} ).(cxlabel)(:,i) = F * W * b;
               end
            end
            
         end % for each post-synaptic neuron
      end % end for each connection

      for li=1:nL
         layer = layernames{li};
         % now we need to ditch the margin around each neuron
         fields.(layer) = fields.(layer)(1:sz_input(1), 1:sz_input(2), :);
      end 
      % Record response details for this iteration (start/end of plasticity)
      rfs.(iterations{iti})  = fields;
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
         % set x & y limits by summing radii across layers 
         radius = sum(  struct2array( outstruct.layerconfig.r ) ) * 3; 
         plotPseudoRFs(rfs, sz, plotlayers, iterations, plast, radius, margin, savefigs);
      catch ME
         cprintf('*errors',sprintf('\nError plotting RFs (%s, line %d), returning RF data...\n\n',...
                  ME.message, ME.stack(1).line));
         return;
      end
   end
end


function plotPseudoRFs(rfs, sz, plotlayers, iterations, plastic, radius, margin, savefigs)
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
         if ~exist( 'fieldsz', 'var' ), fieldsz = size( rfs.(it_name).(layer)(:,:,1) ); end
         if ~isfield( rfs.(iterations{1}), layer )
            str = sprintf('generatePseudoRFs:inputError: %s layer is not known, ignoring\n',layer);
            cprintf('Errors',str);
         else
            if ii==1 || ( ii>1 && any( layer == plastic ) )
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
               
               for ni=1:length(ninds) % for each neuron in the layer
                  % select neuron, extract its f1 score
                  y   = rfs.(it_name).(layer)(:,:,ninds(ni));
                  loc = ind2subv( sz.(layer), ninds(ni) );  % get neuron subscript
                  [~, locRel] = getRelativePosition( fieldsz, ninds(ni), sz.(layer) );
                  locRel = locRel + margin;
                  
                  figure(fh(ff));
                  subplot(N,M,ni),
                     low   = lims.min.(layer);
                     up    = lims.max.(layer);
                     bound = ternaryOp( abs(low) < abs(lims.max.(layer) ), abs(lims.max.(layer)), abs(low) );
                     imagesc( y, [-bound bound] );
                     cmap  = redgreencmap(512);
                     colormap(cmap(end:-1:1,:));
                     
                     % comment out these two lines to zero in on RF
                     ylim( [ max( locRel(1) - radius, 0 ), min( locRel(1) + radius, fieldsz(1)+margin(1) ) ] );
                     xlim( [ max( locRel(2) - radius, 0 ), min( locRel(2) + radius, fieldsz(2)+margin(2) ) ] );
                     
                     % comment out this lines to include whole vis field
                     % axis image;
                     
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
                  saveFigure(fh(ff), figname, true, {'pdf','fig'});
               end
            end % end if layer is plastic
         end % end if field in rfs
      end % end for each layer
   end % end for each iteration
   
end

% Get input image patch for given retinal cell position. Image patch
% assumed to be centred around the cell 
% Inputs:
%  in_sz    - size of visual input
%  kern_sz  - size of receptive field kernel
%  loc      - location of retinal cell given wrt input size
% Outputs:
%  img      - input image patch for retinal cell
%  coords   - visual input corner coords: [ lower left, lower right, upper left, upper right ]
function [img, coords] = getCellInput( input, kern_sz, loc )
   in_sz = size(input);
   upper = @(sz)  ceil((sz-1)/2); % dealing with even kernel sizes --> 
   lower = @(sz) floor((sz-1)/2); % upper & lower halves have unequal sizes
   x     = loc(1);                % x position of cell
   y     = loc(2);                % y position of cell
   sX    = in_sz(1);              % input size in X direction
   sY    = in_sz(2);              % input size in Y direction
   kXu   = upper(kern_sz(1));     % kernel size in upper X direction
   kXl   = lower(kern_sz(1));     % kernel size in lower X direction
   kYu   = upper(kern_sz(2));     % kernel size in upper Y direction
   kYl   = lower(kern_sz(2));     % kernel size in lower Y direction
   
   % if kernel sits within image simply return the image patch
   if 1<=(x-kXl) && (x+kXu)<=sX && 1<=(y-kYl) && (y+kYu)<=sY
      img = input(x-kXl:x+kXu, y-kYl:y+kYu);
      coords = [x-kXl, x+kXu, y-kYl, y+kYu];
      return;
   end
   % if an edge goes over/under then need to wrap around
   ind = allCombos(kern_sz);
   ind = bsxfun(@minus, ind, [kXl kYl]); 
   ind = bsxfun(@plus,  ind, [  x   y]); 
   ind = bsxfun(@mod, ind-1, in_sz) + 1; 
   img = reshape(input(subv2ind(in_sz,ind)),kern_sz);
   % get the 4 edges of the rectangle to determine corner coords (for debugging)
   coords = [ind(1,1), ind(end,1), ind(1,2), ind(end,2)];
   
end





