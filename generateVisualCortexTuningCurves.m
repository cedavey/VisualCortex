%   curves = generateVisualCortexTuningCurves( tuning, tuninglayers, responsefn, plotIterations, doplot, plotType, saveFigs )
%
% E. g. 
% generateVisualCortexTuningCurves(tuning, 'C','f1', 'end', true, 'freq', true)
% generateVisualCortexTuningCurves(tuning, {'B','C'}, {'f1','max'},{'start','perc40','end'}, true, {'orient','freq'})
% 
% For each iteration of processVisualCortexInput for which neural response
% to input has been recorded, extract neural tuning curves to
% orientation.
%
%   curves = generateVisualCortexTuningCurves(tuning,curves,doplot)
% You can also call the function with previously calculated tuning curves,
% just to plot the results.
% 
% Inputs:
%  tuning       - output of processVisualCortexInput, containing responses
%  tuninglayers - layers for which tuning curves are requested (cell for
%                 multiple layers, or char for single layer)
%  responseFn   - name of function to evaluate sinusoidal response at each
%                 spatial freq and angle - options are: 'max', 'f0', 'f1'
%  iterations   - cell array containing which iterations to plot
%                 e.g. {'start', 'perc40', 'end', ...}
%  doplot       - set to true to plot tuning curves
%  plotType     - cell array with list of plot types to include {3d','orientation','freq'}
% 
%  tuning structure contains the following fields:
%       tuning.(iteration_name).input{orientation}.input   = input matrix, with format [M M time];
%       tuning.(iteration_name).input{orientation}.angle   = orientation angle;
%       tuning.(iteration_name).input{orientation}.weights = weights.AB - cell array of postsyn weights
%                                                          = weights.BC
%       tuning.(iteration_name).output{orientation}        = A - matrix of time x neuron rates
%                                                            B
%                                                            C
%       tuning.temporal_freq = temporal_freq;
%       tuning.spatial_freq  = spatial_freq;
%       tuning.grating_angles= grating_angles;
%       tuning.input_phase   = input_phase;
%       tuning.time          = time;
%
% Outputs:
%  curves - for each iteration of presentation of input to the network, and 
%           for each layer for which tuning curves are requested, output,
%           a matrix of neuron x response
%
% See also: processVisualCortexInput
% function curves = generateVisualCortexTuningCurves(tuning,tuninglayers,responsefn,plotIterations,doplot,plotType)
%
%   curves = generateVisualCortexTuningCurves(tuning,tuninglayers,responsefn,plotIterations,doplot,plotType,saveFigs)
%
function curves = generateVisualCortexTuningCurves( tuning, varargin )
   iterations = fieldnames(tuning.iterations);
   layers     = fieldnames(tuning.iterations.(iterations{1}).output{1});
   
   % determine which layers are postsynaptic & plastic, so have changing RF
   plasticLayers = getPlasticLayers(tuning);
      
   % by default plot last iteration, and last layer only
   nargs   = length(varargin);
   optargs = {layers(end), 'f1', iterations{end}, true, {'3d'}, false};
   optargs(1:nargs) = varargin(:);   
   [tuninglayers,responsefn,plotIterations,doplot,plotType,saveFigs] = optargs{:};
   
   curves = [];
   if isempty(tuninglayers)       tuninglayers = layers(end);      end
   if isempty(responsefn),         responsefn  = 'f1';             end
   if isempty(plotIterations), plotIterations  = plasticLayers;    end
   if isempty(doplot),                 doplot  = true;             end
   if isempty(plotType),             plotType  = {'3d'};           end
   % Parse inputs
   if ischar(plotIterations),  plotIterations  = {plotIterations}; end
   if ischar(plotType),              plotType  = {plotType};       end
   
   % if multiple response functions in cell array then call this fn for each
   if iscell(responsefn)
      if length(responsefn)>1
         for ri=1:length(responsefn)
            rfn = responsefn{ri};
            curves.(rfn) = generateVisualCortexTuningCurves( tuning, tuninglayers, rfn, plotIterations, doplot, plotType, saveFigs );
         end
         return;
      else
         responsefn = responsefn{:};
      end
   end
   
   plotFreq        = any(strcmpi(plotType,'freq'));
   plotOrient      = any(strcmpi(plotType,'orientation')) || any(strcmpi(plotType,'orient'));
   plotFeatureMap  = any(strcmpi(plotType,'fm'))          || any(strcmpi(plotType,'feature')) ...
                  || any(strcmpi(plotType,'map'))         || any(strcmpi(plotType,'featuremap'));
   plot3D          = any(strcmpi(plotType,'3d'));
   
   if isstruct(tuninglayers)
      curves       = tuninglayers;
      tuninglayers = fieldnames(curves.(iterations{1}));
   elseif ischar(tuninglayers)
      % if multiple layers in vector, separate into individual layers
      tuninglayers = mat2cell( tuninglayers, ones( numel(tuninglayers), 1 ), 1 );
      tuninglayers = tuninglayers(:)';
   elseif isempty(tuninglayers)
      tuninglayers = layers;
   end
   for li=1:length(layers)
      tmp = sqrt(size(tuning.iterations.(iterations{1}).output{1}.(layers{li}),2));
      sz.(layers{li}) = [tmp tmp];
   end
   orientations    = tuning.grating_angles;
   spatial_freq    = tuning.spatial_freq;
   iterations      = plotIterations;
   if ischar(iterations), iterations = {iterations}; end
   
   % For each weight period that was tested (e.g. start & end of
   % plasticity weight evolution, or percXX), extract a function of its
   % response to the different phases of each orientation presented
   if isempty(curves)
      for ii=1:length(iterations)
         it_name = iterations{ii};
         % if this iteration has been requested to have tuning curves
         if any( strcmpi( it_name, plotIterations ) )
            for li=1:length(tuninglayers)
               layer_name = tuninglayers{li};
               
               % get response for all first iterations, but after that only
               % for plastic postsynaptic layers
               if ii==1 || any(layer_name==plasticLayers)
                  f1 = zeros( length(orientations), length(spatial_freq), prod(sz.(layer_name)) );
      %             f0 = zeros(length(orientations),length(spatial_freq),prod(sz.(layer_name)));

                  % for each spatial frequency
                  for si=1:length(spatial_freq)
                     for ai=1:length(orientations)
                        time = tuning.input{ai,si}.time;
                        outrates = tuning.iterations.(it_name).output{ai,si}.(layer_name);
            %             response(:,ai) = responseFn(outrates,orientations); % cosine fitting
      %                   [f1(ai,si,:),f0(ai,si,:)] = neuronFreqResponse(outrates,tuning.temporal_freq,tuning.time); % fourier fitting
                        [f1(ai,si,:)] = getNeuronResponse(outrates,tuning.temporal_freq,time,responsefn);
                     end
                     curves.(it_name).(layer_name).(responsefn) = f1;
                     if any(isnan(f1(:)))
                        fprintf('\n%g%% of f1 in %s of layer %s at theta %d and sf %.2g is NaN\n', ...
                           sum(isnan(f1(:)))/numel(f1(:))*100, it_name, layer_name, round(orientations(ai)), spatial_freq(si));
                     end
      %                curves.(it_name).(layer_name).f0 = f0;
      %                curves.(it_name).(layer_name) = cosinefit(response,orientations);
                  end
               end
            end
         end
      end
   end
   % get best angle for each cell in tuning layer
   if plotFeatureMap
      ori = ori_tuning_CV(tuning, curves, 2, 0); 
      curves.cv = ori;
      for ii=1:length(iterations)
         it_name = iterations{ii};
         
         for li=1:length(tuninglayers)
            layer_name = tuninglayers{li};
            
            if ii==1 || any(layer_name==plasticLayers)              
               fm = zeros(sz.(layer_name));
               
               for ni=1:prod(sz.(layer_name))
                  % select neuron, extract its f1 score
                  y = curves.(it_name).(layer_name).(responsefn)(:,:,ni);
                  [ymax, imax] = nanmax(y(:)); % index of max y 
                  [~,smax]     = ind2sub(size(y),imax); % get freq ind of max
                  fm(ni) = ori.(it_name).(layer_name).(responsefn).vect_ori_data(ni,smax,1)*180/pi;
               end
               curves.(it_name).(layer_name).fm = fm;
            end
         end
      end
   end
   
   % Plot neurons & their curves determined by responseFn - separate
   % neurons vertically & horizontally for readability
   
   if doplot
      warning off
      % Plot layer's response for selection of cells in figure
      for ii=1:length(plotIterations)
         % do a separate figure for each iteration
         it_name = plotIterations{ii};

         % do a separate figure for each layer
         for li=1:length(tuninglayers)
            layer_name = tuninglayers{li};
            % use same cell indices for each iteration, even if randomised
            if ii==1
               maxDisplay = 5;
               N = min([maxDisplay sz.(layer_name)(1)]);
               M = min([maxDisplay sz.(layer_name)(2)]);
               if sz.(layer_name)>maxDisplay
                  dN    = sz.(layer_name)(1) / maxDisplay;
                  dM    = sz.(layer_name)(2) / maxDisplay;
                  dNM   = round(dN*dM);
                  
                  % select cells evenly across the lattice to plot
                  ninds.(layer_name) = (1:(N*M))*floor(dN*dM);
                  
                  % randomise cells that are plotted
%                   ninds.(layer_name) = sort(randperm(prod(sz.(layer_name)), maxDisplay^2));

                  % get first maxDisplay^2 cells
                  ninds.(layer_name) = 1:(maxDisplay^2); 
               else
                  ninds.(layer_name) = 1:prod(sz.(layer_name));
               end
            end
            % if not a plastic postsyn layer we don't care which iteration
            if ~any(layer_name==plasticLayers), itstr=[]; else, itstr = it_name; end
            
            % plot response for all 1st iterations, but after that only
            % for plastic postsynaptic layers
            if ii==1 || any(layer_name==plasticLayers)
               if plot3D
                  fs = figure;             
                  set(fs,'name',[responsefn ' ' itstr ' 3D tuning - ' layer_name]);
                  set(fs,'visible','off');
               end            
               if plotFreq 
                  % figure for 2D tuning curves so can plot multiple together
                  Ymax     = 0; % want to remember common max for common yaxis
                  fh_freq  = figure; hold on;
                  leg_freq = cell(0); % legend for 2D plot of different layers together
                  set(fh_freq,'name',[responsefn ' ' itstr ' frequency tuning - ' layer_name]);
                  set(fh_freq,'visible','off'); 
               end
               if plotOrient 
                  % figure for 2D tuning curves so can plot multiple together
                  fh_orient  = figure; 
                  Ymax       = 0; % want to remember common max for common yaxis
                  leg_orient = cell(0); % legend for 2D plot of different layers together
                  set(fh_orient,'name',[responsefn ' ' itstr ' orientation tuning - ' layer_name]);
                  set(fh_orient,'visible','off'); 
               end
               if plotFeatureMap
                  fh_fm = figure;
                  set(fh_fm,'name',[responsefn ' ' itstr ' feature map - ' layer_name]);
                  set(fh_fm,'visible','off'); 
               end

   % xx=allCombos(4,4);  + 25;
   % ninds = subv2ind(sz.(layer_name),xx);
               for ni=1:length(ninds.(layer_name)) % for each neuron in the layer
                  % select neuron, extract its f1 score
                  loc   = ind2subv(sz.(layer_name),ninds.(layer_name)(ni));  % get neuron subscript
                  loc   = bsxfun(@minus,loc,sz.(layer_name)/2); % centre cell location 
                  y     = curves.(it_name).(layer_name).(responsefn)(:,:,ninds.(layer_name)(ni));
                  if any(isnan(y))
                     fprintf('\n%d%% of y is Nan for nID %d\n',sum(isnan(y(:)))/numel(y(:))*100,ninds.(layer_name)(ni));
                  else
                     % if we're plotting more than 1 iteration we want to make
                     % sure that the amplitudes of each iteration are set to the
                     % same to enable comparison. Assume any percentage
                     % iteration will be somewhere between start & end
                     if any(layer_name==plasticLayers) && length(plotIterations)>1
                        nIt       = length(plotIterations); 
                        other_ind = 1:nIt; other_ind(ii) = [];
                        y_other   = cell(nIt-1,1);
                        for jj=1:length(other_ind)
                           y_other{jj} = curves.(plotIterations{jj}).(layer_name).(responsefn)(:,:,ninds.(layer_name)(ni));
                        end
%                         switch it_name
%                            case {'start'}
%                               y_other = curves.end.(layer_name).(responsefn)(:,:,ninds.(layer_name)(ni));
%                            case {'end'}
%                               y_other = curves.start.(layer_name).(responsefn)(:,:,ninds.(layer_name)(ni));
%                            otherwise
%                         end
                     else
                        y_other = cell(0);
                     end
                     cmap  = max( [y(:); max( cellfun(@(yo) max(yo(:)), y_other) ) ] );
                     cmap  = ternaryOp(cmap==0, 1, cmap);
                     angle = ternaryOp(loc(:,1)~=0, atand(loc(:,2)./loc(:,1)), 90);

      %                   figure(fh);
      %                   subplot(N,M,ni)
      %                   plot(orientations,y);
      %                   title(sprintf('Loc [%.1g %.1g], \\theta=%d',loc(1),loc(2),round(angle)));
      %                   xlabel('Angle'); ylabel('Response');
                     sfreq = tuning.spatial_freq;
                     theta = tuning.grating_angles;
                     Theta = [theta(:); theta(:)+180; theta(1)]; % 360 degree version

                     if plot3D
                        figure(fs);
                        subplot(N,M,ni)
         %                   imagesc(flipud(y'),[0 cmap]); colormap gray;
                        if cmap>0
                           imagesc(y,[0 cmap]); colormap gray;
                        elseif cmap<0
                           imagesc(y,[cmap 0]); colormap gray;
                        else
                           imagesc(y); colormap gray;
                        end
                        title(sprintf('Loc %s, \\theta=%d',num2str(loc),round(angle)));
                        ylabel('\theta'); xlabel('f');
         %                   set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1]);
                        numlabels = min([6 length(theta) length(sfreq)]); 
                        dlabel = length(sfreq)/numlabels; 
                        freqind = floor(1:dlabel:length(sfreq));
                        dlabel = floor(length(theta)/numlabels); 
                        thetaind = 1:dlabel:length(theta);
                        set(gca,'xtick',freqind,'xticklabels',sfreq(freqind));
                        set(gca,'ytick',thetaind,'yticklabels',theta(thetaind));                  
      %                   set(gca,'xtick',1:length(sfreq),'xticklabels',sfreq);
      %                   set(gca,'ytick',1:length(theta),'yticklabels',theta);
                        colorbar;
                     end
                     if plotFreq || plotOrient
                        [ymax, imax] = nanmax(y(:)); % index of max y 
                        Ymax = ternaryOp(ymax>Ymax, ymax, Ymax);
                     end

                     if plotFreq 
                        % y = [angle, spatial_freq, neuron]
                        % y = curves.(it_name).(layer_name).f1(:,:,ninds(ni));
                        figure(fh_freq);
                        [amax,smax] = ind2sub(size(y),imax);
                        subplot(N,M,ni),                  
                        y2d = y(amax,:);
                        if ymax>0, ylim([0 ymax]); end
                        semilogx(sfreq, y2d); xlabel('Spatial freq');
                        lstr = sprintf( '%s: ID=%d \\theta=%.0f', layer_name, ninds.(layer_name)(ni), round(theta(amax)) );
                        title( sprintf( '%d (\\theta=%.2g)', ninds.(layer_name)(ni), round(theta(amax)) ) );
                        leg_freq(end+1) = {lstr};
                     end

                     if plotOrient 
                        % y = [angle, spatial_freq, neuron]
                        % y = curves.(it_name).(layer_name).f1(:,:,ninds(ni));
                        figure(fh_orient);
                        [amax,smax] = ind2sub(size(y),imax); % get orient & sf index at max
                        subplot(N,M,ni),
                        y2d = y(:,smax);
                        polarplot(deg2rad(Theta), [y2d; y2d; y2d(1)]); hold on;
                        lstr = sprintf( '%s: ID=%d f=%.2f', layer_name, ninds.(layer_name)(ni), spatial_freq(smax) );
                        leg_orient(end+1) = {lstr};
                        title( sprintf( '%d (f=%.2g)', ninds.(layer_name)(ni), spatial_freq(smax) ) );
                     end
                  end
               end

               if plotFeatureMap
                  figure(fh_fm);
                  fm = curves.(it_name).(layer_name).fm;
                  pm = phasewrap(fm,'degrees');
                  imagesc(pm); 
                  axis image; % colorbar;
                  phasemap; 
   %                phasebar('deg','location','northeastoutside');
               end
               % Add legends, set axes limits, & save figures if requested
               % Put figs in arrays so don't have to call for each fig type
               % with separate code
               figlabels = cell(0); figleg = cell(0);
               % I think the new version of matlab doesn't allow creating
               % an empty matrix & then filling it with figure handles,
               % while maintaining the handle type, but rather it stores
               % doubles instead, causing errors later (so init with gojects)
               figs = gobjects(0); 
               if plot3D
                  figs(end+1) = fs; 
                  figlabels(end+1) = {'3D'};     
                  figleg(end+1) = {[]};
               end
               if plotFreq
                  figs(end+1) = fh_freq; 
                  ah = findobj(fh_freq,'type','axes');
%                   set(ah,'ylim',[0 Ymax]);
                  figlabels(end+1) = {'freq'};
                  figleg(end+1) = {leg_freq}; 
               end
               if plotOrient
                  figs(end+1) = fh_orient; 
                  figlabels(end+1) = {'orient'}; 
                  figleg(end+1) = {[]};       
               end
               if plotFeatureMap
                  figs(end+1) = fh_fm;
                  figlabels(end+1) = {'feature_map'};
                  figleg(end+1) = {[]};
               end
               for ff=1:length(figs)
                  figure(figs(ff));
                  set(figs(ff),'visible','on');
                  if ~isempty(figleg{ff}), legend(figleg{ff}); end
                  if saveFigs
                     if isempty(itstr)
                        figname = sprintf('Layer%s_%s_%s',upper(layer_name),responsefn,figlabels{ff});
                     else
                        figname = sprintf('Layer%s_%s_%s_%s',upper(layer_name),responsefn,itstr,figlabels{ff});
                     end
                     set(figs(ff), 'filename', [figname '.fig']);
                     saveFigure(figs(ff),[],1,{'fig','pdf'})
                     close(figs(ff)); 
                  end
               end
            end % end for 1st iteration or plastic postsynaptic layer
         end % end for each layer
      end % end for each iteration    
   end % end if doplot
end % end function


function varargout = getNeuronResponse(outrates,temporal_freq,time,fn)
   switch lower(fn)
      case 'f0'
         [f1,f0] = neuronFreqResponse(outrates,temporal_freq,time);
         varargout{1} = f1./f0;
         varargout{2} = f1;
      case 'f1'
         [f1,f0] = neuronFreqResponse(outrates,temporal_freq,time);
         varargout{1} = f1;
         varargout{2} = f0;
      case 'max'
         m = max(abs(outrates));
         varargout{1} = m;
%       case 'cos'
%          r = responseFn(outrates,orientations);
%          varargout{1} = r;
   end        
end


% Response gives some function of rate response across each orientation's
% phases. This function then fits a cosine to each neurons response across
% the orientations. 
% Inputs:
%  response     - response to each orientation (orientation x neuron)
%  orientations - a vector of all orientations presented to each neuron
% function [amp,offset] = responseFn(response,orientations)
%    % allocate output memory for cosine params, ac=eoss
%    type    = fittype('a*cos(d*x+c)+b');
%    nnfit   = fit(toVec(orientations),toVec(response),type);
%    amp     = nnfit.a;
%    centre  = nnfit.c;
%    offset  = nnfit.b;
%    scale   = nnfit.d;
% end


% Response gives some function of rate response across each orientation's
% phases. This function then fits a cosine to each neurons response across
% the orientations. 
% Inputs:
%  response     - response to each orientation (orientation x neuron)
%  orientations - a vector of all orientations presented to each neuron
function cosfit = cosinefit(response,orientations)
   nn   = size(response,1); 
   % allocate output memory for cosine params, ac=eoss
   amp  = zeros(nn,1); centre = amp; offset = amp; scale = amp;
   type = fittype('a*cos(d*x+c)+b');
   for i=1:nn
      nnfit      = fit(toVec(orientations),toVec(response(i,:)),type);
      amp(i)     = nnfit.a;
      centre(i)  = nnfit.c/length(orientations)*orientations(end);
      offset(i)  = nnfit.b;
      scale(i)   = nnfit.d;
   end
   cosfit.amp    = amp; 
   cosfit.centre = centre;
   cosfit.offset = offset;
   cosfit.scale  = scale;
end










