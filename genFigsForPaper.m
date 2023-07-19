%% visual cortex figures for paper
% FF: indices 9 & 17
%  L: indices 9 & 11

% Convert cycles per pixel axis labels to cycles per degree
% set(gca,'xticklabel', cellfun(@(d) num2str(d), ...
%    cellfun(@(d) d*8, cellfun(@(s) str2num(s), get(gca,'xticklabel'),...
%    'uniformoutput',false) ,'uniformoutput',false),'uniformoutput',false))

% % Convert x-axis from semilogx to linear 
% set(gca,'xscale','linear')

% Convert to cpd by multiplying by 8


%% TO DO:
% Katie - run tuning for L network with natural image input
% Katie - compare network stats to Errol's requested stats
% Katie - run simulation for FF network
% Katie - run tuning for FF network 
% Katie - run tuning for FF network with natural image input
% Errol - generate figures of orientation CV for FF & L networks
% Errol - generate feature map using large layer C 
% Errol - run large dataset for FF network
% Errol - complete the network figure

T = true; F = false;

network = 'L'; % 'FF'; % 

% Simulation data figures
genFreqCurves   = T; % generate spatial frequency data from generateVisualCortexTuningCurves
genOrientCurves = F; % generate orientation data from generateVisualCortexTuningCurves
plotSpatialFreq = F; 
plotOrientation = F; % polar plots of orientation tuning
plotOrientDist  = F; % plot histogram of preferred orientations
plotSchallData  = F;

% Experimental data figures
expSpatialFreq  = F;
expOrientation  = F;

saveFigs        = F;

% can simply write 'all' iterations if you want them all, but then don't
% write over it every time you run because it has to be converted into a
% list of iteration names, which we'll do when we prep our data 
iterations      = 'all'; % {'start', 'perc30', 'end'}; % 'all' %
responseFn      = 'max'; % 'max'; % 'f1'; % 
layers          = {'A','B','C'}; % {'C'}; %  % can choose for spatial freq
layernames      = {'retina', 'LGN', 'V1'};

lineopts   = {'linewidth',3};
fontopts   = {'fontweight','bold','fontsize',10};
if strcmpi(network,'FF') && exist('FFtuning','var')
   sfreq   = FFtuning.spatial_freq;
   plastic = getPlasticLayers( FFtuning );
   if ischar( iterations ) && strcmpi( iterations, 'all' )
      iterations   = fieldnames( FFtuning.iterations );
   end
   % get data for spatial freq, orientation & feature map (for orient dist) all at once
   curves = generateVisualCortexTuningCurves( FFtuning, layers, responseFn, iterations, false, {'orient','freq','feature','3d'}, false );
elseif strcmpi(network,'L') && exist('Ltuning','var')
   sfreq   = Ltuning.spatial_freq;
   plastic = getPlasticLayers( Ltuning );
   if ischar( iterations ) && strcmpi( iterations, 'all' )
      iterations   = fieldnames( Ltuning.iterations );
   end
   % get data for spatial freq, orientation & feature map (for orient dist) all at once
   curves = generateVisualCortexTuningCurves( Ltuning, layers, responseFn, iterations, false, {'orient','freq','feature','3d'}, false );
elseif exist('tuning','var')
   sfreq   = tuning.spatial_freq;
   plastic = getPlasticLayers( tuning );
   if ischar( iterations ) && strcmpi( iterations, 'all' )
      iterations   = fieldnames( tuning.iterations );
   end
   % get data for spatial freq, orientation & feature map (for orient dist) all at once
   curves = generateVisualCortexTuningCurves( tuning, layers, responseFn, iterations, false, {'orient','freq','feature','3d'}, false );
else 
   fprintf('\nNo tuning struct found, exiting...\n');
   return;
end
if strcmpi(network,'FF') && exist('FFstruct','var')
   sz      = FFstruct.layerconfig.sz;
elseif strcmpi(network,'L') && exist('Lstruct','var')
   sz      = Lstruct.layerconfig.sz;
elseif exist('tuning','var')
   sz      = outstruct.layerconfig.sz;
else 
   fprintf('\n outstruct not found, exiting...\n');
   return;
end

%% Plot spatial frequency
if plotSpatialFreq
   maxDisplay = 5; sf = [];
   if genFreqCurves
      if strcmpi(network,'FF') && exist('FFtuning','var')
         curves  = generateVisualCortexTuningCurves(FFtuning,layers,responseFn, ...
                                                    iterations, false, 'freq', false);
         theta   = FFtuning.grating_angles;
      elseif strcmpi(network,'L') && exist('Ltuning','var')
         curves  = generateVisualCortexTuningCurves(Ltuning,layers,responseFn, ...
                                                    iterations, false, 'freq', false);
         theta   = Ltuning.grating_angles;
      else
         curves  = generateVisualCortexTuningCurves(tuning,layers,responseFn, ...
                                                    iterations, false, 'freq', false);
         theta   = tuning.grating_angles;
      end      
   end
   
   fh = zeros(length(layers)+1,1); % +1 for a plot of mean spatial frequency
   for li=1:length(layers)
      fh(li)= figure; hold on;
      layer = layers{li};
      N     = min([maxDisplay sz.(layer)(1)]);
      M     = min([maxDisplay sz.(layer)(2)]);
      if sz.(layer)>maxDisplay
         dN    = sz.(layer)(1) / maxDisplay;
         dM    = sz.(layer)(2) / maxDisplay;
         dNM   = round(dN*dM);
         ninds  = (1:(N*M))*floor(dN*dM);
      else
         ninds = 1:prod(sz.(layer));
      end
      
      % for each neuron in the list, extract spatial frequency data - if
      % the layer is non-plastic, get only the 1st iteration (since they're
      % all the same), and if the layer is plastic, get all iterations
      for ni=1:length(ninds) % for each neuron in the layer
         % plot spatial freq for each iteration
         for ii=1:length(iterations)
            if ii==1 || any( layer==plastic )
               it = iterations{ii};
               y  = curves.(it).(layer).(responseFn)(:,:,ninds(ni));
               % get local max so we know what grating angle to plot the spatial freq of
               [ymax, imax] = max(y(:)); 
               [amax,smax] = ind2sub(size(y),imax);
               sf.(it).(layer)(:,ni) = y(amax,:);
            end
         end
      end
      
      % now plot the difference btwn start & end spatial frequencies if the
      % layer is plastic, else just plot the mean SF for the layer
      if ~any( layer==plastic )
         plot(sfreq, mean( sf.(it).(layer),2 ), lineopts{:} );
         set(gca, fontopts{:});

      else
         % plot plastic layers separately for each neuron
         for ni=1:length(ninds)
            loc = ind2subv(sz.(layer),ninds(ni));  % get neuron subscript
            subplot(N,M,ni), hold on; 
               for ii=1:length(iterations)
                  it = iterations{ii};
                  plot(sfreq, sf.(it).(layer)(:,ni), lineopts{:});
               end
               % give the figure a title
               tstr = sprintf('%s: [%d,%d] \\theta=%.0f',layer,loc(1),loc(2),round(theta(amax)));
               title(tstr);
               % for last neuron add a legend for iterations
               if ni==length(ninds)
                  legend( iterations );
               end
         end
         set(gca, fontopts{:});
         
      end 
      yl = get(gca,'ylim');
      set(gca, 'ylim', [0 yl(end)]);

%       ylim([0 Ymax]); 
      xlim([0 sfreq(end)]);
      xlabel('Spatial frequency (1/pixel)', fontopts{:});
      ylabel('Amplitude', fontopts{:});
      set(fh(li),'name',sprintf('Layer %s',layer));
      yl = get(gca,'ylim');
      set(gca, 'ylim', [0 yl(end)]);
      
      figname = sprintf('Layer%s_%s_spatialfreq', layer, responseFn);
      if saveFigs
         saveFigure(fh(li),figname,1,'fig','pdf');
      end
   end
   
   % plot mean layer spatial frequency results of the first iteration for 
   % non-plastic layers, and all iterations for plastic layers
   fh(li+1) = figure; leg = cell(0);
   for ll=1:length(layers)
      layer = layers{ll};
      for ii=1:length(iterations)
         it = iterations{ii};
         % non-plastic --> only plot 1st iteration (they're all the same)
         if ~any( layer==plastic )
            if ii==1
               plot(sfreq, mean( sf.(it).(layer),2 ), lineopts{:});
               leg{end+1} = layer;
            end
            
         else
            hold on;
            plot(sfreq, mean( sf.(it).(layer),2 ), lineopts{:});
            leg{end+1} = [layer '_{' it '}'];
         end
      end
      
      set(gca, fontopts{:});
   end
   legend( leg );
   xlim([0 sfreq(end)]);
   xlabel('Spatial frequency (1/pixel)', fontopts{:});
   ylabel('Amplitude', fontopts{:});
   set(fh(li+1),'name','All layers');
   yl = get(gca,'ylim');
   set(gca, 'ylim', [0 yl(end)]);
   figname = sprintf('LayeAll_%s_spatialfreq', responseFn);
   if saveFigs
      saveFigure(fh(li+1),figname,1,'fig','pdf');
   end
end

%% Plot orientation
if plotOrientation
   maxDisplay = 5; orient = [];

   if genOrientCurves
      layers  = {'A','C'}; 
      if strcmpi(network,'FF') && exist('FFtuning','var')
         curves  = generateVisualCortexTuningCurves(FFtuning,layers,responseFn, ...
                                                    iterations, false, 'orient', false);
         theta   = FFtuning.grating_angles;
         Theta   = [theta(:); theta(:)+180; theta(1)]; % 360 degree version
      elseif strcmpi(network,'L') && exist('Ltuning','var')
         curves  = generateVisualCortexTuningCurves(Ltuning,layers,responseFn, ...
                                                    iterations, false, 'orient', false);
         theta   = Ltuning.grating_angles;
         Theta   = [theta(:); theta(:)+180; theta(1)]; % 360 degree version
      else
         curves  = generateVisualCortexTuningCurves(tuning,layers,responseFn, ...
                                                    iterations, false, 'orient', false);
         theta   = tuning.grating_angles;
         Theta   = [theta(:); theta(:)+180; theta(1)]; % 360 degree version
      end
   end
   
   % Extract orienation data
   for li=1:length(layers) 
      layer = layers{li};
      N     = min([maxDisplay sz.(layer)(1)]);
      M     = min([maxDisplay sz.(layer)(2)]);
      if sz.(layer)>maxDisplay
         dN    = sz.(layer)(1) / maxDisplay;
         dM    = sz.(layer)(2) / maxDisplay;
         dNM   = round(dN*dM);
         ninds  = (1:(N*M))*floor(dN*dM);
      else
         ninds = 1:prod(sz.(layer));
      end
      % for each neuron in the list, plot start/end/difference of spatial freq
      for ni=1:length(ninds) % for each neuron in the layer
         % select neuron, extract its response score
         loc   = ind2subv(sz.(layer),ninds(ni));  % get neuron subscript
         % plot spatial freq for all iterations - only 1st iteration has
         % data for non-plastic layers since it doesn't change over the
         % simulation
         for ii=1:length(iterations)
            it = iterations{ii};
            if ii==1 || any( layer==plastic )
               it = iterations{ii};
               y = curves.(it).(layer).(responseFn)(:,:,ninds(ni));
               % get local max so we know what grating angle to plot the spatial freq of
               [ymax, imax] = max(y(:)); 
               [amax,smax] = ind2sub(size(y),imax);
% smax = 7; % set all spatial frequency indices to that of the kernel

               orient.(it).(layer)(:,ni) = y(:,smax);
            end
         end
      end
   end
   
   % Plot orientation data
   N = min([maxDisplay sz.C(1)]);
   M = min([maxDisplay sz.C(2)]);
   if sz.C>maxDisplay
      dN    = sz.C(1) / maxDisplay;
      dM    = sz.C(2) / maxDisplay;
      dNM   = round(dN*dM);
      ninds  = (1:(N*M))*floor(dN*dM);
   else
      ninds = 1:prod(sz.C);
   end
   for ii=1:length(iterations)
      it   = iterations{ii};
      fh = figure;
      % now plot retina & v1 orienations together for each iteration
      for ni=1:length(ninds)
         loc   = ind2subv(sz.(layer),ninds(ni));  % get neuron subscript
         % if plastic, get data for each iteration, else just get from 1st
         if ~any( 'A'==plastic )
            y = orient.(iterations{1}).A(:,ni);
         else
            y = orient.(it).A(:,ni);
         end
         z = orient.(it).C(:,ni);
         y = y/max(y)*max(z);
         % plot layer C separately for each neuron
         subplot(N,M,ni),
            polarplot( deg2rad(Theta), [y; y; y(1)] ); hold on;
            polarplot( deg2rad(Theta), [z; z; z(1)] ); hold on;
            % give the figure a title
            tstr = sprintf('%s: [%d,%d] f=%.2f',layer,loc(1),loc(2),sfreq(smax));
            title(tstr);
            if ni==length(ninds)
               legend('retina','v_1');
            end
            set(gca, fontopts{:});
      end
      set(fh,'name',sprintf('Layer %s %s',layer, it));
      figname = sprintf('LayersA_C_%s_%s_orientation', it, responseFn);
      if saveFigs
         saveFigure(fh,figname,1,'fig','pdf');
      end  
   end
   
   % Plot all iterations together
   fc = figure;
   for ni=1:length(ninds)
      loc   = ind2subv(sz.(layer),ninds(ni));  % get neuron subscript
      % plots layer C separately for each neuron
      for ii=1:length(iterations)
         it = iterations{ii};
         y(:,ii) = toVec( orient.(it).C(:,ni) );
      end
      % normalise by dividing by max of iteration, multiply by max of final
      % iteration 
      y = y ./ repmat(max(y),[size(y,1) 1]) * max( y(:,end) );
      subplot(N,M,ni), leg = cell(0);
         for ii=1:length(iterations)
            polarplot(deg2rad(Theta), [y(:,ii); y(:,ii); y(1,ii)]); hold on;
            leg{end+1} = ['v1_{' iterations{ii} '}'];
         end
         % give the figure a title
         tstr = sprintf('%s: [%d,%d] f=%.2f',layer,loc(1),loc(2),sfreq(smax));
         title(tstr);
         % add legend to final neuron 
         if ni==length(ninds)
            legend( leg );
         end
         set(gca, fontopts{:});
   end
   set(fc,'name',sprintf('Layer %s all iterations',layer));
   figname = sprintf('LayerC_start_end_%s_orientation', responseFn);
   if saveFigs
      saveFigure(fc,figname,1,'fig','pdf');
   end

end

%% Plot distribution of orientation
if plotOrientDist
   if genOrientCurves
      if strcmpi(network,'FF') && exist('FFtuning','var')
         curves  = generateVisualCortexTuningCurves(FFtuning, {'A','C'}, responseFn, ...
                                                    iterations, false, 'feature', false);
         theta   = FFtuning.grating_angles;
         Theta   = [theta(:); theta(:)+180; theta(1)]; % 360 degree version
      elseif strcmpi(network,'L') && exist('Ltuning','var')
         curves  = generateVisualCortexTuningCurves(Ltuning, {'A','C'}, responseFn, ...
                                                    iterations, false, 'feature', false);
         theta   = Ltuning.grating_angles;
         Theta   = [theta(:); theta(:)+180; theta(1)]; % 360 degree version
      else
         curves  = generateVisualCortexTuningCurves(tuning, {'A','C'}, responseFn, ...
                                                    iterations, false, 'feature', false);
         theta   = tuning.grating_angles;
         Theta   = [theta(:); theta(:)+180; theta(1)]; % 360 degree version
      end
   end
   Nbins   = 180 / length( theta ); %% TODO: CHECK THIS SHOULDN'T BINS BE THETA ?? %%
   Nbins   = length( theta );
   bins    = theta;
   blims   = [bins(1) - (bins(2) - bins(1))/2, bins(end) + (bins(2) - bins(1))/2];
   
   nA = prod(sz.A); 
   nC = prod(sz.C); 
   
   lw      = 3;
   figtype = @hist; % @hist; % @plot %
   figopts = {'fontweight','bold', 'fontsize',14};
   figure; nr=ceil(sqrt(length(layers))); nc=ceil(length(layers)/nr); fi=1;
   
   for li=1:length( layers )
      layer = layers{li};
      % if the layer is not plastic, just plot orient dist once from 1st iteration
      if ~any( layer==plastic )
         x = curves.(iterations{1}).(layer).fm;
         subplot(nr,nc,fi), fi=fi+1;
            N = hist( x, Nbins );
            N = N / sum(N);  
            bar(bins, N); xlim([0 180]);
            ylabel('Proportion retinal cells', figopts{:});
            xlabel('Preferred orientation', figopts{:});
            xlim( blims );

      % if the layer is not plastic, plot orient dist for all iterations
      else
         for ii=1:length( iterations )
            it = iterations{ii};
            x  = curves.(it).(layer).fm;
            % use the same bins for all iterations
            if ii==1
                N = hist( x, Nbins );
                N = N / sum(N);
            else
               n = hist( x, bins );
               N(:,ii) = n / sum(n);
            end
         end
         subplot(nr,nc,fi), fi=fi+1;
            bar(bins, N); 
            xlim( blims );
            ylabel(['Proportion ' layernames{li} ' cells'], figopts{:});
            xlabel('Preferred orientation', figopts{:});
      end
   end
%    subplot(nr,nc,fi), fi=fi+1;
%       hold on;
%       figtype(yC, sort(xCend(:)), 'linewidth',lw);
%       xlabel('V_1 cell number', figopts{:});
%       ylabel('Preferred orientation', figopts{:});
      lh = legend( iterations );
      set(lh,'box','off');
      set(lh, figopts{:});
      
   ah = findobj(gcf,'type','axes');
   set(ah, figopts{:});
   
   figname = sprintf('LayerC_%s_orientDistribution', responseFn);
   if saveFigs
      saveFigure(gcf,figname,1,'fig','pdf')
   end

end

%% Experimental data figures
ft_size = 16;
l_wid   = 4;
if expSpatialFreq
   % Lateral inhibitory network (from Levin)
   y = 21.*10.^( (-log10(0.04)/43.4).*([ 2.5 17.5 12.6 4.7 45 34 3.5 ]-2.5) + log10(0.04) );

   % x: Spatial Frequency (cpd)
   x = 10.^( (-log10(0.02)/52.33).*[ 17.4 22 31 40 50 59 64 ] + log10(0.02) );

   h = figure;
   semilogx(x, y,'-xk','linewidth',l_wid,'markersize',10)
   hold on
   %legend('vertical grating','horizontal grating');
   ylabel('spikes/s','fontsize',ft_size)
   xlabel('Spatial Frequency (cpd)', 'fontsize',ft_size)
   axis([ 10^(-2) 3 0 25 ])
   set(gca,'fontsize',ft_size)
   hold off
   
   % Feed-forward excitatory network
   x = [5e-2 9e-2 2e-1 1e0 2e+0 2.3e+0];
   y = [0.25 0.4 0.9 1 0.75 0.3];
   figure; 
   semilogx(x, y, '-b.', 'linewidth',l_wid, 'markersize',20);
   xlim([1e-2 3e0]);
   set(gca,'fontsize',ft_size);
end 

if expOrientation
   % Actual orientation bias of a real LGN cell
   theta = -1*( -180:360/36:180 )*2*pi/360;
   rad   = [ 25 26 24 22.5 21.25 18.75 17 14.5 14 13 12.5 12.75 13.5 13 13.5 15.5 19 27 29.75 26.25 23 20.5 18.5 15.5 14.5 14 14.5 15 14.5 14.5 14 15 15.5 17.25 20.5 23 25 ]/29.75;
   % theta_model_array = [ stim_orientations(end) (stim_orientations - 180) stim_orientations ]*2*pi/360;
   % R_model_array = [ L(1).is(16).tuning_curves_bar_ON(end) L(1).is(16).tuning_curves_bar_ON L(1).is(16).tuning_curves_bar_ON ]/max(L(1).is(16).tuning_curves_bar_ON) 
   h     = figure;
   polarplot(theta,rad,'-kx')

   % Lateral inhibitory network
   maxAmp = 20; % spikes/s
   theta  = [0 12.5 25 50 90 ...
            180-25 180-12.5 180 ...
            180+12.5 180+25 270 ...
            360-25 360-12.5 360];
   rad    = [0.5 0.3 0.2 0.05 0 ...
             0.25 0.51 0.75 ...
             0.49 0.22 0 ...
             0.25 0.4 0.5];
   figure; 
   polarplot(deg2rad(theta),rad*maxAmp,'-b.','markersize',20);
   rlim([0 maxAmp]);
   
   % Feed-forward excitatory network
   theta  = [-90 -45 0 45 90];
   rad    = [0.24 0.45 1 0.5 0.23]; 
   theta  = [theta theta(1:end)+180];
   rad    = [rad rad(1:end)];
   figure; 
   polarplot(deg2rad(theta), rad, '-b.', 'markersize',20, 'linewidth',4);
end

if plotSchallData
   % from tracey-mctrace face
   angles90 = 10:10:90;
   angles180= [angles90 angles90+90];
   Ncells   = [17 14 10 6 5 2 6 17 22];
   Nschall  = length(Ncells); % number of bins = 9
   
   % mirror data around 90 degrees
%    angles = [angles angles+90];
%    Ncells = [Ncells Ncells(end:-1:1)];
%    Ncells = Ncells/sum(Ncells);
   
   xAstart = curves.start.A.fm(:);
   xA90    = mod(xAstart, 90);
   xCstart = curves.start.C.fm(:);
   xCend   = curves.end.C.fm(:);
   xCs90   = mod(xCstart, 90);
   xCe90   = mod(xCend, 90);
   Nbins   = 18;
   
   nA = prod(sz.A); 
   nC = prod(sz.C); 
   
   lw      = 3;
   figtype = @hist; % @hist; % @plot %
   figopts = {'fontweight','bold','fontsize',18};
   figure; nr=1; nc=2; fi=1;
   
   subplot(nr,nc,fi), fi=fi+1;
      
      if Nbins==18
         xA = cell2mat(outstruct.layerconfig.RF.angleSet);
         yA = hist(xAstart, xA);
         xl = 180;
         anglesShift = angles180 - 5; % centre rather than edge histogram
      else
         xA = cell2mat(outstruct.layerconfig.RF.angleSet);
         xA(xA>90) = [];
         yA = hist(xA90, angles90);
         xA = angles90;
         xl = 90;
         anglesShift = angles90 - 5; % centre rather than edge histogram
      end
      yA      = yA / nA;
      bar(xA, yA); xlim([0 xl]);
      set(gca,'xtick',0:20:xl);
      ylabel('Proportion retinal cells', figopts{:});
      xlabel('Preferred orientation', figopts{:});

   subplot(nr,nc,fi), fi=fi+1;      
      if Nbins==18
         yCstart = hist(xCstart(:), angles180); 
         yCend   = hist(xCend(:), angles180); 
         xA      = angles180;
         xl      = 180;
         
         Ncells  = [Ncells Ncells(end:-1:1)];
      else
         yCstart = hist(xCs90(:), angles90); 
         yCend   = hist(xCe90(:), angles90); 
         xA      = angles90;
         xl      = 90;
      end
      % normalise
      yCstart = yCstart / sum(yCstart);
      yCend   = yCend / sum(yCend);
      Ncells  = Ncells / sum(Ncells);
      
      lb = bar(xA, [yCstart' yCend']); xlim([0 xl]); 
      hold on;
      lm = plot(angles180, Ncells,'k.','markersize', 20);
      ylabel('Proportion V_1 cells', figopts{:});
      xlabel('Preferred orientation', figopts{:});
      lh = legend([lm lb], 'Experimental data', 'Before plasticity', 'After plasticity');
      set(lh,'box','off');
      set(lh, figopts{:});
      ah = findobj(gcf,'type','axes');
      set(ah, figopts{:});
   
end

% Spatial frequency of cat V1 cells

% x is in log_10 scale, but data traced linearly
x = log10([0.358363326	0.520392809	0.64677579	0.949770906	1.113420585	1.181472959]);
y = [0.256691713	0.421248688	0.971881604	1.03939217	0.733484948	0.305214898];





























