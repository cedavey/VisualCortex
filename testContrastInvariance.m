% Generate tuning curves for a range of sinusoidal amplitudes & see if
% orientation tuning curves exhibit contrast invariance (i.e. tuning width
% on the x-axis does not change with stimuli strength, which impacts firing
% rate on the y-axis)

% Thoughts:
%  - changing amp gives contrast invariance, but not changing dc
%  - changing amp with constant dc can change preferred orientation

T = true; F = false;
loaddata    = F;
runtuning   = F;
savetuning  = F;
plotfigs    = T;
savefigs    = F;
spatialfreqs= [0.1];

clear filename; % clear so that if we want default the var doesn't exist
filename    = 'Amp';
network     = 'FF'; % 'FF' 'L'
freq        = 'opt'; % 'opt' for optimal freq, or any other spatial freq in tuning
% freq        = 0.15; 
plottype    = 'plot'; % 'polar' % 'plot' % 

if runtuning
   scale       = 'log'; % 'linear', 'log'
   numAmps     = 4;   % make odd number coz putting default in centre
   deltaAmp    = 0.3; % proportion change for each jump in amplitude
   numDC       = 1;   % make odd number coz putting default in centre
   deltaDC     = 0.2; % proportion change for each jump in dc
   responsefn  = 'max';
   layer       = 'C';
   
   if ~exist('filename','var')
      if numAmps==1
         filename = 'DC';
      elseif numDC==1
         filename = 'Amp';
      else 
         filename = 'AmpDC';
      end
   end
end

if loaddata
   if strcmpi(network,'L')
      ddir = '/Users/catherinedavey/Documents/code/visual_cortex/results/lateral/contrastInvariance';
   elseif strcmpi(network,'FF')
      ddir = '/Users/catherinedavey/Documents/code/visual_cortex/results/feedforward/contrastInvariance';
   else
      str  = sprintf('\nNetwork type must be ''L'' or ''FF'', not %s, exiting...\n\n', network);
      cprintf('Errors', str);
   end
   fname   = 'outstruct.mat';

   load(fullfile(ddir, fname));
   if ~exist('outstruct','var')
      if exist('Lstruct','var')
         outstruct = Lstruct; 
         clear Lstruct;
      elseif exist('FFstruct','var')
         outstruct = FFstruct;
         clear FFstruct;
      else
         cprintf('Error',sprintf('\nCan''t find outstruct, Lstruct or FFstruct in mat file, exiting...\n'));
         return;
      end
   end
end
% tuning stimuli amplitude defaults to 
inlayer           = outstruct.layerconfig.inlayer;
default_lambda    = outstruct.layerconfig.lambda.(inlayer);
default_dc_offset = default_lambda;
default_amplitude = default_dc_offset/2;  % size of input stimulus (I guess in rate (1/s)?)

if strcmpi(scale, 'linear')
   % put default amplitude in the middle, & create vector of amplitudes
   dcvec   = ( (0:(numDC  -1) ) - (numDC  -1)/2 ) * default_dc_offset*deltaDC  + default_dc_offset;
   ampvec  = ( (0:(numAmps-1) ) - (numAmps-1)/2 ) * default_amplitude*deltaAmp + default_amplitude;
else
   dcvec   = 2.^(0:(numDC  -1)) * default_dc_offset*deltaDC  + default_dc_offset;
   ampvec  = 2.^(0:(numAmps-1)) * default_amplitude*deltaAmp + default_amplitude;
end

fprintf('\nDC: ');  fprintf('\t%d', round(dcvec )); fprintf('\n');
fprintf('\nAmp: '); fprintf('\t%d', round(ampvec)); fprintf('\n');

% put default amplitude at the start, & create vecctor of amplitudes
% ampvec  = (0:(numAmps-1) ) * default_amplitude*deltaDC + default_amplitude;
% dcvec   = (0:(numDC-1) ) * default_dc_offset*deltaDC + default_dc_offset;

Ndc     = length(dcvec);
Namp    = length(ampvec);

if any(ampvec)<0
   cprintf('errors', sprintf('testContractInvariance: amplitude has negative components, exiting...'));
   disp(ampvec);
   return;
end
if any(dcvec)<0
   cprintf('errors', sprintf('testContractInvariance: dc_offset has negative components, exiting...'));
   disp(dcvec);
   return;
end

if runtuning
   tuning = cell(Namp,Ndc);
   curves = cell(Namp,Ndc);
   config = {'iterations','end', 'rectify',true, 'max_phases',8, 'sfreq',spatialfreqs, 'num_angles',10};
   for ai=1:length(ampvec)
      amp = ampvec(ai);
      if ai==1
         config{end+1}    = 'amplitude';
         config{end+1}    = amp;
         ampindex         = length(config);
      else
         config{ampindex} = amp;
      end
      
      for di=1:Ndc
         % if only a single dc, set dc t0 0
         if Ndc==1
            dc_offset = 0;
         else
            dc_offset = dcvec(di);
         end
         cprintf('Errors',sprintf('\tComputing tuning curves for amp %d and dc offset %d\n',round(amp), dc_offset));
         if di==1
            config{end+1}   = 'dc_offset';
            config{end+1}   = dc_offset;
            dcindex         = length(config);
         else
            config{dcindex} = dc_offset;
            config{ampindex} = amp;
         end
         
         outname = sprintf('tuning_%dcoffset_%damp', dc_offset, amp);
         tuning{ai,di} = processVisualCortexInput(outstruct, false, [], config);
         
         % fyi: curves.end.C.max = [20 x 15 x 25]; % angles x freq x cell
         curves{ai,di} = generateVisualCortexTuningCurves(tuning{ai,di}, layer, responsefn, {'end'}, false, {'orient'}, false);
         curves{ai,di}.amplitude = amp;
         curves{ai,di}.dc_offset = dc_offset;

      end
   end
   
   if savetuning
      save(fullfile(ddir,[filename '_contrastInvar.mat']),...
         'ampvec', 'curves', 'dcvec', 'deltaAmp', 'deltaDC', 'layer', ...
         'numAmps', 'numDC', 'outstruct', 'responsefn', 'tuning', '-v7.3');
   end
end

if plotfigs   
   sz     = outstruct.layerconfig.sz.(layer);
   angles = tuning{1}.grating_angles;
   freqs  = tuning{1}.spatial_freq;
   Nf     = length(freqs);
   Na     = length(angles);
   Nc     = prod(sz);
   R      = zeros(Na,Nc,Namp,Ndc); % record response for each 
   Theta  = [angles(:); angles(:)+180; angles(1)]; % 360 degree version
   angcentre = ceil(Na/2); % need for centring tuning curve on pref orientation
   
   fh = figure; maxDisplay = 5; 
   N  = min([maxDisplay sz(1)]);
   M  = min([maxDisplay sz(2)]);
   
   if any(sz > maxDisplay)
      dN  = sz(1) / maxDisplay;
      dM  = sz(2) / maxDisplay;
      dNM = round(dN*dM);

      % select cells evenly across the lattice to plot
      ninds = (1:(N*M))*floor(dN*dM);

      % randomise cells that are plotted
%       ninds = sort(randperm(prod(sz.(layer_name)), maxDisplay^2));
   else
      ninds = 1:prod(sz); N = sz(1); M = sz(2);
   end
   
   cv = cell(Namp, Ndc); % circ var to centre figures on preferred orient.
   lh = zeros(Namp, Ndc, Nc); % line handles for all lines in all subplots
   for ai=1:Namp
      for di=1:Ndc
         y   = curves{ai,di}.end.(layer).(responsefn);
         tmp = generateVisualCortexTuningCurves(tuning{ai,di}, layer, responsefn, {'end'}, false, {'feature'}, false);
         tmp = tmp.cv.end.(layer).(responsefn).vect_ori_data;
         
         if ischar(freq) 
            Y = reshape(y, [Na*Nf Nc]);    % collapse freq & orient to find all time max response
            [ymax, imax] = nanmax(Y,[],1); % index of max response at any freq/angle combo
            [amax, fmax] = ind2sub([Na Nf],imax); % get angle & freq index of each cell's max
            % get response at optimal freq: angles x cell
            r = cell2mat(arrayfun(@(i) y(:,fmax(i),i), 1:25, 'uniformoutput',0)); 
            R(:,:,ai,di) = r;
            
            cv{ai,di} = cell2mat(arrayfun(@(i) tmp(i,fmax(i),1), 1:25, 'uniformoutput',0));
            
         else
            [~,fi] = min(abs(freq - freqs));
            r = permute(y(:,fi,:), [1 3 2]); % get response at optimal freq: angles x cell
            R(:,:,ai,di) = r;

            cv{ai,di} = tmp(:,fi,1);
            
         end
         cv{ai,di} = cv{ai,di}/pi*180;
         
         for ni=1:length(ninds)
            ah = subplot(N,M,ni); hold on;
            switch plottype
               case 'plot'
                  % centre plot on cell's preferred orientation
                  % get index of angle the cell's preferred orientation is closest to
                  prefori = cv{ai,di}(ni);
                  rate    = r(:,ni);
                  [~,pp]  = min(abs(angles - prefori)); % preferred index
                  % put preferred orientation in the centre
                  if pp >= ceil(Na/2)
                     si  = pp - angcentre + 1;
                  else
                     si  = pp + angcentre ;
                  end
                  rate   = [rate(si:end); rate(1:si-1)];
                  ang    = [angles(si:end) angles(1:si-1)];
                  lh(ai,di,ni) = plot(1:Na, rate); xlim([1 Na]);
                  if ai==1 && di==1
                     set(ah, 'xtick',angcentre,'xticklabel',num2str(round(angles(pp))));
%                      xt = get(ah, 'xtick');
%                      set(ah, 'xticklabel', num2str(ang(xt)'));
                     title(sprintf('\\theta=%d', round(prefori)));
                  end
                  
               case 'polar'
                  lh(ni) = polar(deg2rad(Theta), [r(:,ni); r(:,ni); r(1,ni)]); hold on;
            end
         end
      end % end for each dc offset
   end % end for each amp
   
   % set unique colours so can disambiguiate lines properly - need lh to be
   % unique for ai, di also!
   cols  = getColourMatrix( Namp*Ndc );
   lh    = reshape(lh, [Namp*Ndc Nc]);
   for ni=1:length(ninds)
      ah = subplot(N,M,ni);
      cc = mat2cell(cols,ones(Namp*Ndc,1),3);
      try
         set( arrayfun(@(i) set(lh(i,ninds(ni)),'color',cols(i,:)), 1:(Namp*Ndc) ,'uniformoutput',false) );
      catch ME
      end
   end
   
   legDC  = str2legend('dc=',dcvec);
   legAmp = str2legend('amp=',ampvec,',  ');
   [di,ai]= allCombos([Namp, Ndc]);
   leg    = [legAmp(ai,:) legDC(di,:)];
   legend(leg);
   if ischar(freq)
      fstr = freq;
   else
      fstr = sprintf('%.2g',freq);
   end
   fstr(fstr=='.')='_';
%    fname  = sprintf('layer%s_contrastInvar_freq%s', layer,ternaryOp(ischar(freq), freq, fstr));
      fname  = sprintf('%s_freq%s', filename,ternaryOp(ischar(freq), freq, fstr));
   set(gcf,'name',fname);
   
   if savefigs
      if ~exist('ddir','var')
         layers = outstruct.layerconfig.layernames;
         if any(strcmpi(layers,'L'))
            ddir = '/Users/catherinedavey/Documents/code/visual_cortex/results/lateral/contrastInvariance';
         elseif strcmpi(network,'FF')
            ddir = '/Users/catherinedavey/Documents/code/visual_cortex/results/feedforward/contrastInvariance';
         else
            str  = sprintf('\nNetwork type must be ''L'' or ''FF'', not %s, exiting...\n\n', network);
            cprintf('Errors', str);
         end
      end
     set(gcf,'filename', fullfile(ddir, [fname '.fig']));
     saveFigure(gcf,[],0,{'fig','pdf'});
   end
end % end plotfigs






