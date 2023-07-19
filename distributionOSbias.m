%  curves = distributionOSbias(tuning, layernames, biastype, responsefn, iterations, doplot, savefigs)
%     or 
%  curves = distributionOSbias(tuning, layernames)
% 
% Determine distribution of orientation bias & the preferred orientation at
% all spatial frequencies, plus the optimal spatial frequency
% 
% Inputs:
%  biastype    - 'levick' or 'cv' (defaults to levick)
%  responsefn  - 'f1' or 'max' (defaults to max)
%  iterations  - 'start', 'end', or {'start','end'} (defaults to {'start','end'})
%  doplot      - true or false (defaults to true)
%  savefigs    - true or false (saved in to working directory, defaults to false)
%
% Outputs:
%  curves.iteration.layer.bias  - degree of orientation bias
%  curves.iteration.layer.ori   - preferred orientation
%  Each of these has format: cell x spatial_freq+1 because it gives the
%  bias & orientation for all spatial frequencies plus the optimal freq

% Note: angle of Levick's bias exactly coincides with angle of CV
function curves = distributionOSbias(tuning, layernames, varargin)
   optargs = {'levick', 'max',  {'start','end'}, true, false};
   nargs   = length(varargin);
   Nbins   = 10;    % number of histogram bins (Levick used 10)
   optargs(1:nargs) = varargin;
   [biastype, responsefn, iterations, doplot, savefigs] = optargs{:};

   curves     = generateVisualCortexTuningCurves(tuning,layernames,responsefn,iterations,0,{'orient','feature'},0);
   anglesDeg  = tuning.grating_angles(:); Na = length(anglesDeg);
   anglesRad   = deg2rad(anglesDeg); 
   freqs      = tuning.spatial_freq(:);   Nf = length(freqs);
% Nf = 10;

   % for each spatial frequency, to calculate strength of bias sum the
   % response at each frequency (as a vector) and divide by summing the
   % response as a scalar, so that the vectors in the numerator cancel each 
   % other out but the denominator is like the sum if they're all pointing the
   % same direction

   if ischar(layernames), layernames = {layernames}; end
   if ischar(iterations), iterations = {iterations}; end
   % determine which layers are postsynaptic & plastic, so have changing RF
   plasticLayers = getPlasticLayers(tuning);
   
   for ii=1:length(iterations)
      it = iterations{ii};

      for li=1:length(layernames)
         layer = layernames{li};
         
         % for 1st iteration or for plastic layers, calculate bias dist.
         if ii==1 || any(plasticLayers==layer)
            Nc    = size(curves.(it).(layer).(responsefn), 3); % angles x freqs x cells

            switch biastype
               case 'cv'
                  % vect_ori_data has format (n cells x n spat_freqs x modes+CV), where the
                  % mode will be 2 since we're folding 0-180 to get 360, & the 3rd is the CV
                  bias = curves.cv.(it).(layer).(responsefn).vect_ori_data(:, :, 3); % cell x freq 

               case 'levick'
               bias = zeros(Nc, Nf+1); % +1 cuz plot for optimal spatial freq also
               psi   = zeros(Nc, Nf+1); % +1 cuz plot for optimal spatial freq also
               
               for fi=1:Nf
                  % response data in curves has format: angles x freq x cell
                  rad  = permute( curves.(it).(layer).(responsefn)(:,fi,:), [1 3 2]); % angles x cell
                  psi  = rad2deg(curves.cv.(it).(layer).(responsefn).vect_ori_data(:, :, 1)); % cell x freq x mode

%                   % calculate in cartesian coordinates
%                   [x, y] = pol2cart(repmat(angleRad, [1 Nc]), rad); %  angles x cells 
%                   num  = abs( x' * ones(Na, 1) + 1i *  y' * ones(Na, 1) ); % cells x 1
%                   den  = rad' * ones(Na, 1);      % cells x 1
%                   bias(:,fi) = num ./den;
                  
                  % calculate in polar coordinates
                   B   = ( rad' * exp(2i*anglesRad) ) ./ sum(rad)';
                   bias(:,fi) = abs(B);
                   psi(:,fi)  = mod(rad2deg(angle(B)/2)+180, 180);

               end
            end

            % Now get bias at optimal spatial frequency of the cell
            for ci=1:Nc
                       y = curves.(it).(layer).(responsefn)(:,:,ci); % angles x freq
               [~, imax] = nanmax(y(:)); % index of max y 
               [~, fmax] = ind2sub(size(y),imax); % get orient & sf index at max

                switch biastype
                   case 'cv'
                      bias(ci,end+1) = cv(ci, fmax);
                      psi(ci,end+1)  = rad2deg(curves.cv.(it).(layer).(responsefn).vect_ori_data(ci, fmax, 1)); % cell x freq x mode
                      
                   % see Levick 1982 for formula for orientation bias
                   % if Rvec  = R exp(j*2*theta) then 
                   %    Bvec  = sum( Rvec) / sum(R) and 
                   %    B     = abs(Bvec) gives bias while 
                   %    theta = angle(Bvec) gives preferred orientation
                   % Q: why the extra factor of 2 ?
                   case 'levick'
                      % note: to test bias vs polar plot use
%                       Theta = [angles(:); angles(:)+180; angles(1)]; 
%                       figure; polar(deg2rad(Theta), [rad; rad; rad(1)]); 
                      rad = y(:,fmax);
%                       % sum in cartesian coordinates
%                       [x, y] = pol2cart(anglesRad, rad);
%                       bias(ci,end) = abs( sum(x + 1i*y) ) / sum(rad);
%                       bias(ci,end) = sqrt( sum(x)^2 + sum(y)^2 ) / sum(rad);

                      % sum in polar coordinates
                      B   = sum( rad .* exp(1i*2*anglesRad) ) / sum(rad);
                      bias(ci,end) = abs(B);
                      psi(ci,end)  = mod( rad2deg( angle(B) / 2 ) + 180, 180);
                      
                end
            end

            curves.(it).(layer).bias = bias; 
            curves.(it).(layer).ori  = psi;
            curves.spatial_freqs     = freqs;

         end % end for 1st iteration or plastic layer
      end % end for each layer
   end % end for each iteration

   if doplot
      fh = []; nr = ceil(sqrt(Nf+1)); nc = ceil((Nf+1) / nr);
      switch biastype
         case 'cv'
            maxbin = 1;
         case 'levick'
            maxbin = 0.5;
      end
      dbins = 0.5/(Nbins-1); 
      bins  = (0 : dbins: 0.5) + dbins/2;
      for ii=1:length(iterations)
         it = iterations{ii};
         for li=1:length(layernames)
            layer = layernames{li};

            % for 1st iteration or for plastic layers, plot bias distribution
            if ii==1 || any(plasticLayers==layer)
               fh(end+1) = figure; ah = zeros(Nf+1,1);
               for fi=1:(Nf+1)
                  ah(fi) = subplot(nr,nc,fi); 
   %                   [y, x] = hist(curves.(it).(layer).bias(:,fi) , bins);
                     [y, x] = hist( curves.(it).(layer).bias(:,fi) , Nbins );
                     y = y/sum(y);
                     bar(x, y);
                     ylim([0 max(y)]); 
                     if fi<=Nf
                        title(sprintf('SF %0.2g', freqs(fi)));
                     else
                        title('Optimal spatial freq');
                     end
                     xlabel( 'Bias' ); 
                     ylabel( 'Probability' );
               end
               set(fh(end), 'name', sprintf('%s %s', it, layer));

               yl = max(toVec(cell2mat(get(ah, 'ylim'))));
               switch biastype
                  case 'cv'
                     set(ah, 'xlim', [0 1]);
                  case 'levick'
                     set(ah, 'xlim', [0 0.5]);
               end
               if savefigs
                  figname = sprintf('Layer%s_%s_%s_%sbias_%dbins', layer, responsefn, it, biastype, Nbins);
                  saveFigure(fh(end),figname,0,{'fig','pdf'});
               end
            end
         end
      end
   end % end if plot figs

   if any(strcmpi('A',layernames))
      % Levick's paper
      Lbins = (0.05:0.05:0.5) - 0.025; % the bin positions are edges not centres
      Lbars = [23 54 55 59 28 20 6 3 1 1];
      Lbars = Lbars / sum(Lbars); 
      Lmean = sum( Lbars .* Lbins ); 
      Lvar  = sum( (Lbins - Lmean).^2 .* Lbars );

      Amean = mean(curves.(it).A.bias(:,fi));
      Avar  = var(curves.(it).A.bias(:,fi));

      fprintf('\n\tLevick mean=%.2g, var=%0.2g\n', Lmean, Lvar);
      fprintf('\n\tUs mean=%.2g, var=%0.2g\n', Amean, Avar);
   end
end
%% Levick's paper
% Lbins = (0.05:0.05:0.5) - 0.025; % the bin positions are edges not centres
% Lbars = [23 54 55 59 28 20 6 3 1 1];
% from isabelle's tracy mc-traceface used on Levick's paper
% 23.03554365	54.20632318	55.23393126	59.25873011	28.1735872	20.20962143	6.251280407	3.168456162	1.284510102








