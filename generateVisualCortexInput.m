% rates = generateVisualCortexInput(input_type, lambda, sz_input, margin, params)
% Generate noisy rates using a given distribution. For multiple neurons
% the noise is generated independently for each neuron.
%
% E.g. 
%   rates = generateVisualCortexInput(lambda, sz_input, input_type, params, margin)
% 
%   rates = generateVisualCortexInput(lambda, sz_input, input_type, params)
%
%   rates = generateVisualCortexInput(lambda, sz_input)
%
% Inputs:
%  lambda:     mean neuronal output rate for layer
%  sz_layer:   2D size of input layer to generate values for
%  input_type: 'gaussian' (noise --> assumed mean of zero)
%              'uniform' (noise --> assumed mean of zero)
%              'grating'
%              'poisson' (default)
%  params:     parameters for input type, 
%                gaussian parameter  - {mean, std dev, rectify}
%                  if only mean provided -> gaussian noise, 
%                  if std dev provided -> assume hyperparam & gen rayleigh
%                    std dev (mean of rayleigh = sigma * sqrt(pi/2) with
%                    sigma gen from rayleigh & mean shifted by mu
%                gratings parameters - {t, angle, temp_freq, spatial_freq, amp, dc_offset, rectify}
%                poisson parameter   - lambda (mean & variance, defaults to 10)
%                uniform parameter   - {upper_bound, rectify}
%                  even upper/lower bounds assumed [-upper_bound upper_bound])
%  margin:     size of margin if needed to ensure that the kernel doesn't wrap
%              (defaults to 0)
% Outputs:
%  rates:  resulting rates, where rates = input_type + lambda
%          rates has the size of the input layer, wrapped with an edge the
%          size of each layer's input from the input layer. E.g. for an
%          input layer size of [50 50], & layer size of [10 10], the rates
%          will be [50 50] + an edge of 5 around the whole outside
function rates = generateVisualCortexInput(lambda, sz_input, varargin)
% function rates = generateVisualCortexInput(lambda,sz,sz_ratio,varargin)
   nargin  = length(varargin);
   optargs = {'poisson',10,0};
   optargs(1:nargin) = varargin;
   [input_type, params, margin]  = optargs{:};
   
   switch input_type
      %% Gaussian noise
      case {'gauss','gaussian','norm','normal','wgn'}
         if isscalar(params)
            fn = @(N) randn(N)*params; % params = std dev
            
         elseif length(params)==1
            i  = 1;
            std = params{i}; i=i+1; % mean of gaussian
            rectify = false;
            
            % If hyperparameter is applied to mean
            % A gaussian mean is conjugate prior whereas gaussian variance
            % is not, therefore use the hyperparameter on the mean instead.
            % If params is vector, first parameter is mean of the mean & 
            % second is hyperparameter giving variance of the noise mean
            
            % Got confused - if single param assume mean of noise is 0 &
            % param is std dev of gaussian noise
            fn = @(N) randn(N).*std;
            
         elseif length(params)==2
            i = 1;
            mu      = params{i}; i=i+1; % mean of gaussian
            std     = params{i}; i=i+1; % std dev of gaussian
            rectify = false;

            % sigma   = randn(N,1)*std + mu; % stupid cuz taking abs below
            sigma   = raylrnd(std, sz_input + margin*2) + mu; 
            % If hyperparameter is applied to variance
            % if params is vector, first parameter is mean of std dev of the
            % noise, & second is hyperparameter, giving std dev in std dev
            % of the noise
            fn = @(N) randn(N) .* sigma;
            
         elseif length(params)==3
            i = 1;
            mu      = params{i}; i=i+1; % std dev of gaussian
            std     = params{i}; i=i+1; % mean of gaussian
            rectify = params{i}; i=i+1; % enforce signal being > 0

            sigma   = randn(N)*std + mu;
            % If hyperparameter is applied to variance
            % if params is vector, first parameter is mean of std dev of the
            % noise, & second is hyperparameter, giving std dev in std dev
            % of the noise
            fn = @(N) randn(N).*abs(sigma);
            
         end
         rates = lambda + fn(sz_input + margin*2);
         if rectify, rates(rates<0) = 0; end
         
      %% Sinusoidal grating
      case {'grate','grating','gratings'}
         % create grating stimulus with orientation 'angle'
         if ~isempty(params)
            i = 1;
            t       = params{i}; i=i+1; % timestep
            angle   = params{i}; i=i+1; % angle of grating
            tphase  = params{i}; i=i+1; % temporal frequency
            sfreq   = params{i}; i=i+1; % spatial frequency
            amp     = params{i}; i=i+1; % stimulus amplitude (not incl. background)
            offset  = params{i}; i=i+1; % dc offset for rates sinusoid
            rectify = params{i}; i=i+1; % do not allow negative rates
         else
            margin  = [0 0];
            t       = 1;
            angle   = 45; 
            tphase  = 1; 
            sfreq   = 1; 
            amp     = 30;
            offset  = 50;
            rectify = false;
         end
         % theta = 0 --> grating is constant across the x-axis
         theta  = angle;
         % when creating input, wrap it with an extra cell to avoid edge effects
         [X,Y]  = meshgrid(1:(sz_input(1)+margin(1)*2),1:(sz_input(2)+margin(2)*2));
%          X      = (X-margin(1)) / sz_input(1)*sz_layer(1)*2*pi; % convert freq from 0...num_pixels to 0...num_neurons
%          Y      = (Y-margin(2)) / sz_input(2)*sz_layer(2)*2*pi;
         X      = (X-margin(1));
         Y      = (Y-margin(2));
         ramp   = -X*sind(theta) + Y*cosd(theta);
         rates  = cos(2*pi*sfreq*ramp + tphase*t); % + lambda
         rates(abs(rates)<(1e-8)) = 0;
         if max(abs(rates(:)))==0
            rates = zeros(size(rates)); 
         else
            rates = rates / max(abs(rates(:)));
         end
%          rates  = amp*rates + amp + offset; % make sin wave sit on axis
         rates  = amp*rates + offset;
         
         if rectify, rates(rates<0) = 0; end
         % What Levin did 
%          theta = angle;
%          [X,Y]  = meshgrid(1:sz(1),1:sz(2));
%          ramp   = X*cosd(theta) + Y*sind(theta);
%          rates  = lambda + cos(sfreq*ramp + tfreq*t);
%          figure; imagesc(rates);

     case {'natural','image','images'}
         % noiseparams: [timesteps/img margin img_list curr_img]
         margin = params{2};
         img    = params{4};
         sz     = sz_input + margin*2;
         [hn wn]= size(img);                % height & width of natural image
         ho     = sz_input(1) + margin*2;
         wo     = sz_input(2) + margin*2;   % height & width of required output img
         centre = ceil([hn/2 wn/2]);
         
         % extract output img from input img - can do 1 of several possible
         % things here since we most likely want a much smaller img 
         % - rescale natural img to our size --> gets blurry
         % - just extract a portion of the natural image
         
         % resize - use min ratio of height to width of natural vs output
         rates = imresize(img, max(ho/hn, wo/wn)); 
         rates = rates(1:ho, 1:wo) * lambda; % extract required size & scale rate
         
         % get patch from img centre of required size
%          rates = img(centre(1)-ceil(ho/2)+1:centre(1)+ceil(ho/2), ...
%                      centre(2)-ceil(wo/2)+1:centre(2)+ceil(wo/2));
        
      %% Poisson noise
      case {'poisson','poiss','exponential','exp','wpn'}         
         rates = poissrnd(lambda,sz_input + margin*2);
         
      %% Uniform noise
      case {'u','uni','uniform'}
         % assume noise on [-upper upper]
         % fyi: uniform var = 1/12*(b-a)^2
         i = 1;
         if length(params==1)
            upper   = params; % not cell array if length is 1
            rectify = false;
         elseif length(params==2)
            i = 1;
            upper   = params{i}; i=i+1;
            rectify = params{i}; i=i+1; % timestep
         end
         
         fn    = @(N) rand(N)*upper*2 - upper;
         rates = lambda + fn(sz_input + margin*2);
         if rectify, rates(rates<0) = 0; end
         

   end
end

% % Background on generating a grating:
% We generated an oriented grating by taking the sine of an oriented ramp. 
% The ramp was generated as a linear combination of the 'x' and 'y' matrices, 
% which shows us that the cross sections of an orientated grating along the 
% x or y dimensions are sinusoids modulating at frequencies based on this 
% linear combination. 
% E.g., if the orientation is 30 degrees, the grating is made like this:
% 
% orientation = -30;  %deg (counter-clockwise from vertical)
% sf = 5; % spatial frequency (cycles/deg)
% 
% ramp = cos(orientation*pi/180)*x - sin(orientation*pi/180)*y;
% 
% grating = sin(2*pi*sf*ramp);
% The spatial frequency along the x and y-dimensions are:
% 
% sfx = sf*cos(orientation*pi/180)
% sfy = sf*sin(orientation*pi/180)
% 








