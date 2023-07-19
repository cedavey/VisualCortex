% kernel = retina_kernel(I,doplot)
% Generate a kernel for an on/off receptive field in the retina. 
% Inputs:
%  I  - structure containing standard deviation values for centre and
%  surround, horizontal & vertical. An example struct is given below
%    RF.kernel_size     = 10;           % size of receptive field
%    RF.sd_surr_hor     = 5;            % surround horizontal stnd dev
%    RF.sd_surr_ver     = 3;            % surround vertical stnd dev
%    RF.sd_cent_hor     = 5;            % centre horizontal stnd dev
%    RF.sd_cent_ver     = 3;            % centre vertical stnd dev
%    RF.coeff_cent      = 1.00e-4;      % scale of centre Gaussian 
%    RF.coeff_surr      = 0.30e-4;      % scale of surround Gaussian 
%    RF.ratio_surr_cent = 2.5;          % ratio of surround size compared to centre
%                                         i.e. 2.5:1 in this example
% Outputs:
%  kernel - time domain kernel modelling on/off receptive field
function dog = retina_kernel(RF,doplot)
   if nargin<2
      doplot      = true;
   end
   if nargin<1
      RF.kernel_size  = 10;
      RF.sd_surr_hor  = 4.5; 
      RF.sd_surr_ver  = 4.5;
      RF.sd_cent_hor  = 5;
      RF.sd_cent_ver  = 3; 
      RF.coeff_cent   = 1.20e-4;
      RF.coeff_surr   = 0.50e-4;
      RF.ratio_surr_cent = 1.75;
   end
   ratio_sur_cent = RF.ratio_surr_cent;
   kernel_size    = RF.kernel_size;
   spatial_extent = kernel_size; % normalise stnd dev & scale by this amount

   [xind,yind]    = meshgrid((0:(kernel_size-1))-(kernel_size-1)/2);

   % Determine actual stnd dev from image size & ratio of cent/surr
   sd_surr_hor    = ratio_sur_cent*spatial_extent*RF.sd_surr_hor/(RF.sd_surr_hor + RF.sd_surr_ver)/(1+ratio_sur_cent);
   sd_surr_ver    = ratio_sur_cent*spatial_extent*RF.sd_surr_ver/(RF.sd_surr_hor + RF.sd_surr_ver)/(1+ratio_sur_cent);
   sd_cent_hor    = spatial_extent*RF.sd_cent_hor/(RF.sd_cent_hor + RF.sd_cent_ver)/(1+ratio_sur_cent);
   sd_cent_ver    = spatial_extent*RF.sd_cent_ver/(RF.sd_cent_hor + RF.sd_cent_ver)/(1+ratio_sur_cent);

   dog_centre     = @(x,y) exp(-(x.^2/sd_cent_hor^2 + y.^2/sd_cent_ver^2));
   dog_surround   = @(x,y) exp(-(x.^2/sd_surr_hor^2 + y.^2/sd_surr_ver^2));
   dog_kernel     = @(c,s) c-s; % - (c(1)-s(1));

   dog_cent       = RF.coeff_cent*dog_centre(xind,yind);
   dog_surr       = RF.coeff_surr*dog_surround(xind,yind);
   dog            = dog_kernel(dog_cent,dog_surr);
   scale          = sum(abs(dog(:)));
   
   switch RF.normfn
      case 'mean'
         norm_fn  = @(d) d - mean(d(:));
       case 'norm'
         % linear transformation from old [min max] to new [min max]
         max_old = max(dog(:)); min_old = min(dog(:)); 
         max_new = 1;           min_new = -0.1;
         a = (max_new - min_new)/(max_old - min_old); 
         b = max_new - a*max_old;
         norm_fn  = @(d) a*d + b;
      case 'power'
         %% 2 dimension power spectral density
         % Integral of PSD is not quite equally variance, but not sure why
%          dx  = 1; dy = 1; [Nx, Ny] = size(dog);
%          Kx  = 1/(2*dx); Ky  = 1/(2*dy);    % Nyquist in x & y
%          Nx2 = Nx/2+1;   Ny2 = Ny/2+1;      % half plane in x & y direction
%          X   = Nx*dx;    Y   = Ny*dy;       % x & y ranges
%          dKx = Kx/(Nx/2);dKy = Ky/(Ny/2);   % spacing in x & y freqs
%          kernel   = dog - mean(dog(:));     % 0 mean kernel
%          KERNEL   = fft2(kernel)*dx*dy;     % 2D fourier transform of kernel
%          KERNEL   = KERNEL(:,1:Nx2);        % keep only right half plane
%          PSD      = (2/(X*Y)) * abs(KERNEL).^2; % power spectral density
%          PSD2     = dKx*dKy*sum(PSD(:));    % integral of PSD 
%          power    = var(kernel(:));         % temporal variance
         
         %% 1 dimension process for calculating power spectral density
%          T        = length(dog(:));
%          tpts     = (0:1:T-1)';         % length of Hanning window = length(X)
%          hann     = 1/2*(1-cos(2*pi/1024*tpts)); % Hanning window
%          filt_dog = dog.*hann;          % filtered kernel - apply Hanning window
%          DOG      = fft(filt_dog);      % take fft of X with Hanning window applied
%          NumUnique= ceil((T)/2);        % Calc num unique points 
%          DOG      = DOG(1:NumUnique);   % FFT is symmetric, throw away second half   
%          freq     =(0:NumUnique-1)*Fs/T;% evenly spaced freq vector with NumUniquePts points. 
%         % scale the fft so that it is not a function of the length of x 
%          DOG      = abs(DOG)/length(filt_dog); 
%          PSD      = DOG.^2;             % power spectral density - square of mag of fft
%          % Since we dropped half the FFT, we multiply mx by 2 to keep the 
%          % same energy.  The DC component and Nyquist component, if it 
%          % exists, are unique and should not be mulitplied by 2. 
%          if rem(T, 2) % odd nfft excludes Nyquist point
%            PSD(2:end) = PSD(2:end)*2;
%          else
%            PSD(2:end -1) = PSD(2:end -1)*2;
%          end

         norm_fn  = @(d) d./sqrt(sum(d(:).^2));
      case 'std'
         norm_fn  = @(d) reshape(stnd(d(:)),[kernel_size kernel_size])/max(stnd(d(:))); % d-mean(d(:)); % d/scale; % 
      case 'sum'
         norm_fn  = @(d) d./sum(abs(d(:)));
      otherwise
         fprintf('RF: unknown normalisation function (%s), defaulting to power...\n',RF.normfn);
         norm_fn  = @(d) d/sum(d.^2);
   end
   dog            = norm_fn(dog);

   if doplot
      image       = randn(kernel_size*10)*3;
      IMAGE       = fft2(image);
      DOG         = fft2(image,size(image,1),size(image,2));
      OUTPUT      = IMAGE.*DOG;
      output      = ifft2(OUTPUT);

      figure; i=1;
      ax = subplot(2,2,i); i=i+1;
      cc = fix(kernel_size)/2;
%       plot(yind(:,fix(kernel_size/2)),[dog_cent(:,cc)/scale dog_surr(:,cc)/scale dog(:,cc)]);
      plot(yind(:,fix(kernel_size/2)),[dog_cent(:,cc) dog_surr(:,cc) dog(:,cc)]);
      hold on;
      ax.ColorOrderIndex = 1;
      plot(yind(:,fix(kernel_size/2)),[dog_cent(cc,:)' dog_surr(cc,:)' dog(cc,:)'], '--');
      title('dog');
      legend('centre','surround','DoG');

      subplot(2,2,i), i=i+1;
      imagesc(flipud(dog),[-max(abs(dog(:))) +max(abs(dog(:)))]);
      title('dog');
      cmap = redgreencmap(128);
      colormap(cmap(end:-1:1,:));
      colorbar;
      axis image;

      subplot(2,2,i), i=i+1;
      plot(dog);
      title('dog');
%       imagesc(flipud(dog_cent));
%       title('centre kernel');
%       colorbar;
%       axis image;

      subplot(2,2,i), i=i+1;
      plot([sd_surr_hor sd_surr_ver sd_cent_hor sd_cent_ver RF.coeff_cent RF.coeff_surr],'o');
      xlim([0 7]);
      set(gca,'xtick',1:6,'xticklabel',{'surr_h','surr_v','cent_h','cent_v','A','B'});
      title('DOG - \sigma');

   end
end
