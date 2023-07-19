doplot = true;
% Note: lower case = temporal, upper case = spectral
lambdaA        = 20;
image_size     = [100 100];

image          = randn(image_size)*3;

spatial_extent = 10;
size_sur_cent  = 2.5;
kernel_size    = 10;

[xind,yind]    = meshgrid((0:kernel_size)-kernel_size/2);

I.sd_surr_hor  = 5; 
I.sd_surr_ver  = 3;
I.sd_cent_hor  = 5;
I.sd_cent_ver  = 3; 
I.coeff_cent   = 1.00e-4;
I.coeff_surr   = 0.30e-4;

% Determine actual stnd dev from image size & ratio of cent/surr
sd_surr_hor    = size_sur_cent*spatial_extent*I.sd_surr_hor/(I.sd_surr_hor + I.sd_surr_ver)/(1+size_sur_cent);
sd_surr_ver    = size_sur_cent*spatial_extent*I.sd_surr_ver/(I.sd_surr_hor + I.sd_surr_ver)/(1+size_sur_cent);
sd_cent_hor    = spatial_extent*I.sd_cent_hor/(I.sd_cent_hor + I.sd_cent_ver)/(1+size_sur_cent);
sd_cent_ver    = spatial_extent*I.sd_cent_ver/(I.sd_cent_hor + I.sd_cent_ver)/(1+size_sur_cent);

dog_centre     = @(x,y) exp(-(x.^2/sd_cent_hor^2 + y.^2/sd_cent_ver^2));
dog_surround   = @(x,y) exp(-(x.^2/sd_surr_hor^2 + y.^2/sd_surr_ver^2));
dog_kernel     = @(c,s) c-s; % - (c(1)-s(1));

dog_cent       = I.coeff_cent*dog_centre(xind,yind);
dog_surr       = I.coeff_surr*dog_surround(xind,yind);
dog            = dog_kernel(dog_cent,dog_surr);
IMAGE          = fft2(image);
DOG            = fft2(image,image_size(1),image_size(2));
OUTPUT         = IMAGE.*DOG;
output         = ifft2(OUTPUT);

if doplot
   figure; i=1;
   subplot(2,2,i), i=i+1;
   plot(yind(:,fix(kernel_size/2)),[dog_cent(:,fix(kernel_size/2)) dog_surr(:,fix(kernel_size/2)) dog(:,fix(kernel_size/2))]);
   title('dog');
   legend('centre','surround','DoG');
   colorbar;

   subplot(2,2,i), i=i+1;
   imagesc(dog,[-max(dog(:)) +max(dog(:))]);
   title('dog');
   cmap = redgreencmap(128);
   colormap(cmap(end:-1:1,:));
   colorbar;
   axis image;

   subplot(2,2,i), i=i+1;
   imagesc(output,[-max(dog(:)) +max(dog(:))]);
   title('output');
   colorbar;
   axis image;
   
   subplot(2,2,i), i=i+1;
   plot([sd_surr_hor sd_surr_ver sd_cent_hor sd_cent_ver I.coeff_cent I.coeff_surr]);
   xlim([0 7]);
   set(gca,'xtick',1:6,'xticklabel',{'surr_h','surr_v','cent_h','cent_v','A','B'});
   title('DOG');

end

