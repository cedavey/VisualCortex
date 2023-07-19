%   [ output, retina, ret_input ] = generateRetinaInput( retina_sz, input_sz, margin, cell_loc, ...
%                                                        lambda,    input,    RF,     negRF,    DEBUG)
% Generate visual cortex input from input & receptive field kernels. 
% For each retinal cell convolve the input using the given receptive field
%
% Inputs:
%  ret_sz   - 2D size of retina layer (format: [sz_x, sz_y])
%  in_sz    - size of input to the retina (format: [sz_x, sz_y])
%  margin   - size of margin around full layer to avoid wrapping causing artifacts
%  cell_loc - location of retinal cells (don't have to be a grid)
%  lambda   - background/spontaneous rate of retinal cells
%  input    - visual input to the retina
%  RF       - receptive field structure containing the following fields:
%  negRF    - off centre RF so have to multiply by -1
%    RF.kernel_size  = 10; (example values)
%    RF.sd_surr_hor  = 5; 
%    RF.sd_surr_ver  = 3;
%    RF.sd_cent_hor  = 5;
%    RF.sd_cent_ver  = 3; 
%    RF.coeff_cent   = 1.00e-4;
%    RF.coeff_surr   = 0.30e-4;
%    RF.ratio_surr_cent = 2.5;
%    RF.retina_kernel= retina_kernel(RF,false); % receptive field kernel
%    RF.input_var    = 5;
%    RF.kernels      = cell array of kernel indices, with format [x y]
%    RF.theta        = float array of angle values
%
% Outputs:
%  output    - output of the retina, with a single value for each cell, having
%              size as specified by ret_sz
%  retina    - a cell array containing each retinal cell's patch of input * kernel 
%  vis_layer - visual field for each retinal cell, i.e. a copy of the total
%              visual field input for each RGC, with non-zero values only
%              where its kernel is
function [ output, retina, visfield ] = generateRetinaInput( ret_sz, in_sz, margin, cell_loc, Ra, input, RF, negRF, DEBUG )
   if nargin<9 || isempty(DEBUG)
      DEBUG = false;
   end
   if nargin<8 || isempty(negRF)
      negRF = false;
   end
   % TO DO: 
   % - generate WGN & convolve with RF kernels for each retina cell
   N_ret   = prod(ret_sz);                            % number of retinal cells
   N_input = length( toVec( input ) );
   kern_sz = ones(1,2)*RF.kernel_size;                % size of receptive field kernel
   input   = flipud( reshape(input, in_sz+margin*2) );% convert from vector to 2D image
   retina  = cell( ret_sz );
   output  = zeros( ret_sz );
%    output  = zeros(prod(ret_sz),1);

   retina_fn = @(rf,input) conv2( rf, input, 'same' ) / N_input;
   retina_fn = @(rf,input) rf .* input;

   if DEBUG
      fh=figure; title('Input');
      f2=figure; title('2D kernel');
      fr=figure; title('Kernel');
      fc=figure; title('Input convolved with kernel');
      fp=figure; title('Sum of convolution');
      inputmin = min(input(:)); 
      inputmax = max(input(:));
   end
   
   % entire visual field for each RGC, with non-zero values only where it has its kernel
   % (to calc pseudo-rf when input kernels potentially overlap)
   visfield = zeros( [ in_sz + margin*2 N_ret ] ); 
   for ni=1:N_ret
      % calculate retina cell location with respect to input size (i.e.
      % centre input over retina & stretch to equivalent size)
      [indVec,retRel] = getRelativePosition( in_sz, ni, ret_sz );
      retRel      = periodicCellPosition(retRel + cell_loc(ni,:), in_sz); % add jitter to retinal cell location wrt input
      retRel      = retRel + margin;
      % input layer had a margin around the edges to avoid artifacts - add in
      [img,coords]= getCellInput(input, kern_sz, retRel); % get input to this cell using periodic boundaries
      if negRF
         rfi =-RF.kernels{ni}; % receptive field kernel for this cell
      else
         rfi = RF.kernels{ni}; % negative of the RF kernel for this cell
      end
      img_conv    = retina_fn( rfi, img );
      % images are too small for freq domain convolution - too much error
      retina{ni}  = img_conv;
      output(ni)  = sum(img_conv(:)) + Ra; % (lambda was added in genVisCortexInput)
      
      % to calculate pseudo-rf it's useful to show the actual kernel values
      % across the input layer: [ lower left, lower right, upper left, upper right ]
      ni_ind      = repmat(ni, [prod(kern_sz) 1]);
      inds        = [coords ni_ind]; % 3D indices - kernel position & this pixel
      inds_vec    = subv2ind(size(visfield), inds );
      visfield( inds_vec ) = img_conv(:);
      
      if DEBUG
         indTranspose = subv2ind(ret_sz,fliplr(ind2subv(ret_sz,ni)));
         lowerx(ni) = coords(1); upperx(ni) = coords(2); 
         lowery(ni) = coords(3); uppery(ni) = coords(4); 
         figure(fh);
         subplot(ret_sz(1),ret_sz(2),indTranspose), 
            imagesc((img),[inputmin inputmax]); axis image; 
            title(sprintf('%d ',retRel),'fontsize',8);
            axis off;
            
         figure(f2);
         subplot(ret_sz(1),ret_sz(2),indTranspose), 
            plot(img(:,round(ret_sz(1)/2))); 
            title(sprintf('%d ',retRel),'fontsize',8);
%             axis off;
         figure(fr);
         subplot(ret_sz(1),ret_sz(2),ni), 
            imagesc((rfi)); axis image; % colorbar; 
            title(sprintf('%d ',indVec),'fontsize',8);
            axis off;
         figure(fc);
         subplot(ret_sz(1),ret_sz(2),ni), 
            imagesc((img_conv)); axis image; % colorbar; 
            title(sprintf('%d ',indVec),'fontsize',8);
            axis off;
      end
   end
   if DEBUG
      figure(fp); imagesc(output);
      figure(fh); set(fh,'name','Input');
      figure(f2); set(f2,'name','2D input (over y)');
      figure(fr); set(fr,'name','Kernel');
      figure(fc); set(fc,'name','Input convolved with kernel');
      figure(fp); set(fp,'name','Sum of convolution'); colorbar;
   end
end

% Get input image patch for given retinal cell position. Image patch
% assumed to be centred around the cell 
% Inputs:
%  in_sz    - size of visual input
%  kern_sz  - size of receptive field kernel
%  loc      - location of retinal cell given wrt input size
% Outputs:
%  img      - input image patch for retinal cell
%  ind      - visual input indices (for debugging)
function [img, ind] = getCellInput(input, kern_sz, loc)
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
%       ind = [x-kXl:x+kXu, y-kYl:y+kYu];
      ind = allCombos(kern_sz);
      ind = ind + [kXu kYu];
      return;
   end
   % if an edge goes over/under then need to wrap around
   ind = allCombos(kern_sz);
   ind = bsxfun(@minus, ind, [kXl kYl]); 
   ind = bsxfun(@plus,  ind, [  x   y]); 
   ind = bsxfun(@mod, ind-1, in_sz) + 1; 
   img = reshape(input(subv2ind(in_sz,ind)),kern_sz);
   
   % return indices for debugging
   
   
end






