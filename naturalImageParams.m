% Convert user input for natural image input into a parameter structure
% that we can pass to updateNaturalImageParams.m and
% generateVisualCortexInput.m
%
% If inputting natural images then get names of images to use now
% use user input param of number of timesteps to display each img for,
% with total number of timesteps in simulation, to determine num imgs
function noiseparams = naturalImageParams( layerconfig, noiseparams )
   tperimg    = noiseparams;     % number of timesteps per img
   blurred    = false;           % user blurred or linear images
   
   % if noiseparams already cell array then we've got img names already
   if isnumeric(noiseparams) 
      if ~blurred
         % van hateren linear images
         img_dir    = fullfile(userpath, 'hpc', 'vanhateren_linear');
         img_names  = dir([img_dir '/*.iml']);     % get list of image names in dir
      else
         % van hateren deblurred images
         img_dir    = fullfile(userpath, 'hpc', 'vanhateren_deblurred');
         img_names  = dir([img_dir '/*.imc']);     % get list of image names in dir
      end
      
      totalimgs  = size(img_names,1);           % total number of images in dir
      nimg       = round(layerconfig.T/layerconfig.dt) / tperimg; 
      nimg       = min(totalimgs, nimg+1);      % max number of imgs is totalimgs (+buffer!)
      if nimg==0
         str = sprintf('\nNo images found in %s, can''t continue\n\n', img_dir);
         cprintf('Keywords*', str);
          noiseparams = [];
         return;
      end
      img_inds   = randperm(totalimgs, nimg);   % get random selection of imgs
      images     = cell(nimg, 1);
      [images{:}]= img_names(img_inds).name;
      images{end+1}= img_dir;                   % put directory at end of list
      first_img  = readNaturalImageFile( fullfile(img_dir, images{1}) ); 
      if isempty(first_img)
         ntwkconfig = [];
         return;
      end

      % calculate margin around image so that we don't need to wrap
      % input, which would create wacky input images
      if isfield( layerconfig, 'retina' )
         RF = layerconfig.RF;
         margin = RF.kernel_size(1);
      else
         margin = ones(1,2)*max(struct2array(r));
      end
      margin = 0; % causing problems with rate.I sizes, so try get away without it
      % noiseparams: [timesteps/img margin img_list first_img]
      noiseparams= {tperimg margin images first_img};
      
   end
   
end