% theta = initVisualCortexOrientation(retina_kernel,sz,angles,prob)
%
% Initialise orientation across 
%    RF.kernel_size  = 10;
%    RF.sd_surr_hor  = 3; 
%    RF.sd_surr_ver  = 5;
%    RF.sd_cent_hor  = 3;
%    RF.sd_cent_ver  = 5; 
%    RF.coeff_cent   = 1.00e-4;
%    RF.coeff_surr   = 0.30e-4;
%    RF.ratio_surr_cent = 2.5;
%    RF.kernel       = retina_kernel(RF,doplot); % receptive field kernel

function [kernels,theta] = initVisualCortexOrientation(kernel, sz, angles, prob)
   % For each neuron in the layer, randomly allocate it a receptive field
   % angle from the angle set, according to the given probabilities
   
   N     = prod(sz);      % num cells
   u     = rand(N,1);     % uniform random number for each cell
   P     = cumsum(prob);  % cumulative probability giving cut-off prob of each angle
   a     = arrayfun(@(r) find(r<=P,1,'first'), u); % compare random num to cum prob
   isang = find(cellfun(@isnumeric,angles));       % assign numeric angles immediately
   theta = arrayfun(@(i) ternaryOp(any(i==isang),angles{i},NaN), a);
   
   % for radial angles we need to calculate them according to cell position
   rind  = find(isnan(theta));          % radial cell indices
   loc   = ind2subv(sz,rind);           % cell location vectors
   loc   = bsxfun(@minus,loc,sz/2);     % centre cell location 
   % angle from cartesian coordinates - if x coord is 0 then assume
   % position is along the y-axis, & give an angle of 90 degrees
   ang   = ternaryOp(loc(:,1)~=0, atand(loc(:,2)./loc(:,1)), 90);
   theta(rind) = ang;                   % populate angle vector with radial angle
   
   [szx,szy] = size(kernel);
   x = allCombos(szx,szy); 
   y = x(:,2); x=x(:,1);% image indices
   % if 1 kernel for all retinal cells, rotate this kernel for each theta
   if isnumeric(kernel)
      kernels   = arrayfun(@(th) imrotate(kernel,th,'bicubic','crop'), theta, 'uniformoutput',false);
   % can provide a different kernel for each retinal cell
   elseif iscell(kernel)
      kernels   = arrayfun(@(i) imrotate(kernel{i}, theta(i),'bicubic','crop'), 1:N, 'uniformoutput',false)';
   end
end






