function ori_tun = ori_tuning_CV(tn, crv, modes, fit_curve)

% Returns preferred orientations and circular variance (CV) values for
% responses to various spatial frequencies and orientations at start and 
% end of plasticity simulation.
% 
% Takes output of [tun, curv] = processVisualCortexInput(outstruct)
% 
% For each cell, calculates preferred ori and CV for cell's response to
% each spatial frequency.
% 
% As for which layers, at which stage (start/end) and which response(s)
% (f1/f0) are used, is up to the user to select - see below 'filter data to
% analyse'.
% 
% Example
% -------------------------------------------------------------------------
% [tun, curv] = processVisualCortexInput(outstruct);
% ori_tun = ori_tuning_CV(tun, curv, outstruct, 2, 0);
% 
% ARGUMENTS
% -------------------------------------------------------------------------
%   tn:     'tuning' output from processVisualCortexInput(). Details 
%           parameters of the input stimuli used to generate cellular 
%           responses
% 
%   crv:    'curves' output from processVisualCortexInput(). Contains
%           output rates for each cell in each layer etc.
% 
%   modes:  Defaults to 2.  Passed into the VC_ori_tuning() function, which
%           actually does the calculating.  See docs for that function.
%           Defines how many preferred orientations are calculated.
%           Defaults to 2 as we're generally interested in a preferred
%           axis, with an angle between 0 and 180 degs and a corresponding
%           angle that is + 180degs around, forming an axis.  Ie, we
%           presume the response(orientation) function is bi-modal.
%           Theoretically, you could calculate the preferred orientations
%           for any amount of modalities, provided the responses are
%           sampled at the appropriate orientations.
% 
%   fit_curve: Defaults to false.  Generally unnecessary.  Fits von_mises
%               curves (standard orientation tuning curve function) to the
%               orientation responses.  May be useful when writing a paper
%               and want to cover all bases, but generally not substantially
%               different information.
%               Also rather slow - lots of curves to fit.
% 
%   other:  Under comment 'filter data to analyse' below, the particular
%           stages, layers and responses to be analysed are selected in the
%           dumb way of commenting or not commenting the lines the assign
%           the desired stages, layers etc.
% 
%           Currently defaults to f1 only, layers B and C, at start & end.
% 
% RETURNS
% -------------------------------------------------------------------------
% ori_tun:      Struct.
%               Fields:
%                   orientations: orientations of stimuli (rads)
%                   spatial_freqs: spatial frequencies of stimuli (cyc
%                   per ?)
%                   
%                   Then start>layer>f1>vect_ori_data for each of the
%                   stages, layers and responses selected to be analysed.
%               
% vect_ori_data:    leaf of the struct tree.  see docs for VC_ori_tuning
% 
%                   Matrix - 3D - (n cells x n spat_freqs x modes+CV)
%                   In third dimension, the last value is the Circular
%                   Variance (CV); preceding are preferred orientations,
%                   one for each mode.  Generally, only need the first, as
%                   modes defaults to 2, so second will always be
%                   first+180degs.
% 
% 
% DEPENDANCIES
% -------------------------------------------------------------------------
% VC_ori_tuning.m
% curve fit toolbox (if fit_curve)??


if nargin < 4
    fit_curve = 0;
end

% filter data to analyse
% create variable curv to store all curve output from 
% processVisualCortexInput that will be analysed
% Could be done better than just commenting out lines!!!

iterations = fieldnames(crv);
layers     = fieldnames(crv.(iterations{1}));
responsefn = fieldnames(crv.(iterations{1}).(layers{1}));

plasticLayers = getPlasticLayers(tn);

for it=1:length(iterations)
   for li=1:length(layers)
      layer = layers{li};
      % if 1st iteration, or if a postsynaptic layer of plastic cx
      if it==1 || any(layer==plasticLayers)
         for fi=1:length(responsefn)
            curv.(iterations{it}).(layer).(responsefn{fi}) = crv.(iterations{it}).(layer).(responsefn{fi});
         end
      end
   end
end

% get states of the curv data
iterations = fieldnames(curv);

% get layer names
layers = fieldnames(curv.(iterations{1}));

% get responses (f1 and/or f0)
responsefn = fieldnames(curv.(iterations{1}).(layers{1}));

orientations = tn.grating_angles./180 * pi;
spatial_freq = tn.spatial_freq;

if isrow(orientations)
    orientations = orientations';
end

ori_tun.orientations = orientations;
ori_tun.spatial_freqs = spatial_freq;


% define von mises function 
% a-max amplitude; k-width factor; phi-preferred orientation
von_mises = @(a,k,phi,x) a*exp( k* ( (cos(x-phi).^2) -1));

% generate fitting function
fit_vm = fittype(von_mises);


for st=1:length(iterations)
    
    for lay=1:length(layers)
       layer = layers{lay};
        
        if st==1 || any(layer==plasticLayers)

           for resp=1:length(responsefn)

   %             l_sz = prod(outstruct.layerconfig.sz.(layers{lay}));
               l_sz = size(tn.iterations.(iterations{1}).output{1}.(layer),2);

               vect_ori_dat = zeros(l_sz, length(spatial_freq), modes+1);

               curv_ori_dat = zeros(l_sz, length(spatial_freq), 6); % three parameters of vm curve, and two half-widths and rsquare

               for cell=1:l_sz

                   for sf=1:length(spatial_freq)

                       r = curv.(iterations{st}).(layers{lay}).(responsefn{resp});

                       [pref_ori, CV] = VC_ori_tuning(orientations, r(:,sf,cell), modes);

                       % There are as many preferred orientations as there are
                       % modes
                       vect_ori_dat(cell, sf,1:modes) = pref_ori;
                       vect_ori_dat(cell, sf,modes+1) = CV;

                       if fit_curve
                           [val, idx] = max(r(:, sf, cell));
                           [fo, gof] = fit(orientations, r(:,sf,cell), fit_vm, 'StartPoint', [val, 1, orientations(idx)], 'Lower', [0, 0.35, 0]);


                           curv_ori_dat(cell, sf, 1:3) = coeffvalues(fo); % order in definition of vm anon function (a, k, phi)
                           % half widths, to be 4:5, are done on the whole matrix
                           curv_ori_dat(cell, sf, 6) = gof.rsquare;

                       end


                   end
               end

               ori_tun.(iterations{st}).(layers{lay}).(responsefn{resp}).vect_ori_data = vect_ori_dat;

               if fit_curve
                   % calculate half widths
                   ks = curv_ori_dat(:,:,2);
                   curv_ori_dat(:,:,4) = 0.5 * acos( (log(0.5)+ks)./ks); % std
                   eta = 1 - ((1-exp(-ks))./2);
                   curv_ori_dat(:, :, 5) = acos( sqrt( (log(eta)./ks) + 1)); % real

                   ori_tun.(iterations{st}).(layers{lay}).(responsefn{resp}).curv_ori_data = curv_ori_dat;
               end

           end
        end % end 1st iteration or plastic layer
    end % end for each layer
end % end for each iteration



end




                
                
            
            

