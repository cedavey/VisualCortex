% fit a Gaussian RF to response from white noise
function gaussRF = fitGaussToSpatialRF( rfs, outstruct, iterations, layers )
      inlayer    = outstruct.layerconfig.inlayer;
      layernames = outstruct.layerconfig.layernames;
      sz         = outstruct.layerconfig.sz; 

      possible_it= fieldnames(rfs);
      % ditch unrecognised iterations
      if isempty( iterations ), iterations = possible_it; end
      if isempty( layers ),     layers = layernames(end); end
      if ischar(iterations), iterations = {iterations}; end
      for iti=1:length(iterations)
         it = iterations{iti};
         if ~any( strcmpi( possible_it, it ) )
            iterations(iti) = [];
         end
      end
      if isempty( iterations)
         str = sprintf( 'fitGaussToSpatialRF - iterations not found in rfs, exiting...\n' );
         cprintf( 'error*', str );
         return;
      end
      % ditch unrecognised layers
      if ischar(layers), layers = {layers}; end
      for li=1:length(layers)
         lname = layers{li};
         if ~any( strcmpi( layernames, lname ) )
            layers(li) = [];
         end
      end
      if isempty( layers)
         str = sprintf( 'fitGaussToSpatialRF - layers not found in outstruct, exiting...\n' );
         cprintf( 'error*', str );
         return;
      end
      num_angles = 10; 

      d_angle    = 180 / num_angles;
      angles     = 0 : d_angle: (180-d_angle); 
      
   % take the position of each C neuron in layer B & calculate a line
   % emanating at a certain angle from this position. Extract grid
   % coordinates of this line, and get data from those coordinates. Fit a
   % Gaussian to the data. 
   insize = sz.(inlayer);
   line   = toVec( -max(insize):max(insize) );
   for iti=1:length(iterations)
      it = iterations{iti};

      for li=1:length(layers)
         post    = layers{li};

         for ai=1:num_angles
            theta = angles( ai ); % angle of line to calc Gauss RF across

            for i = 1:prod(sz.(post))
               pos = getRelativePosition( sz.(post), i, sz.(inlayer) );
               [~,mind] = max( toVec(rfs.(it).(post)(:,:,i)) ); 
               pos = ind2subv( sz.(inlayer), mind );
               x   = round( pos(1) + line * cosd(theta) ); % x-coord of line
               y   = round( pos(2) + line * sind(theta) ); % y-coord of line
               val = 0<x & x<=sz.(inlayer)(1) & 0<y & y<=sz.(inlayer)(2);
               x   = x(val);
               y   = y(val);
               l   = line(val);
               ind = sub2ind( size( rfs.(it).(post) ), x, y, ones(size(x))*i );
               r   = rfs.(it).(post)(ind);            % response at line coords
               try
                  f   = fit(l, r, 'gauss1'); 
                  a1  = f.a1; b1 = f.b1; c1 = f.c1;
                  sigma = c1/sqrt(2);
                  % FWHM = 2*sqrt(2*log(2)) * sigma;
                  FWHM = 2*sqrt(log(2)) * c1;
                  gaussRF.(it).(post).a1(ai,i)  = a1;
                  gaussRF.(it).(post).b1(ai,i)  = b1;
                  gaussRF.(it).(post).c1(ai,i)  = c1;
                  gaussRF.(it).(post).fwhm(ai,i)= FWHM;
               end
               
%                gaussRF.(it).(post).mu(i)  = mean(r);
%                gaussRF.(it).(post).std(i) = std(r);
%                gaussRF.(it).(post).pos(i,:) = pos;
            end
         end


      end
   end

end

