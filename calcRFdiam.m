% Calculate RF diameter using pseudo receptive fields generated using
% generatePseudoRFs, which gives a struct of RFs for each layer & for each
% iteration

% TO DO: input from entire layer so each postsynaptic neuron has a much
% larger RF, else we've biased the solution 
function diam = calcRFdiam( rfs, doplot )
   if nargin<2 || isempty( doplot ), doplot = false; end
   iterations = fieldnames( rfs );
   ind = strcmpi( iterations, 'debug' );
   iterations( ind ) = [];

   % rotate RF & then find 0 crossings
   for iti=1:length( iterations )
      iter   = iterations{iti};
      layers = fieldnames( rfs.(iter) );
      
      for li=1:length(layers)
         layer = layers{li};
         [szRF, N] = size( rfs.(iter).(layer) );
         szRF = round( sqrt( szRF ) );
         diam.(iter).(layer) = zeros( N, 1 );
         for i=1:N
            rf    = rfs.(iter).(layer)(:,i);
            rf    = reshape( rf, [szRF szRF] );
            bw    = rf; 
            mu    = mean( rf(:) ); 
            bw( rf>mu ) = 1; bw( rf<=mu ) = 0;
            stats = regionprops( bw, 'Orientation' );
            bwrot = imrotate( bw, -stats.Orientation );
            bwdif = diff( bwrot ); % find 0-crossings for each column
            dcol  = arrayfun( @(i) find( bwdif(:,i) < 0, 1 ) - find( bwdif(:,i) > 0, 1 ), 1:szRF, 'uni', false );
            % set diameter to max width 
            diam.(iter).(layer)(i) = max( cellfun( @(d) ternaryOp( isempty(d), 0, d ), dcol ) );
         end
      end
   end
   
   if doplot
      nbins   = 10;
      nlayers = length( layers ); niter = length( iterations );
      figure; nr = niter; nc = nlayers;
      for iti=1:length( iterations )
         iter   = iterations{iti};
         layers = fieldnames( rfs.(iter) );

         for li=1:length(layers)
            layer = layers{li};   
            
            subplot( nr, nc, (iti-1)*nlayers + li ),
               histogram( diam.(iter).(layer), nbins, 'normalization', 'probability' );
               title( sprintf( 'RF diam: %s %s', iter, layer ) );
         end
      end
   end
end