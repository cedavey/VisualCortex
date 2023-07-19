%% Calculate average orientation of retina (layer A) to output (layer C) cells
T = true;  F = false;
genData      = T;
plotOrient   = F; % plot retina orientations in v1
plotRetinaID = F; % plot retina cell id in v1
plotWgtAng   = F; % debugging info
plotWgts     = F; % plot retina weights or connections in v1
plotRetHist  = T; % plot histogram of number of unique retinal inputs to v1
plotRetProp  = T; % plot proportion/strength of each retinal cx to v1 
plotHist     = F; 
saveFigs     = F;

network     = 'L'; % 'FF'; % 
iterations  = {'start','perc30','end'};

% Extract orientations and weights etc coming from retina through lgn to v1
if genData || ~exist( 'data', 'var' )
   if strcmpi(network,'FF') && exist('FFstruct','var')
      outstruct = FFoutstruct;
   elseif strcmpi(network,'L') && exist('Lstruct','var')
      outstruct = Loutstruct;
   elseif ~exist('outstruct','var')
      cprintf('Error','\nNeed to get data from somewhere! Exiting...\n\n');
      return;
      % do nothing - keep outstruct;
   end
   niter   = length( iterations );

   % Assumes that data comes from outstruct
   retina  = 'A';
   lgn     = 'B';
   v1      = 'C';
   sz      = outstruct.layerconfig.sz;
   nV1     = prod(sz.(v1));
   plastic = getPlasticLayers( outstruct );

   synloc  = outstruct.ntwkconfig.synloc;
   conns   = outstruct.ntwkconfig.conns;
   sz      = outstruct.layerconfig.sz;
   theta   = outstruct.layerconfig.RF.theta;
   radius  = outstruct.layerconfig.r;
   nR      = prod(sz.(retina));                  % num retinal cells
   nV      = prod(sz.(v1));                % num V1 cells
   nC      = outstruct.layerconfig.N.([lgn v1]); % num cxs from lgn->v1

   % Calculate avg retina input orientation to v1 before & after plasticity
   allocRetV1    = @() zeros( nR, nV );
   allocV1       = @() zeros( nV, 1 );
   wgtdRetCx     = allocRetV1();
   LGNThetaCx    = cell( nV, 1 );
   
   % Calculate avg retina input orientation to v1 before & after plasticity
   numRetConns   = allocRetV1();   % number of cxs from retina -> v1
   avgRetAngle   = allocV1();
   varRetAngle   = allocV1();
   
   maxDisplay    = 5; sf = [];
   resz          = @(x) reshape(x,[nC nC]);

   % go through every layer C cell, & work backwards to extract the number of
   % connections from each layer A cell to it, & calculate the average
   % orientation - first allocate memory for each iteration struct
   data.maxT  = length( outstruct.outtime );
   for ii=1:niter
      it = iterations{ii};
      switch it
         case 'start'
            tpt  = 1;
         case 'end'
            tpt  = length( outstruct.outtime );
         otherwise
            if strfind( iterations{ii}, 'perc' )
               tmp  = regexp( iterations{ii}, '\d+', 'match' );
               perc = str2double( tmp{1} );
               if perc <= 1 || 100 <= perc
                  str = sprintf( 'Iteration percentage (%d) must be between 0 and 100 \n',...
                                   perc );
                  cprintf( 'Errors*', str );
                  return;
               end
            end
            % get index of the requested percentage of time
            tpt = round( length( outstruct.outtime ) *perc/100 ); 
      end
      it = iterations{ii};
      data.(it).timept      = tpt;
      data.(it).wgtdRetCx   = allocRetV1(); 
      data.(it).v1orient    = allocV1();
   end
   
   for vi=1:nV1
      cxsLGNv1   = conns.([lgn v1]){vi}; % conns from lgn to outlayer
      retV1Theta = theta( cell2mat( conns.([retina lgn])(cxsLGNv1) ) );
      retV1Theta( retV1Theta==180 ) = 0;
      
      avgRetAngle(vi) = ( numRetConns(:,vi)'  * theta ) / nC; 
      varRetAngle(vi) = std( numRetConns(:,vi) .* theta ); 
      
      for ii=1:niter
         it = iterations{ii}; tpt = data.(it).timept;

         wgtsLGNv1  = outstruct.outweights.([lgn v1]){vi}(tpt,:)';       
         data.(it).v1orient(vi) = ( wgtdRetCx(:,vi)' * theta ) / sum( wgtsLGNv1 ); 
         
         retV1Theta = theta( cell2mat( conns.([retina lgn])(cxsLGNv1) ) );
         cxsAllRetToV1 = cell2mat( conns.([retina lgn])(cxsLGNv1) );
         
         % Get weighted angle coming in from retina to v1, weighted by
         % plastic weights from ret->lgn & lgn->v1, given per lgn cell
         % Number of connections from retina -> LGN, and from LGN -> V1 may be
         % different for each cell
         for ri=1:prod(sz.(retina))
            if ~any( lgn==plastic )
               % get non-plastic weights of any lgn cell with a connection from this retinal cell
               wgtsRetLGN = arrayfun( @(id) outstruct.ntwkconfig.weights.([retina lgn]){id}( conns.([retina lgn]){ id } == ri ) , cxsLGNv1, 'uni', false);
            else
               % get weights of any lgn cell with a connection from this retinal cell
               wgtsRetLGN = arrayfun( @(id) outstruct.outweights.([retina lgn]){id}( tpt, conns.([retina lgn]){ id } == ri ) , cxsLGNv1, 'uni', false);
            end
            % sum number of cxs from each retina cell to each v1 cell
            numRetConns(ri,vi)  = sum( arrayfun(@(id) sum( conns.([retina lgn]){id} == ri ), cxsLGNv1) );
            % get weighted sum of cxs from each retina cell to each v1 cell,
            % weighted by plastic weights from retina->lgn, and lgn->v1, with 
            % this iteration's weights
            wgtdRetCx(ri,vi)   = sum( arrayfun( @(id) sum( wgtsRetLGN{id} * wgtsLGNv1(id) ), 1:length(cxsLGNv1) ) );
            % get weighted angle coming in from retina to v1, weighted by
            % post-plasticity weights from retina->lgn & lgn->v1
            LGNThetaCx{vi}     = arrayfun( @(id)  sum( theta( conns.([retina lgn]){cxsLGNv1(id)} ) ...
                                                            * wgtsLGNv1(id) ) ...
                                                            / length( conns.([retina lgn]){cxsLGNv1(1)} ), ...
                                                            1:length(cxsLGNv1) )';

            % endWgtdRetCx(mi,ni)  = sum(arrayfun(@(id) sum(conns.([retina lgn]){lgnCxs(id)}==mi)*endWgts(id)  , 1:length(lgnCxs)));
         end
         % endRetOrient(vi)   = ( endWgtdRetCx(:,vi)'  * theta ) / sum(endWgts); 
         data.(it).wgtdRetCx(:,vi)  = wgtdRetCx(:,vi);
         data.(it).v1orient(vi)     = ( wgtdRetCx(:,vi)' * theta ) / sum( wgtdRetCx(:,vi) ) ; % / sum( wgtsLGNv1 ) 
         data.(it).LGNThetaCx{vi}   = LGNThetaCx{vi};
         data.(it).avgRetOrient(vi) = ( wgtdRetCx(:,vi)'* theta ) / sum( wgtsLGNv1 ); 
      end

   end
end

if plotRetProp
   Nstart    =  outstruct.layerconfig.N.([lgn v1]);
   Nend      = -outstruct.layerconfig.k1.([lgn v1]) / outstruct.layerconfig.k2.([lgn v1]) * Nstart;
   numV1     =  prod( outstruct.layerconfig.sz.(v1) );
   maxT      = data.maxT;
   
   Nprop     = cell(numV1,1);
   
   % calculate total weight, which is sum( weight * synapse )
   % the total number of connections will be proportion of total weight left *
   % initial number of connections
   for ii=1:niter
      it = iterations{ii};
      data.(it).totalWeight = sum( data.(it).wgtdRetCx ); 
      data.(it).N = Nstart * data.(it).totalWeight / data.start.totalWeight;
   end
   
   for ii=1:niter
      it  = iterations{ii};
      tpt = data.(it).timept;
      wgtdRetCx = data.(it).wgtdRetCx;
      w_min     = mean( wgtdRetCx( wgtdRetCx(:)>0 ) ); 
      
      N = data.(it).N; 
      for vi=1:numV1
         Nprop{vi} = ( wgtdRetCx( wgtdRetCx(:,vi) > w_min, vi ) / data.(it).totalWeight(vi) );
%          Nprop{vi} = ( wgtdRetCx( wgtdRetCx(:,vi)>0, vi ) / sum( wgtdRetCx(:,vi) > 0 ) );
      end
      data.(it).Nprop = toVec( cell2mat( Nprop ) );
   end
   
   nBins     = 10;
   % number of unique RGC cells at the end will be fewer, & hence the
   % proportion of total input will be higher for these cells --> get bins
   % from end data, since they will be more dispersed
   [~, bins] = hist( data.end.Nprop, nBins );
   
   figure; nr = 1; nc = niter; 
   fopts = {'fontsize',16,'fontweight','bold'};
   for ii=1:niter
      it = iterations{ii};
      itHist = hist( data.(it).Nprop, bins );
      itHist = itHist / ( sum( itHist ) );
      subplot(nr,nc,ii),
         bar(bins, itHist);
         xlabel('Proportion of synapses from RGC cell');
         ylabel('Proportion of RGC cells');
         title( it );
   end
end

if plotRetHist

   figure; nr=niter; nc=1; fi=1; 
   fopts   = {'fontsize',16,'fontweight','bold'};
   xmin    = Inf; xmax = -Inf;
   minProp = 0.05; % not interested in super small values
   for ii=1:niter
      it = iterations{ii};
      % num retina inputs before learning
      N = sum( data.(it).wgtdRetCx > minProp );
N(N>100) = [];
      if length(unique(N))<10
         bins = unique( N );
      else
         bins = 10; % number of bins so they're found automatically
      end
      [y, x] = hist( N, bins );
      y = y/sum(y);
      y( y <= minProp ) = 0;
      if max(x)>xmax, xmax = max(x); end
      if min(x)<xmin, xmin = min(x); end
      subplot(nr,nc,fi), fi=fi+1;
         bar(x, y);
         xlabel('Number of unique retinal inputs',fopts{:});
         ylabel('Proportion of V_1 cells', fopts{:});
         % xlim(xl); ylim(yl);
         % set(gca,'xtick',xtick);
         title( it );
   end
   ah = findobj( gcf,'type','axes' );
   set(ah,fopts{:});
   set( ah ,'xlim', [0 xmax] );
end

numFigs = plotWgts*niter + plotOrient*niter + plotRetinaID*niter + plotWgtAng*1; ff = 1;
fh      = zeros( numFigs, 1 );
figname = cell( numFigs, 1 );


% Plot true/false connections from retina to v1
if plotWgts
   for ii=1:niter
      fh(ff) = figure; 
      WgtFigs(ii) = ff; ff=ff+1;
   end
   
   set( fh(WgtFigs), 'visible','off' )
   maxDisplay = 5;
   minCxs     = 2;
   % Plot the data
   N = min([maxDisplay sz.C(1)]);
   M = min([maxDisplay sz.C(2)]);
   if sz.C>maxDisplay
      dN     = sz.C(1) / maxDisplay;
      dM     = sz.C(2) / maxDisplay;
      dNM    = round(dN*dM);
      ninds  = (1:(N*M))*floor(dN*dM);
   else
      ninds  = (1:prod(sz.C))';
   end
   for ii=1:niter
      it   = iterations{ii};
      tpt  = data.(it).timept;
      
      for vi=1:length(ninds)
      cxsLGNv1 = conns.([lgn v1]){ninds(vi)}; % conns from lgn to outlayer
      
      % include 1 entry for each retinal cell
         wgts = outstruct.outweights.([lgn v1]){vi}(tpt,:)';
         x    = NaN(sz.(lgn));
         x(cxsLGNv1) = wgts;
         x(x >= mean(wgts)) = true;
         x(x < mean(wgts))  = NaN;
         x    = genCentredWeightImage( sz.(v1), sz.(lgn), ninds(vi), x );
         
         figure( fh( WgtFigs(ii) ) );
         subplot(N,M,vi)
            pcolor( x ); 
            colormap hsv; axis image; shading flat; % set(gca,'clim',[0 1]);
            xlim([sz.(lgn)(1)/2 - radius.([lgn v1])*2  sz.(lgn)(1)/2 + radius.([lgn v1])*2]); 
            ylim([sz.(lgn)(2)/2 - radius.([lgn v1])*2  sz.(lgn)(2)/2 + radius.([lgn v1])*2]); 
   %          xlim([sz.(retina)(1)/2 - radius.([lgn v1])*1.5  sz.(retina)(1)/2 + radius.([lgn v1])*1.5]); 
   %          ylim([sz.(retina)(2)/2 - radius.([lgn v1])*1.5  sz.(retina)(2)/2 + radius.([lgn v1])*1.5]); 
            set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1]);
            if vi==length(ninds), colorbar; end

      end
      
      figure( fh(WgtFigs(ii)) ); set(gcf,'name', sprintf( '%s plasticity', it ));
      figname{ WgtFigs(ii) } = sprintf( '%sWeights_BC', it );
   end
   if saveFigs
      saveFigure( fh(WgtFigs), figname(WgtFigs), 'fig','eps','pdf' )
   end
   
end

% Plot retinal cell number that ultimately connects to v1
if plotRetinaID
   retIDFigs = zeros(niter,1);
   for ii=1:niter
      fh(ff) = figure; 
      retIDFigs(ii) = ff; ff=ff+1;
   end   
   set( fh(retIDFigs), 'visible','off' )

   maxDisplay = 5;
   % Plot the data
   N = min([maxDisplay sz.C(1)]);
   M = min([maxDisplay sz.C(2)]);
   if sz.C>maxDisplay
      dN     = sz.C(1) / maxDisplay;
      dM     = sz.C(2) / maxDisplay;
      dNM    = round(dN*dM);
      ninds  = (1:(N*M))*floor(dN*dM);
   else
      ninds  = (1:prod(sz.C))';
   end
   for ii=1:niter
      figure( fh( retIDFigs(ii) ) );
      it        = iterations{ii};
      tpt       = data.(it).timept;
      wgtdRetCx = data.(it).wgtdRetCx;
      
      for vi=1:length(ninds)
      % include 1 entry for each retinal cell
         wgts      = outstruct.outweights.([lgn v1]){vi}(tpt,:)';
         x         = NaN(sz.(retina));
         x(wgtdRetCx(:,vi)>0) = find( wgtdRetCx(:,vi) > 0 );
         x         = genCentredWeightImage( sz.(v1), sz.(retina), ninds(ni), reshape(x, sz.(retina)) );

         figure(fh(retIDFigs(ii)));
         subplot(N,M,vi)
   %          imagesc(xstart,[0 180]); axis image; 
            pcolor( x ); 
            colormap hsv; axis image; 
            colorbar;
   %          set(gca,'clim',[1 900]); 
            shading flat;
            xlim([sz.(retina)(1)/2 - radius.([lgn v1])*2  sz.(retina)(1)/2 + radius.([lgn v1])*2]); 
            ylim([sz.(retina)(2)/2 - radius.([lgn v1])*2  sz.(retina)(2)/2 + radius.([lgn v1])*2]); 
            set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1]);
            if vi==length(ninds), colorbar; end
         
      end
         
      set( fh(retIDFigs(ii)), 'visible','on' )
      figure( fh(retIDFigs(ii)) ); set(gcf,'name', sprintf( 'Retina input ID at %s of plasticity', it) );
      figname{ retIDFigs(ii) } = sprintf( 'WgtdInputAngle_%s', it);

   end
   
   if saveFigs
      saveFigure( fh(orientFigs), figname(orientFigs), 'fig','eps','pdf' )
   end
end

% Plot orientations coming from retina through lgn and into v1
if plotOrient 
   orientFigs = zeros(niter,1);
   for ii=1:niter
      fh(ff) = figure; 
      orientFigs(ii) = ff; ff=ff+1;
   end   
   set( fh(orientFigs), 'visible','off' )

   maxDisplay = 5;
   % Plot the data
   N = min([maxDisplay sz.C(1)]);
   M = min([maxDisplay sz.C(2)]);
   if sz.C>maxDisplay
      dN     = sz.C(1) / maxDisplay;
      dM     = sz.C(2) / maxDisplay;
      dNM    = round(dN*dM);
      ninds  = (1:(N*M))*floor(dN*dM);
   else
      ninds  = (1:prod(sz.C))';
   end
   nTheta = length(outstruct.layerconfig.RF.angleSet);
   cmap   = hsv(nTheta+1);    % circulare colour map
   cmap   = cmap(1:nTheta,:); % don't want 1st & last to be the same
   
   for ii=1:niter
      figure( fh(orientFigs(ii)) );
      it         = iterations{ii};
      tpt        = data.(it).timept;
      wgtdRetCx  = data.(it).wgtdRetCx;
      LGNThetaCx = data.(it).LGNThetaCx;

      for vi=1:length(ninds)
         cxsLGNv1 = conns.([lgn v1]){ninds(vi)}; % conns from lgn to outlayer
         x        = NaN(sz.(lgn)); 
         x(cxsLGNv1) = round( LGNThetaCx{vi} );
         x        = genCentredWeightImage( sz.(v1), sz.(lgn), ninds(vi), x );

         figure( fh( orientFigs(ii) ) );
         subplot(N,M,vi)
   %          imagesc(xstart,[0 180]); axis image; 
            pcolor(x); set(gca, 'ydir', 'reverse');
            colormap(cmap); 
            axis image; 
            set(gca,'clim',[0 162]); 
            shading flat;
            xlim([sz.(lgn)(1)/2 - radius.([lgn v1])*2  sz.(lgn)(1)/2 + radius.([lgn v1])*2]); 
            ylim([sz.(lgn)(2)/2 - radius.([lgn v1])*2  sz.(lgn)(2)/2 + radius.([lgn v1])*2]); 
            set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1]);
            if vi==length(ninds), colorbar; end

      end
      set( fh( orientFigs(ii) ), 'visible', 'on' )
      figure( fh( orientFigs(ii)) ); set( gcf, 'name', sprintf( '%s plasticity', it ) );
      figname{ orientFigs(ii) } = sprintf( 'WgtdInputAngle_%s', it );
   end
   if saveFigs
      saveFigure( fh(orientFigs), figname(orientFigs), fig','eps','pdf' );
   end
end

if plotWgtAng
   fh(ff) = figure; WgtAngFig = ff; ff=ff+1; 
   nr=ceil(sqrt(niter+1)); nc=ceil((niter+1)/nr); fi=1;
%    subplot(nr,nc,fi), fi=fi+1;
%       plot( avgRetAngle, 'o', 'linewidth', 3);
%       title('Mean retina to V_1 input orientation');
%       xlabel('Neuron ID');
%       ylabel('Orientation');

   subplot(nr,nc,fi), fi=fi+1;
      plot( numRetConns, 'o' );
      title('Number of retina inputs');
      xlabel('Retina ID');

   for ii=1:niter
      it = iterations{ii};
      subplot(nr,nc,fi), fi=fi+1;
         plot( data.(it).avgRetOrient, 'o' );
         title( sprintf( '%s weighted retina to V_1 orientation', it ) );
         xlabel( 'V_1 neuron ID' );
         ylabel( 'Orientation' );
   end
   
end

if plotHist
   maxDisplay = 3;
   minCxs     = 2;
   % Plot the data
   N = min([maxDisplay sz.C(1)]);
   M = min([maxDisplay sz.C(2)]);
   if sz.C>maxDisplay
      dN     = sz.C(1) / maxDisplay;
      dM     = sz.C(2) / maxDisplay;
      dNM    = round(dN*dM);
      ninds  = (1:(N*M))*floor(dN*dM);
   else
      ninds  = (1:prod(sz.C))';
   end
   fhist     = zeros(length(iterations),1);
   for ii=1:length(iterations)
      it = iterations{ii};
      fhist(ii) = figure;
      for vi=1:length(ninds)
         y = data.(it).wgtdRetCx(:,vi);
         
         figure( fhist(ii) );
         subplot(N,M,vi),
            nx    = sum( y >= minCxs);
            x     = find(y >= minCxs);
            bounds= [mean(x)-3*std(x), mean(x)+3*std(x)];
            out   = x<bounds(1) | bounds(2)<x;
            x(out)= [];
            th    = theta(x);
            bar(1:length(x), y(x));
               % xl    = [min(x)*0.9 max(x)*1.1];
               % xlim(xl);
               % give labels to the top nxLab bars
               nxLab = min([length(x) N]);
               if strcmpi(it,'end')
                  [~,xind] = sort(y(x),'descend');
                  xtick = sort( xind(1:nxLab) );
                  xlab  = sort( round( th(xind(1:nxLab)) ) );
               else
                  dx    = length(x) / nxLab;
                  xtick = round(1:dx:length(x));
                  xlab  = x(xtick);
               end
               set(gca,'xtick', xtick, 'xticklabel', str2legend([],xlab,[],'cell'));
               loc = ind2subv( sz.C, ninds(vi) );  % get neuron subscript
               title( sprintf( '%s [%d,%d]', it, loc(1), loc(2)) );

      end
      figname = sprintf( 'InputAngleHistogram_%s', it );
      if saveFigs
         saveFigure( fhist(ii), figname, 'fig','eps','pdf' )
      end
   end
end

