%% TO DO:
% Cull edges of final weight figures for large layers, so cx's are not a
% tiny spec in the middle of the layer

%% Plot figures
% This code is designed to plot the weight timeseries as the simulation is
% evolving, and the final weights etc at the very end
warning off;
maxFigs   = 36;
fopts     = {'fontsize',6};
isplastic = toVec(struct2array(plastic)); 

for ci=1:size(layerconns,1)
   pre     = layerconns{ci,1};
   post    = layerconns{ci,2};
   cxlabel = [pre post];
   if plastic.([pre post])
      % Prepare figures - number of figs, titles, subplots etc (set invisible)
      np   = prod(sz.(post)); % num postsynaptic subplots to include
      np   = min([np maxFigs]);    % upper bound on number of subplots 
      delta= floor(prod(sz.(post))/np);
      nIDs = 1:delta:round(prod(sz.(post)));
      nr   = ceil(sqrt(np));  % num subplot rows in figure
      nc   = ceil(np/nr);     % num subplot cols in figure
      % if 1st output iteration, or we didn't run the simulation & we're 
      % just plotting, ft won't exist yet
      if ~exist('ft', 'var'), ft = NaN( size(isplastic) ); end      
      if outr == 1 || ~ishandle(ft(ci))
         ft(ci) = figure('visible','on');   % fig for t.series of weights for all neurons
         dstr   = ['[' pre ' ' post '] delay = [ ' sprintf('%g ',delparams.([pre post])) ']'];
         if any( delparams.([pre post]) > 0 )
            dstr = sprintf( '[%s %s] delay = [%g]', pre, post, delparams.([pre post]) );
         else
            dstr = sprintf( '[%s %s]', pre, post );
         end
         set(ft(ci),'name',['Weight over time: ',dstr]);
      end
      
      % Update weight figure as the simulation evolves
      for ni=1:np
         nID  = nIDs(ni); 
         wN   =  toVec(outweights.(cxlabel){nID}(end,:));
         cN   = conns.(cxlabel){nID};
         [~,col] = genCentredWeightImage(sz.(post), sz.(pre), nID, wN, cN);
         % summarise each neuron on its own figure - wgts/corr/lags etc
         set(0, 'currentfigure', ft(ci));
         subplot(nr,nc,ni),
            hold off; % draw line each time, else matlab won't connect the samples
            plot( outtime(1:outr), outweights.(cxlabel){nID}(1:outr,:) );
            lh = findobj(gca,'type','line');
            set(lh,{'color'},col);
            hold on; 
            plot( outtime(1:outr), mean( outweights.(cxlabel){nID}(1:outr,:), 2 ),...
                 'r', 'linewidth', 2.5 );
            xlim( [0 outtime(end)] );
            if thresholdWgts
               ylim( [n.(cxlabel)-1 n.(cxlabel)] );
            end
         % mult (post) layer neurons --> have separate fig. for final img & final wgts
      end
      refresh(ft(ci));
      set(ft(ci), 'units', 'normalized', 'outerposition', [0 0 1 1]);
      figure(ft(ci));
      
      % if finished simulation, plot final weights & cov matrices etc
      if ~outr == length(outtime), return; end
      fi      = NaN( size(isplastic) );
      fi(ci)  = figure('visible','off');   % fig for final weight of all neurons
      set( fi(ci), 'name', sprintf('Final weights: %s', dstr) );
      figs    = [ft; fi];
      layfigs = [ft(ci); fi(ci)];
      if calcCorr
         fc     = NaN( size(isplastic) );
         fc(ci) = figure( 'visible', 'off' ); % figure for correlation info
         layfigs= [ft(ci); fi(ci); fc(ci)];
         set( fc(ci), 'name', sprintf('Covariance: %s', dstr) );
         figs   = [figs; fc];
         numCov=4; % number of covariance related figures
         if np==1
            nrr = ceil(sqrt(np*numCov));
            ncc = ceil(np*numCov/nrr);
         else
            nrr = numCov; ncc = numCov; % only 4 neurons to display cov 
         end
      end

      noCxVal   = NaN; % value of neurons with no connections
      for ni=1:np
         nID    = nIDs(ni);
         wN     =  toVec(outweights.(cxlabel){nID}(end,:));
         cN     = conns.(cxlabel){nID};
         img    = zeros(sz.(pre));
         img(conns.([pre post]){nID}) = wN;

         [img_c,col] = genCentredWeightImage(sz.(post), sz.(pre), nID, ...
                        wN, cN, noCxVal);
         nanCol = n.(cxlabel)-2;
         img_c(isnan(img_c)) = nanCol;
         set( 0, 'currentfigure', fi(ci) );
         subplot(nr,nc,ni),
         if thresholdWgts
            imagesc( flipud(img_c), [nanCol n.(cxlabel)] ); colormap('gray'); axis image; 
         else
            imagesc( flipud(img_c) ); axis image; 
         end
         set(gca,'xcolor',[1 1 1],'ycolor',[1 1 1])
         % sc(img,[n.C-1 n.C],'gray',[0.5 0.5 0.5]);

         % covariance matrix only has entries where there's a connection
         if calcCorr
            if ni<=nrr
               nC      = N.(cxlabel);
               if ni==1, end
               xrho = zeros(nC);
               rho  = zeros(nC);
               mlags = zeros(nC);
   %             shared  = zeros(nC);
               numlags = min([T/dt 10]);
               for n1=1:nC
                  for n2=1:n1
                     s1 = outrates.(pre)(:,n1); 
                     s2 = outrates.(pre)(:,n2);
                     [cc,lags]   = crosscorr(s1,s2,numlags);
                     [~,ind]     = max(abs(cc));
                     mlags(n1,n2)= lags(ind); % cc has [-numlags numlags]
                     xrho(n1,n2) = cc(ind);
                     rho(n1,n2)  = (corr(s1,s2));
                  end
               end

%                clim = max(abs(rho(rho(:)<(1-eps*100))));
               clim = max(abs(rho(:)));

               set( 0, 'currentfigure', fc(ci) );
               subplot(nrr,ncc, (ni-1)*numCov+1),
               set(gca,fopts{:})
               q = squeeze(outQ.(cxlabel){nID}(end,1:nC,1:nC));
               imagesc(q); axis image; 
               title('Cov');

               subplot(nrr,ncc,(ni-1)*numCov+2),
               set(gca,fopts{:})
               imagesc(xrho,[-clim clim]); axis image; % colorbar; 
               title('Max \rho');

               subplot(nrr,ncc,(ni-1)*numCov+3),
               set(gca,fopts{:})
               imagesc(rho,[-clim clim]); axis image; % colorbar;
               title('\rho');

               subplot(nrr,ncc,(ni-1)*numCov+4),
               set(gca,fopts{:})
               imagesc(mlags); axis image; % colorbar;
               title('Lags');
            end
         end
      end
      set( layfigs, 'visible', 'on');
      figname = sprintf('finalWeights_%s', cxlabel);
      set( fi(ci), 'units', 'normalized', 'outerposition', [0 0 1 1] );
      saveFigure( figure(fi(ci)), figname, 1, {'fig','pdf'} );

      figname = sprintf('weightsOverTime_%s', cxlabel);
      set( ft(ci), 'units', 'normalized', 'outerposition', [0 0 1 1] );
      saveFigure( figure(ft(ci)), figname, 1, {'fig','pdf'} );
   end
end
figs( isnan(figs) ) = []; % get rid of invalid figs from non-plastic layers
set( figs, 'visible', 'on' );

