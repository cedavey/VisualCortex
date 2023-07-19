% Draw fig for paper with ellipse from each retina cell to each v1 cell to
% which it has a connection, with the size of the ellipse being
% proportional to the number of connections it has to the v1 cell

% You MUST run calculateAvgLGNInputOrientation first to gen data structs
% with genData set to TRUE
T = true; F = false;

nr = 7; nc = 7;       % number of subplots / V1 cells to show RGC inputs to V1 for

% Params for RGC to V1 connections figures
genRetCxs   = T;      % generate figure showing RGC conns to V1 cells & RGC orientation
minSynCxs   = 2;      % must have at least this many synapses to count as connected
encodeWgt   = 'radius'; % encode cx strength from RGC to V1 by 'radius' or line 'width'
incDeadSyn  = true;   % to include dead synapses as greyed out ellipses set to true

% Params for legends
genScaleLeg = T;      % generate ellipse legend showing ellipse size per num synapses
genOSLeg    = F;      % generate legend for colour coded orientations
ellLegVert  = F;      % plot ellipse sizes vertically (vertical list of ellipses)
nEllLegend  = 1;      % number of ellipse scales to include in legend

saveFigs    = F;      % save all figures generated
closeOnSave = F;      % close figures after saving them (not used if save=F)

if ~exist('startWgtdRetCx','var')
   calculateAvgLGNInputOrientation;
   startWgtdRetCx  = data.start.wgtdRetCx;
   perc30WgtdRetCx = data.perc30.wgtdRetCx;
   endWgtdRetCx    = data.end.wgtdRetCx;
end
% V cell ID & params
szV1Layer   = sz.(v1);

% Ellipse parameters
Nellipse    = 50; % number of samples to make the ellipse circumf. from
majorAxis   = 10; % length of principal axis of ellipse 
minorAxis   = 5;  % length of minor axis of ellipse 

% Retina cell grid parameters
szRetLayer  = sz.(retina);
sizeRetCell = majorAxis * 2 + 1;   % space needed for single retina cell
numRetCells = radius.([lgn v1])*3; % num retina cells in figure 
if ~isodd(numRetCells), numRetCells = numRetCells + 1; end % enable centring of V1 on grid
ret_startCx = startWgtdRetCx;             % starting # cxs from each rgc to each v1 
ret_P30Cx   = perc30WgtdRetCx;            % connections at 30%
ret_endCx   = endWgtdRetCx;               % ending # cxs from each rgc to each v1 
ret_theta   = outstruct.layerconfig.RF.theta; % bias of rgc 
ncxs        = sort(ret_startCx(:));
% get all possible angles, & for each make an ellipse
thetaSet    = cell2mat(outstruct.layerconfig.RF.angleSet);
% thetaSet(thetaSet==180) = 0; 
thetaSet = sort(thetaSet);
[ellX, ellY]= ellipse(majorAxis, minorAxis, thetaSet, [0 0], 50, 0);
cmap        = hsv(length(thetaSet)+1);    % circular colour maps
cmap        = cmap(1:length(thetaSet),:); % don't want 1st & last colour to be the same
% dead synapses depicted by a small, grey dot
[deadSynX, deadSynY]= ellipse(minorAxis, minorAxis, 0, [0 0], 50, 0);
nTheta      = length(thetaSet); 
thetaAxis   = (1:(sizeRetCell+1):((sizeRetCell+1) * (nTheta+1)));

% For the v1 cell in question, plot a square for each retina, & then in the
% square put a rotated ellipse, according to orientation of retina cell, &
% make the linewidth proportional to the number of connections to v1.

% make the grid - retinal cell positions depend on position number & grid size
% need an extra line to make numRetCells squares --> numRetCells+1
retinaAxis  = (1:sizeRetCell:(sizeRetCell * (numRetCells+1))) - floor(numRetCells/2)*sizeRetCell;
[gridX, gridY] = meshgrid(retinaAxis, retinaAxis);

% prep for 2 figs - 1 for ret -> v1 before plasticity, & 1 after
if genRetCxs
   fa = figure; fb = figure; f30 = figure; 
   figsV1 = [fa f30 fb]; % figures of V1 synaptic inputs from retina (via LGN)
end
% legend gets made in another figure since it requires drawing ellipses
if genScaleLeg
   fls     = figure;    % legend for number of synapses vs ellipse size before plasticity
   fle     = figure;    % legend for number of synapses vs ellipse size after plasticity
   fl30    = figure; 
   figsLeg = [fls fl30 fle]; % put legend figures in array so can loop through them
end
grey = 0.7; grey = [grey grey grey]; % colour of mesh grid lines, & synapses with 0 weight
LW   = 1; % default line width

for vID=1:(nr*nc) % prod(szV1Layer)
   % location of V1 cell relative to retina layer
   [~,~,loc_rel]= getGroupsLocationArray( szV1Layer, vID, szRetLayer );
   % gotta flip x axis to match imagesc, which is spatial layout
   loc_rel(1)    = sz.(retina)(1) - loc_rel(1);

   centreShift  = loc_rel; % - sizeRetCell/2;

   % for this v1 cell, find all retina inputs & draw ellipse in approp. spot
   retStart  = ret_startCx(:, vID);
   retPerc30 = ret_P30Cx(:, vID);
   retEnd    = ret_endCx(:, vID);
   retCxs    = [retStart retPerc30 retEnd];

   % get all connections at start of simulation - plot all ellipses
   % coloured according to angle for start of sim, & for end of sim
   % plot as grey if no longer connected, & colour coded if still connected
   retIDs   = find(retStart > minSynCxs);

   % for before (start), perc30, & after (end) plasticity, plot retinal connections
   for ff=1:3
      % get max & min number of connections before or after plasticity
      totalCxs   = sum(retCxs(retIDs,ff)); % sums up all synapses before/after
      maxPropCxs = max(retCxs(:,ff)) / totalCxs; % max percentage of connections
      
      % Generate plot of retinal connections to V1 if requested by the user
      if genRetCxs
         figure(figsV1(ff)); 
         subplot(nr, nc, vID), hold on;

         plot(gridX, gridY, 'color', grey, 'linewidth', 1);
         plot(gridY, gridX, 'color', grey, 'linewidth', 1);
         axis square;

         for ri=1:length(retIDs)

            % get retina ID & convert to position in grid
            retID        = retIDs(ri);    % get retina ID
            retPos       = ind2subv(sz.(retina), retID);
            % gotta flip x axis to match imagesc, which is spatial layout
            retPos(1)    = sz.(retina)(1) - retPos(1);
            % gotta account for periodic boundaries
            retPosCent   = getWrappedCentredCellPosition(sz.(retina), loc_rel, retPos);
            retGridPos   = (retPosCent-1) * sizeRetCell + 1 + (sizeRetCell/2);
            retAngle     = mod(ret_theta(retID), 180);
            angleID      = find(thetaSet == retAngle);
            
            % determine if synapse weights from this rgc have all gone to 0
            if retCxs(retID,ff)<= minSynCxs
               synIsDead = true;
            else
               synIsDead = false;
            end

            % Before plasticity, or after but still connected --> colour
            % according to angle, & use before/after number of conns
            if ~synIsDead
               ellcolour  = cmap(angleID, :); % colour of ellipse
               edgecolour = [0 0 0];          % colour of edge of ellipse
               propRetCxs = retCxs(retID, ff) / totalCxs;% num conns before/after plasticity
               ell        = [ellX(:, angleID) ellY(:, angleID)];
               ellScale   = propRetCxs / maxPropCxs;  % scale ellipse proportional to num cxs

            % After plasticity & no longer have any conns --> make grey,
            % with size from before plasticity if user wants to include in
            % after plasticity figure
            elseif incDeadSyn
               ellcolour  = grey;   % colour of ellipse
               edgecolour = grey;   % colour of edge of ellipse
               ell        = [deadSynX(:) deadSynY(:)];
               ellScale   = 0.3;    % set dead syn to constant size for small dot
               
               % for dead syn give before plasticity ellipse scale but make it grey
%                propRetCxs = retCxs(retID, 1) / outstruct.layerconfig.N.([lgn v1]); 

            % After plasticity with dead synapses, but user doesn't want to plot it
            else
               ellScale   = 0;
            end

            if strcmpi(encodeWgt, 'width')
               lw        = LW * ellScale / 3;
               plot(ell(:,1) + retGridPos(1), ell(:,2) + retGridPos(2), ...
                              'color', ellcolour, 'linewidth', lw);
            elseif strcmpi(encodeWgt, 'radius')
               lw        = LW;
               ell       = ell * ellScale;
               patch( ell(:,1) + retGridPos(1), ell(:,2) + retGridPos(2), ellcolour, ...
                           'EdgeColor', edgecolour); % 'none'
            end

         end % end for each rgc with synapse to the current v1 cell
         
         xlim([retinaAxis(1) retinaAxis(end)]);
         ylim([retinaAxis(1) retinaAxis(end)]);
         set(gca,'xcolor',[1 1 1]);
         set(gca,'ycolor',[1 1 1]);
      end

      % Generate plot of legend for ellipse scale, showing percentage 
      % synapses wrt ellipse size, if requested by the user
      if genScaleLeg
         figure(figsLeg(ff)); 
         subplot(nr, nc, vID), hold on;

         maxPercCxs = round(max(retCxs(:,ff)) / totalCxs * 100); % convert from num synapses to % of synapses
         minPercCxs = min([round(retCxs(retCxs(:,ff) > minSynCxs,ff) / totalCxs * 100); 1]);
         ellSz      = (1:nEllLegend)*(maxPercCxs - minPercCxs)/nEllLegend + minPercCxs;
         ellAxis    = (1:(sizeRetCell+1):((sizeRetCell+1) * (nEllLegend+1)));
         [gridThX, gridThY] = meshgrid([1 sizeRetCell], ellAxis);

         axis image;
         set(gca,'xcolor',[1 1 1]);
         set(gca,'ycolor',[1 1 1]);

         % gotta add 180 degree ellipse
         for ti=1:nEllLegend
            ellScale      = ellSz(ti) / maxPercCxs;
            ellGridPos   = (ti-1) * (sizeRetCell+1)*2 + 1 + (sizeRetCell/2);
            lw           = 3;

            % make 180 degree the same as 0 degrees but give it's own grid spot
            ellind       = 1;
            ell          = [ellX(:,ellind) ellY(:,ellind)];

            % plot diff ellipse sizes vertically 
            if ellLegVert
               x  = sizeRetCell/2; y = ellGridPos;
               xl = -3;           yl = ellGridPos;
            % plot diff ellipse sizes horizontally 
            else
               y = sizeRetCell/2; x = ellGridPos;
               yl = -3;          xl = ellGridPos + sizeRetCell/4;
            end

            if strcmpi(encodeWgt, 'width')
               lw        = LW * ellScale / 3;
               plot(ell(:, 1) + x, ell(:, 2) + y, 'color', grey, 'linewidth', lw);
            elseif strcmpi(encodeWgt, 'radius')
               lw        = LW;
               ell       = ell * ellScale;
               patch(ell(:, 1) + x, ell(:, 2) + y, grey, 'EdgeColor', grey); % 'none'
            end

            text(xl, yl, [num2str(ellSz(ti)) '%'],'fontsize',13,'fontname','Times','HorizontalAlignment','right');
         end

         % set y-axis to thetaSet values - scale theta's by max ellipse height
         set(gca,'ytick', thetaSet/thetaSet(end)*thetaAxis(end)-(sizeRetCell/2))
         set(gca,'yticklabel',num2str(thetaSet))
         ah = gca;
         
         % resize axis so scale of ellipse is same as for rgc cxs figure
%          set(gca,'position',[p(1) p(2) p(3)*nEllLegend/numRetCells p(4)*nEllLegend/numRetCells])

      end % end if genLegend

   end % end for each iteration - before/after plasticity

end % end for each V1 cell

% give figures a title & save if requested
for ff=1:3
   if ff == 1
      phase = 'start';
   elseif ff == 2
      phase = 'perc30'; 
   elseif ff == 3
      phase = 'end';
   end
      
   if genRetCxs
      set(figsV1(ff), 'name', sprintf('RGC connections to V_1 %s plasticity', phase ));
      fname = sprintf('retinalOrientsToV1_%s', phase );

      if saveFigs
         % for each phase of plasticity (before & after), save the figure
         fname = sprintf('retinalOrientsToV1_%s', phase );
         saveFigure(figure(figsV1(ff)),fname,1,{'fig','pdf'});
         if closeOnSave
            close(figure(figsV1(ff)));
         end
      end   
   end
   
   if genScaleLeg
      set(figsLeg(ff), 'name', sprintf('Ellipse scales %s plasticity', phase ));
      fname = sprintf('EllipseScaleLegend_%s', phase );

      if saveFigs
         % for each phase of plasticity (before & after), save the figure
         fname = sprintf('EllipseScaleLegend_%s', phase );
         saveFigure( figure(figsLeg(ff)), fname, 1, {'fig','pdf'});
         if closeOnSave
            close(figure(figsLeg(ff)));
         end
      end 
   end
end

% generate legend of ellipses
if genOSLeg
   % make it circular, so starts with 0, ends with 180
   if ~any(thetaSet==180), thetaSet = [thetaSet; 180]; end
%    nTheta      = length(thetaSet); 
%    thetaAxis   = (1:(sizeRetCell+1):((sizeRetCell+1) * (nTheta+1)));
   [gridThX, gridThY] = meshgrid([1 sizeRetCell], thetaAxis);

   fl = figure; hold on; % legend figure handle
   % Plot the
%    arrayfun(@(i) plot(gridThX(i,:), gridThY(i,:), 'color', grey, 'linewidth', 1), 1:(nTheta+1));
%    plot([gridThX(1) gridThX(1)], [gridThY(1) gridThY(end)], 'color', grey, 'linewidth', 1);
%    plot([gridThX(end) gridThX(end)], [gridThY(1) gridThY(end)], 'color', grey, 'linewidth', 1);
   axis image;
   set(gca,'xcolor',[1 1 1]);
   set(gca,'ycolor',[1 1 1]);

   % gotta add 180 degree ellipse
   for ti=1:nTheta
      ellGridPos   = (ti-1) * (sizeRetCell+1) + 1 + (sizeRetCell/2);
      theta        = thetaSet(ti);
      angleID      = find(thetaSet == retAngle);
      lw           = 3;

      % make 180 degree the same as 0 degrees but give it's own grid spot
      ellind       = ternaryOp(ti==nTheta, 1, ti);
      ell          = [ellX(:,ellind) ellY(:,ellind)];
      if strcmpi(encodeWgt, 'width')
         lw        = LW * propRetCxs / 3;
         plot(ell(:, 1) + (sizeRetCell/2), ell(:, 2) + ellGridPos, ...
                     'color', cmap(ti, :), 'linewidth', lw);
      elseif strcmpi(encodeWgt, 'radius')
         lw        = LW;
         patch(ell(:, 1) + (sizeRetCell/2), ell(:, 2) + ellGridPos, ...
                     cmap(ellind, :), 'EdgeColor', 'k'); % 'none'
      end

      text(-3, ellGridPos, num2str(theta),'fontsize',13,'fontname','Times','HorizontalAlignment','right');
   end
   
   % set y-axis to thetaSet values - scale theta's by max ellipse height
   set(gca,'ytick',thetaSet/thetaSet(end)*thetaAxis(end)-(sizeRetCell/2))
   set(gca,'yticklabel',num2str(thetaSet))
   ah = gca;
      
   if saveFigs
      fname = 'OrientationLegend';
      saveFigure(fl,fname,1,{'fig','pdf'});
      if closeOnSave
         close(fl);
      end
   end
end


















