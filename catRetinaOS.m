% [theta, pdf] = catRetinaOS(N)
% 
% Generate the distribution for orientation selectivity in the cat retina,
% based on the number of required angles, and return an angle set & its pdf
function [theta, pdf] = catRetinaOS(N, doplot)
   if nargin<2 || isempty(doplot), doplot = false; end
   dir      = fullfile(userpath,'visual_cortex','Errol');
   fname    = 'orientation_probability_curve.csv'; 
   fullname = fullfile(dir,fname);
   
   OSdata   = csvread(fullname,2); % ignore header (data must be numeric to use this fn)
   Ndata    = size(OSdata,1);
   
   if nargin<1 || isempty(N), N = Ndata; end
   
   
   if N ~= Ndata   
      % resample fn assumes points before & after sample vector are 0, which
      % can cause large inaccuracies at end points. To avoid this, make period
      % copies of the data, & then resample at required points
      theta    = OSdata(:,1);
      pdf      = OSdata(:,2);
      theta    = [theta - 2*pi;  theta;  theta + 2*pi];
      pdf      = [pdf; pdf; pdf];
      
      pdfNew   = resample(pdf, N*3, Ndata*3, 1);
      theta    = linspace(0, 180, N+1)';
      theta    = theta(1:end-1);
      pdf      = pdfNew(N+1:2*N);

   else % N == Ndata
      theta = OSdata(:,1)/pi*180;
      pdf   = OSdata(:,2);
   end
   
   pdf      = pdf / sum(pdf); % ensure if sums to 1
   
   if doplot
      figure;
      plot2y(OSdata(:,1)/pi*180, OSdata(:,2), theta, pdf);
%       plot(OSdata(:,1)/pi*180, OSdata(:,2)); 
%       hold on; 
%       plot(theta, pdf);
      legend('Original','Resampled');
   end
   
   % Original data from Schall
   x = (5:10:85)'; 
   y = [0.173952, 0.203471, 0.119205, 0.110095, 0.016680, 0.080624, 0.100368, 0.085423, 0.110181]';
   x = [x; x+90];
   y = [y; y(end:-1:1)];
end












