% Q = calculateLinskerCov(rates.(pre),dt)
% Calculate covariance between input rate signals, assuming stationarity
% Calculates mean, var & Q in real-time, so needs the current time index
% Note: when you change f0 you need to change how you're scaling/normalising Q
function [Q, sig, mu] = calculateVisualCortexCov(Q_prev, sig_prev, mu_prev, rates, tindex, normQ)
   % Calculate real-time average of rate mean
   rates  = toVec(rates);
   
   if isempty(mu_prev)
      mu       = ones(size(rates)) * mean(rates);
      mu_prev  = ones(size(rates)) * mean(rates);
      sig_prev = ones(size(rates)) * mean(rates);
   else
      mu  = mu_prev + (rates - mu_prev)/tindex;
   end

   % Calculate total sum of distance from mean, squared, then convert to
   % variance. This is a more stable real-time calculation that doesn't
   % suffer from catastrophic cancellation (errors getting big wrt values)
   if tindex<=1
      mse = (rates - mu_prev).^2;
   else
      mse = sig_prev.^2*(tindex-1) + (rates - mu_prev) .* (rates - mu);
   end
   % calculate real-time average of rate variance
   sig    = sqrt(mse/tindex); % convert mse to var, and then take square root
   
   % For Q normalised, it should tend to 1 on the diagonals, else it should
   % tend to the variance, which is the mean for a Poisson process, so
   % whatever mu was initialised to
   if tindex<=1
      dQ = (rates - mu_prev) * (rates - mu_prev)';
   else
      dQ = ((rates - mu)*(rates - mu_prev)');
   end
   dQ(isnan(dQ)) = 0;

   % If normalising, convert covariance estimate to correlation
   % estimate, which is cov divided by the stnd dev of each neuron's
   % rate, which is estimated by the square root of the diagonal of the 
   % covariance structure
   if normQ
      if isempty(Q_prev) || tindex<=1
         Q = dQ./(sig*sig');
      else
         Q = (Q_prev*(tindex-1) + (dQ./(sig*sig'))) / tindex;
      end
   else
      if tindex<1
         Q = dQ;
      else
         Q = (Q_prev*(tindex-1) + dQ) / tindex;
      end
   end

   % We don't need to scale by f02 IF we're taking standardised cov above
   %    Q = Q/f02; % scale by learning rate
   
   if isnan(max(Q(:)))
      error('problem with max Q');
   end
   if any(isnan(Q(:)))
      disp('problem with Q - 1 or more NaNs in covariance estimate');
   end
end



