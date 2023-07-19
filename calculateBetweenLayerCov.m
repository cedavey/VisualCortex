% Q = calculateLinskerCov(rates.(pre),dt)
% Calculate coveriance between input rate vector & an output neuron, 
% assuming stationarity. Fills in one column of the cov(A_i,B_j) matrix 
% - i.e. covariance between output j (in layer B) & all inputs i in layer A
% (at a delay specific to output j, which is why it's done for a single 
% postsynaptic neuron per function call.)
% Note: when you change f0 you need to change how you're scaling/normalising Q
function [Q, sig_pre, sig_post, mu_in, mu_out] = ...
               calculateBetweenLayerCov(Q_prev, sigin_prev, sigout_prev,...
                                        muin_prev, muout_prev,...
                                        inrates, outrate, tindex, dt, normQ)
   % TO DO: 
   % - figure out how to have Q & sig from multiple postsyn neurons, &
   %   a single presyn neuron
   % - do we need inrate mean to be calculated separately for each output
   %   neuron? Mathematically no, since signal's stationary.
   
   % Estimate covariance between cells in different layers --> need mean
   % rates etc from each layer
   % Cells have diff delays so need to use loops

   if any(isinf([inrates(:); outrate])) || any(isnan([inrates(:); outrate]))
      display('problem with input rates being Inf or NaN');
   end
   
   time = tindex * dt;
   
   % get mean of input & output rates: mean = old_mean*(nT-1)/nT + rates/nT
   inrates = toVec(inrates); 
   mu_in   = muin_prev  + (inrates - muin_prev )/time;
   mu_out  = muout_prev + (outrate - muout_prev)/time;
   
   % calculate total sum of distance from mean, squared - doesn't suffer
   % from catastrophic cancellation
   if tindex<=1
      mse_in  = (inrates - muin_prev).^2;
      mse_out = (outrate - muout_prev).^2;
   else
      mse_in  = sigin_prev.^2 *(time-dt) + (inrates - muin_prev ) .* (inrates - mu_in) * dt;
      mse_out = sigout_prev.^2*(time-dt) + (outrate - muout_prev) .* (outrate - mu_out) * dt;
   end
   % calculate real-time average of rate variance
   sig_pre    = sqrt(mse_in /time); % convert mse to var, and then take square root
   sig_post   = sqrt(mse_out/time);
   
   % For Q normalised, it should tend to 1 0 in most parts, and 1 where
   % there's a connection btwn pre * post layers, or the variance if not
   % normalise, which is the mean for a Poisson process, i.e.
   % whatever mu was initialised to
   
   % calculate total mse for each required value, then scale by num samples
   if tindex<=1
      mse_Q      = (inrates - muin_prev)*(outrate - muout_prev)'; 
      mse_norm   = sig_pre*sig_post';
      
   else
      mse_Q      = Q_prev*(time-dt) + (inrates - mu_in)*(outrate - mu_out)' * dt; 
      mse_norm   = sigin_prev*sigout_prev'*(time-dt) + sig_pre*sig_post' * dt;
   end

   if normQ
      norm    = mse_norm / time; 
      if any(isinf(norm) | isnan(norm) | norm==0)
         cprintf('*Errors','norm is inf or nan, will give bad Q\n');
      end
      Q       = mse_Q./norm / time; 
   else
      Q       = mse_Q / time; 
   end

   if any(isinf(Q) | isnan(Q) | Q==0)
      cprintf('*Errors','Q is inf or nan\n');
      pause
   end

%    if tindex<=1
%       dQ = (inrates - muin_prev)*(outrate - muout_prev)';
%    else
%       dQ = (inrates - mu_in)*(outrate - mu_out)';
%    end
%    dQ(isnan(dQ)) = 0;
%    
%    % If normalising, convert covariance estimate to correlation
%    % estimate, which is cov divided by the stnd dev of each neuron's
%    % rate, which is estimated by the square root of the diagonal of the 
%    % covariance structure
%    if normQ
%       if tindex<=1
%          Q = dQ./(sig_pre*sig_post');
%       else
%          Q = (Q_prev*(tindex-1) + dQ./(sig_pre*sig_post')) / tindex;
%       end
%    else
%       if tindex<1
%          Q = dQ;
%       else
%          Q = (Q_prev*(tindex-1) + dQ) / tindex;
%       end
%    end

   if any(isnan(Q(:)))
      display('problem with Q - 1+ NaNs in covariance estimate');
      pause
   end
end



