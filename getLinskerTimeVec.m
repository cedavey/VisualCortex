% [time,nT,nQ] = getLinskerTimeVec(T,dt,interval,endpoints,Qinterval)
function [time,nT,nQ] = getLinskerTimeVec(T,dt,interval,endpoints,Qinterval)
   if nargin<5 || isempty(Qinterval)
      Qinterval = interval;
   end
   if nargin<4 || isempty(endpoints)
      endpoints = false;
   end

   nT   = fix((T/dt)/interval);
   nQ   = fix((T/dt)/Qinterval);
   nQ   = ternaryOp(isnan(nQ) || isinf(nQ), 0, nQ);
   nT   = ternaryOp(isnan(nT) || isinf(nT), 0, nT);
   
   % Change interval to mean number of samples requested
   time = toVec((interval:interval:(nT*interval))*dt);
%    time = toVec((nT:nT:(nT*interval))*dt);

   % include endpoints
   if endpoints
      nT   = nT + 1; % don't need cov timept here cuz starts as identity matrix
      time = [0; time];
      % if endpoint is not already included
      if time(end)~=T
         nT   = nT + 1;
         nQ   = nQ + 1;
         time = [time; T];
      end
   end
