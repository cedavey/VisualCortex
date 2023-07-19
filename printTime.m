% s = printTime(t,label)
% Print time to a string or, if no output variables provided, to standard
% out. 
function varargout = printTime(t,label)
   if ischar(t)
      tmp   = t;
      t     = label;
      label = tmp;
   end
   if nargin<2
      label = [];
   end
   if nargin<1 || isempty(t)
      try
         t = toc;
      catch
         fprintf('\nMust provide a valid time or preset tic\n');
         return;
      end
   end
   hours = fix(t/60^2);
   t     = t - hours*60^2;
   mins  = fix(t/60);
   t     = t - mins*60;
   secs  = round(t);
   t     = t - secs;
   msecs = t*1e3;
   sh    = ternaryOp(hours~=1,'hours',  'hour');
   sm    = ternaryOp(mins~=1, 'minutes','minute');
   ss    = ternaryOp(secs~=1, 'seconds','second');
   sms   = ternaryOp(msecs~=1, 'milliseconds','millisecond');
   
   if hours>0
      s  = sprintf('%s %d %s, %d %s, %d %s\n',label,hours,sh,mins,sm,secs,ss);
   elseif mins>0
      s  = sprintf('%s %d %s, %d %s\n',label,mins,sm,secs,ss);
   elseif secs>0
      s  = sprintf('%s %d %s\n',label,secs,ss);
   elseif msecs>0
      s  = sprintf('%s %2d %s\n',label,round(msecs),sms);
   else
      s  = sprintf('%s 0 time\n', label);
   end
   if nargout==0
      fprintf('%s',s);
   else
      varargout{1} = s;
   end
end
