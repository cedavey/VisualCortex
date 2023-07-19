%       outvariables = struct2v(s, fields)
% Extract field values from a structure, s, & return to output. 
%       outvariables = struct2v(s)
% If no output variables provided, save the specified field values to the
% calling workspace, in variable names that match the fieldnames
% Inputs:
%   s      - struct containing field names & field values
%   fields - name of fields in struct that you want returned
% Outputs:
%   outvariables - names of output variables to assign field values to
function varargout = struct2v(s,varargin)
   fnames    = fieldnames(s);
   nfields   = length(fnames);
   varargout = cell(nargout,1);
   
   % if we haven't requested any specific fields, just fill up output until
   % we run out of space
   if nargin==0
      for fi = 1:nargout
         varargout{fi} = struct.(fnames{fi});
      end
      return;
   end
   
   % if no outputs save requested fields to caller workspace
   if nargout == 0 
      % if no inputs or outputs specified, get all fieldnames in s
      if nargin==1
          varargin = fieldnames(s);
      end
      for fi = 1:nargin
         try
            assignin('caller',varargin{fi},s.(varargin{fi}));
         end
      end
      return;
   end

   % if we've requested certain fields then find them & output them until
   % we run out of space
   varargout = cell(nargin,1);
   for fi=1:length(varargin)
      if any(strcmp(varargin{fi},fnames))
         varargout{fi} = s.(fnames{strcmp(varargin{fi},fnames)});
      else
        cprintf('*comments',fprintf('\nstruct2v: can''t find field name %s, returning empty cell\n',...
                                     varargin{fi}));
        trace = dbstack;
        cprintf('*comments',fprintf('\tstruct2v called by %s line %d\n\n',...
                                     trace(2).file,trace(2).line));
      end
   end
   if nargout < nargin
      varargout = varargout(1:nargout);
   end
%    nout = 1; nin = 1;
%    while nout <= nargout && nin <= length(fnames)
%       if any(strcmpi(fnames{nin},varargin))
%          varargout{nout} = s.(fnames{nin});
%          nout = nout + 1;
%       end
%       nin = nin + 1;
%    end
end
