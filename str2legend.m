% c = str2legend(str,vec,suffix,type)
% Returns a cell array for a legend entry, such that each legend entry in
% vec has 'str' appended to it. You can request a cell array in return, or 
% a character array. The cell array is to 3d.p.
% Inputs:
%    str    - string to prepend to numbers in vec
%    vec    - numbers to append to string
%    suffix - add suffix to each legend entry (optional)
%    type   - output 'char', or 'cell' (dep. on if strings are same length)
function c = str2legend(str,vec,suffix,type)
    if nargin<4
        type = 'char';
    end
    if nargin<3
       suffix = [];
    end
    if nargin==1, vec=str; str=[]; end
    if isempty(str)
       if isint(vec)
        c = sprintf('%d\t',vec);
       else
        c = sprintf('%0.2f\n',vec);
       end
    end
    
    switch type
        case 'cell'
           if isint(vec)
              c = arrayfun(@(d) sprintf('%s%d%s',str,d,suffix), vec,'uniformoutput',0);
           else
              c = arrayfun(@(d) sprintf('%s%0.2f%s',str,d,suffix), vec,'uniformoutput',0);
           end
        case 'char'
            vec = vec(:);
            c = [repmat(str,length(vec),1) num2str(vec)];
            if ~isempty(suffix)
               c = strcat(c,suffix);
            end
    end
end
