function params = fgGetParams(fg,val)
% fgGetParams returns the params.stat of fg.
%
% val: could be (1) string with name of parameter to return.
%            or (2) index of the parameter
%
% e.g.
% t1var = fgGetParams(fg,'T1_Var');
%
% Note: it might be faster to input fg.params and not the whole fg
% structure.

if isnumeric(val)
    idx = val;
else
    paramNames = cell(1,length(fg.params));
    for pI = 1:length(paramNames)
        if ~isfield(fg.params{pI}, 'name')
            paramNames{pI} = 'MissingNameField';
        else
            paramNames{pI} = fg.params{pI}.name;
        end
    end
    idx = find(strcmp(paramNames,val));
end

if length(idx)==0
     warning(['No field called "' val '"']);
     params = ones(size(fg.params{1}.stat));
     return
end
if length(idx)>1
    warning(['More than one params field called ' val ' exists. Use numeric value to get the params field you want (' num2str(idx) '). Taking last one by default.']);
    idx = idx(end);
end

params = fg.params{idx}.stat;