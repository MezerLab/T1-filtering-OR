function fgOut = fgRetainIndices(fg,indices)
% fgRetainIndices returns the fg with only the fibers indicated by indices.
% indices: a list of indices, or a binary vector the same length as
% fg.fibers.
fgOut = fg;
if ~(isnumeric(indices)) && ~(islogical(indices) && length(indices) == length(fg.fibers))
    error('The "indices" argument must be numeric');
end
if isempty(indices)
    fgOut.fibers = {};
    warning('Returning empty fg');
    return
end
fgOut.fibers = fgOut.fibers(indices);

for i = 1:length(fgOut.params)
    if ~isfield(fgOut.params{i},'stat')
        continue
    end
    lengthSmallDim = sort(size(fg.params{i}.stat));
    lengthSmallDim = lengthSmallDim(1);
    if lengthSmallDim == 1
        fgOut.params{i}.stat = fgOut.params{i}.stat(indices);
    elseif lengthSmallDim == 2
        fgOut.params{i}.stat = fgOut.params{i}.stat(indices,:);
    else
        error('No support for stat field with more than 2 dimensions')
    end
    
end

if isfield(fgOut,'pathwayInfo')
    fgOut.pathwayInfo = fgOut.pathwayInfo(indices);
end