function [indicesLowz zDist]= fgGetLowZIndices(cfg, fg)
% fgGetLowZIndices gets an fg (usually a fg with low T1 STD values, like fgLowStdR_notMidsag),
% and returns the indices of fibers with extremely low z coordinates. In some
% subjects with a bad WM mask, fibers go down through the Pons, to the
% Cerebellum, and from there "jump" to V1.
%
%
% cfg is a configuration struct that can include the following:
%
%   cfg.method: 'minz' (minimal z coordinate along the fiber), or 'meanz' (the mean z coordinate along the fiber)
%   cfg.minDist: the minimal distance to be considered a low-z fiber (between the LGN and the lowest z voordinate) (default 3 mm)
%   cfg.percentLength: percent of fiber length to use, between its most anterior and most posterior coordinates (default 0.5)

% Get configuration parameters (or set default ones)

if isfield(cfg,'method')
    method = cfg.method;
else
    method = 'minz';
end
if isfield(cfg,'minDist')
    minDist = cfg.minDist;
else
    minDist = 20;
end
if isfield(cfg,'percentLength')
    percentLength = cfg.percentLength;
else
    percentLength = 0.5;
end

fgOriginal = fg;
fg = getFgBetweenCoords(fg, percentLength);

% Get maximal number of nodes to save
nodesNum = 1;
for fI = 1:length(fg.fibers)
    nodesNum = max(nodesNum, length(fg.fibers{fI}));
end

x = nan(length(fg.fibers),nodesNum);
y = nan(length(fg.fibers),nodesNum);
z = nan(length(fg.fibers),nodesNum);

for fI = 1:length(fg.fibers)
    if abs(fg.fibers{fI}(2,1)) < abs(fg.fibers{fI}(2,end)) % If fibers coordinates start at the LGN and end at V1
        fInodesNum = size(fg.fibers{fI},2);
        
        x(fI,1:fInodesNum) = fg.fibers{fI}(1,:);
        y(fI,1:fInodesNum) = fg.fibers{fI}(2,:);
        z(fI,1:fInodesNum) = fg.fibers{fI}(3,:);
    else
        fInodesNum = size(fg.fibers{fI},2);
        
        x(fI,1:fInodesNum) = fg.fibers{fI}(1,end:-1:1);
        y(fI,1:fInodesNum) = fg.fibers{fI}(2,end:-1:1);
        z(fI,1:fInodesNum) = fg.fibers{fI}(3,end:-1:1);
    end
end
z1Median = nanmedian(z(:,1)); % the median z coordinate of all fibers (an approximation to the LGN)


% Create a mask on the coordinates, to only look at coordinates
% posterior to the LGN (otherwise the Meyer's loop segment might mess
% things up in the clustering)

mask = ones(size(x));
mask = y>=(repmat(y(:,1),1,nodesNum)-20);
zMasked = z;
zMasked(mask) = nan;

if strcmp(method,'minz')
    Z = nanmin(zMasked'); 
elseif strcmp(method,'meanz')
    Z = nanmean(zMasked');
end


zDist = z1Median-Z; % distance between the LGN (approximately) and the computed Z of each fiber

indicesLowz = find(zDist>minDist);
