function lateralIndices = fgGetLateralIndices(cfg, fg)
% fgGetLateralIndices gets an fg (usually a fg with low T1 STD values, like fgLowStdR_notMidsag),
% separates it into a lateral and a medial bundle, and returns the indices of the
% lateral bundle.
%
%
% cfg is a configuration struct that can include the following:
%
%   cfg.method: 'xmin' (minimal x of each fiber in the segment),
%               'xmin_highstd' (same, but only for a small segment in y coordinates where the x coordinates are wide spread)
%   cfg.clusterMethod = 'kmenas', 'kmedoids' or 'gmm' (EM algorithm for Gaussian Mixture Model)
%   cfg.iterNum: number of iterations (default 1)
%   cfg.minDist: the minimal distance to separate medial and lateral bundles (default 3 mm)
%   cfg.percentLength: percent of fiber length to use, between its most anterior and most posterior coordinates (default 0.5)
%   cfg.plotFlag: 1 ('yes') or 0 ('no') (default 0)
% TODO: add about medialToLgnFlag

% Get configuration parameters (or set default ones)
if isfield(cfg,'method')
    method = cfg.method;
else
    method = 'xmin';
end
if isfield(cfg,'clusterMethod')
    clusterMethod = cfg.clusterMethod;
else
    clusterMethod = 'kmeans';
end
if isfield(cfg,'iterNum')
    iterNum = cfg.iterNum;
else
    iterNum = 1;
end
if isfield(cfg,'minDist')
    minDist = cfg.minDist;
else
    minDist = 3;
end
if isfield(cfg,'percentLength')
    percentLength = cfg.percentLength;
else
    percentLength = 0.5;
end
if isfield(cfg,'medialToLgnFlag')
    medialToLgnFlag = cfg.medialToLgnFlag;
else
    medialToLgnFlag= 1;
end

if isfield(cfg,'plotFlag')
    plotFlag = cfg.plotFlag;
else
    plotFlag = 0;
end


fgOriginal = fg;

%(1) Cut the fibers at percentLength
fg = getFgBetweenCoords(fg, percentLength); % Don't worry. it can take fibers starting either at LGN or V1, and understand where LGN is

%(2) Take only that part of the fibers which is posterior to the LGN
nodesNum = 1;
for fI = 1:length(fg.fibers)
    nodesNum = max(nodesNum, length(fg.fibers{fI}));
end

y = nan(length(fg.fibers),nodesNum);
x = nan(length(fg.fibers),nodesNum);

for fI = 1:length(fg.fibers)
    fInodesNum = size(fg.fibers{fI},2);
    if abs(fg.fibers{fI}(2,1)) < abs(fg.fibers{fI}(2,end)) % If fibers coordinates start at the LGN and end at V1
        x(fI,1:fInodesNum) = fg.fibers{fI}(1,:);
        y(fI,1:fInodesNum) = fg.fibers{fI}(2,:);
    else
        x(fI,1:fInodesNum) = fg.fibers{fI}(1,end:-1:1);
        y(fI,1:fInodesNum) = fg.fibers{fI}(2,end:-1:1);
    end
end
y1Median = nanmedian(y(:,1)); % the median first y coordinate of all fibers (an approximation to the LGN)
yThreshold = y1Median - 5; % This used to be 10, but I lowered it to 2. Needed to ignore coordinates at the LGN itself
x1Median = nanmedian(abs(x(:,1))); % the median x coordinate of all fibers (an approximation to the LGN)
x1Std= nanstd(abs(x(:,1))); % the median x coordinate of all fibers (an approximation to the LGN)

for fI = 1:length(fg.fibers)
    fInodesNum = size(fg.fibers{fI},2);
    
    if abs(fg.fibers{fI}(2,1)) < abs(fg.fibers{fI}(2,end)) % If fibers coordinates start at the LGN and end at V1
        y(fI,1:fInodesNum) = fg.fibers{fI}(2,:);
    else
        y(fI,1:fInodesNum) = fg.fibers{fI}(2,end:-1:1);
    end
end

for fI = 1:length(fg.fibers)
    indicesToEliminate = find(fg.fibers{fI}(2,:) > yThreshold);
    fg.fibers{fI}(:,indicesToEliminate) = [];
end


flag = 1; % a flag to continue the cleaning procedure

for k = 1:iterNum
    
    if flag == 0
        disp(['Continuing on iteration: ' num2str(k)])
        continue
    end
    
    % Get maximal number of nodes to save (each fiber has a different
    % number of nodes, so to put them all in a matrix, it must be the size
    % of the maximal number of nodes, "nodesNum")
    nodesNum = 1;
    for fI = 1:length(fg.fibers)
        nodesNum = max(nodesNum, length(fg.fibers{fI}));
    end
    
    % Initialize x,y,z coordinates of each node in each fiber as NaN
    x = nan(length(fg.fibers),nodesNum);
    y = nan(length(fg.fibers),nodesNum);
    z = nan(length(fg.fibers),nodesNum);
    
    for fI = 1:length(fg.fibers)
        if isempty(fg.fibers{fI}) % ROEY: added this to avoid error when fiber{fI} is empty
            continue
        end
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
    
    if nanmean(x(:)) > 0
        rightHemisphere = 1;
    else
        rightHemisphere = 0;
    end
    
    % Create a mask on the coordinates, to only look at coordinates
    % posterior to the LGN (otherwise the Meyer's loop segment might mess
    % things up in the clustering)
% % %     mask = ones(size(x));
% % %     mask =   y>=(repmat(y(:,1),1,nodesNum)-20); % mask of coordinates to *ignore*, about 2 cm posterior to LGN, and anteriorly
    xMasked = abs(x);
% % %     xMasked(mask) = nan;
    
    %% Get xMin (minimal x coordinate) of each fiber, according to the wanted method
    if strcmp(method, 'xmin')
        [xMin, xMinIndices] = nanmin(xMasked'); % y<y(1) is used to only get x values between thalamus and cortex, and allow fibers to stray in x in Meyer's loop
    elseif strcmp(method, 'xmin_highstd') || strcmp(method, 'xmin_highinterquartilerange') || strcmp(method, 'xmin_highprctile')
        % Divide the length of Y to 10 segments, and use the one with
        % greatest x spread (greatest x STD, greatest x 1-99 percentile etc)
        yMasked = y;
% % %         yMasked(mask) = nan;
        
        ymin = min(nanmin(yMasked));
        ymax = max(nanmax(yMasked));
        ydist = ymax - ymin;
        
        ySegNum = 6;
        continueFlag = 1;
        
        while ySegNum>=1 & continueFlag==1 % 6 should be good for most subjects, but if it fails, try reducing the number of segments "ySegNum"
            yseg = [ymin + [0:ySegNum-1]*ydist./ySegNum; ymin + [1:ySegNum]*ydist./ySegNum]; % each column is the two borders of its y segment
            
            xStd = 0;
            
            for i = 1:size(yseg,2)
                segMask = zeros(size(yMasked));
                segMask(y>=yseg(1,i) & y<=yseg(2,i)) = 1; % mask of coordinates to *use*, so ignore all points with y outside the segment
                
                xInSeg = abs(x);
                xInSeg(~logical(segMask)) = nan;
                if all(all(isnan(xInSeg')))
                % % %                 if any(all(isnan(xInSeg')))
%                     error('There are fibers that have no coordinates in this Y-segment. Consider using larger segments (e.g. lower ySegNum = 5 (instead of 6)), or call fgGetLateralIndices with a larger cfg.percentLength (e.g. 0.6)');
% % %                     ySegNum = ySegNum-1;
% % %                     disp(['Reducing ySegNum to ' num2str(ySegNum)]);
% % %                     break
                    xInSeg = 0;
                end
                
                if strcmp(method, 'xmin_highstd')
                    xStdInSeg = nanstd(xInSeg(:));
                elseif strcmp(method, 'xmin_highinterquartilerange')
                    xStdInSeg = iqr(xInSeg(:));
                elseif strcmp(method, 'xmin_highprctile')
                    xStdInSeg = diff(prctile(xInSeg(:),[1 99]));
                end
                
                if xStdInSeg > xStd % if this segment has higher STD in x coordinates than in previous y segments
                    [xMin, xMinIndices] = nanmin(xInSeg');
                    
                    % Save y of xMin for plot
                    yOfXMin = zeros(size(xMin));
                    for t = 1:length(yOfXMin)
                        yOfXMin(t) = y(t,xMinIndices(t));
                    end
                    
                    xStd = xStdInSeg;
                else
                    continue
                end
                
            end
            continueFlag = 0; % Stop the while loop and continue
        end
        if continueFlag == 1
            error('Could not find any y segment with fibers in it');
        end
        
        
        %%% Start of new method
        elseif strcmp(method, 'xmin_correlation')
          % Divide the length of Y to 10 segments, and use the one with
        % greatest x spread (greatest x STD)
        yMasked = y;
        
        ymin = min(nanmin(yMasked));
        ymax = max(nanmax(yMasked));
        ydist = ymax - ymin;
        yseg = [ymin + [0:9]*ydist./10; ymin + [1:10]*ydist./10]; % each column is the two borders of its y segment
        
        xMinMat = zeros(size(yseg,2), size(abs(x),1));
        for i = 1:size(yseg,2)
            segMask = zeros(size(yMasked));
            segMask(y>=yseg(1,i) & y<=yseg(2,i)) = 1; % mask of coordinates to *use*, so ignore all points with y outside the segment
            
            
            xInSeg = abs(x);
            xInSeg(~logical(segMask)) = nan;
            
            [xMin, xMinIndices] = nanmin(xInSeg');
            xMinMat(i,:) = xMin;
                        
        end
        %%% End of new method
            
            
    end
    
    
    %% Cluster
    
    switch clusterMethod
        case 'histogramClustering'
            % this usually works, but in some cases some tweak is needed, which is why there is another option
            xMinHist = histogram(xMin,'BinWidth',0.2);
            values = xMinHist.Values;
            values(values<=3) = 0; % This is a threshold on the histogram values - we ignore X coordinated with only 3 fibers, because we want to look only at robust peaks in the histogram
            % A SECOND VERSION USES PERCENTILE INSTEAD OF AN ABSOLUTE THRESHOLD: values(values<=prctile(values(values>0),10)) = 0; % This is a threshold on the histogram values - we ignore X coordinated with only 3 fibers, because we want to look only at robust peaks in the histogram
            [~, maxIdx] = max([zeros(1,ceil(length(values)/2-1)), values(ceil(length(values)/2):end)]); % I use this instead of max(values) to fix cases in which the medial fiber is actually where the maximum is obtained, like in subject 101. This is like looking for the point of maximal number of fibers, but only in the lateral half of the possible coordinates
            binaryValues = istrue(values);
            binaryValues = binaryValues(1:maxIdx); % CTruncate the histogram to bin "maxIdx"
            onesIndices = find(binaryValues == 1); % fill in first sequence of zeros with ones
            if isempty(onesIndices) % Trying to fix the case of subject 11, who gets an empty onesIndices vector
                xMinHist = histogram(xMin,'BinWidth',0.2);
                values = xMinHist.Values;
                values(values<=3) = 0; % This is a threshold on the histogram values - we ignore X coordinated with only 3 fibers, because we want to look only at robust peaks in the histogram
                % Remove most lateral zeros (this fixes cases where one might get an empty onesIndices vector later, like subject 11)
                lastNonZero = find(values>0);
                if isempty(lastNonZero)
                   disp('Too few fibers.');
                   lateralIndices = [];
                   return 
                end
                lastNonZero = lastNonZero(end);
                if length(values)>lastNonZero
                    values(lastNonZero+1:end) = [];
                end
                [~, maxIdx] = max([zeros(1,ceil(length(values)/2-1)), values(ceil(length(values)/2):end)]); % I use this instead of max(values) to fix cases in which the medial fiber is actually where the maximum is obtained, like in subject 101. This is like looking for the point of maximal number of fibers, but only in the lateral half of the possible coordinates
                binaryValues = istrue(values);
                binaryValues = binaryValues(1:maxIdx);
                onesIndices = find(binaryValues == 1); % fill in first sequence of zeros with ones
            end
            
            if isempty(onesIndices)
                error('The histogramClustering method failed. I should think of a solution for this')
            end
            binaryValues(1:onesIndices(1)) = 1;
            
            diffSig = diff(binaryValues);
            
            if all(diffSig == 0)
                Cbest = [minDist+1 0];
                idxbest = ones(size(xMin));
            else
                
                startIndex = find(diffSig < 0); % Where segments of 0's start
                endIndex = find(diffSig > 0)-1; % Where segments of 0's end
                duration = endIndex-startIndex+1; % Length of the different 0's segments
                
                [~ ,maxDurationIdx] = max(duration);
%                 binIdx = ceil( (startIndex(maxDurationIdx)+endIndex(maxDurationIdx))/2); % Middle of the longest 0's segment. WHY AM I NOT USING JUST "endIndex(maxDurationIdx)"?
                binIdx = endIndex(maxDurationIdx); % Lateral end of the longest 0's segment.
                % TODO: maybe change here to "binIdx = endIndex(maxDurationIdx)" or startIndex, in case it's the other hemisphere
                threshold = xMinHist.BinEdges(binIdx); % threshold on the x coordinate
                
                idxbest = ones(size(xMin));
                idxbest(xMin<threshold) = 2;
                
                Cbest = [];
                Cbest(1) = nanmedian(xMin(idxbest==1)); % Representative median value of one bundle
                Cbest(2) = nanmedian(xMin(idxbest==2)); % Representative median value of the other bundle
            end
        case 'kmedoids'
            disp(['Clustering using method ' clusterMethod]);
            [idxbest, Cbest, ~, ~] = kmedoids(xMin', 2,'Algorithm','small');
        case 'kmeans'
            disp(['Clustering using method ' clusterMethod]);

            [idxbest, Cbest, ~, ~] = kmeans(xMin', 2,'Replicates',1000);
        case 'gmm'
            disp(['Clusterin using method ' clusterMethod]);

            gmfit = fitgmdist(xMin', 2);
            idxbest = cluster(gmfit,xMin');
            Cbest = gmfit.mu';
       
            if Cbest(1)<Cbest(2) % if the centroid has smaller x
                lateralIdx = 2;
                medialIdx = 1;
            else
                lateralIdx  = 1;
                medialIdx = 2;
            end
            P = posterior(gmfit,xMin');
            indicesOfLowProbMedial = find(idxbest==medialIdx & P(:,medialIdx)<0.95);
            idxbest(indicesOfLowProbMedial) = lateralIdx;
        case 'kmedoids_correlation'
            disp(['Clustering using method ' clusterMethod]);
            [idxbest, Cbest, ~, ~] = kmedoids(xMinMat', 2,'Distance','correlation');
            
        otherwise
            disp('No clustering method chosen');
    end
    Cbest = abs(Cbest);
    
    
    disp(['Distance between bundles: ' num2str(abs(Cbest(1)-Cbest(2)))])
    
    if Cbest(1)<Cbest(2) % the smaller one (in absolute value) is the medial one
        lateralIdx = 2;
        medialIdx = 1;
    else
        lateralIdx  = 1;
        medialIdx = 2;
    end
    
    if medialToLgnFlag % All fibers lateral to LGN in the coordinates used for classification will be treated as lateral
%         idxbest(xMin > x1Median+x1Std) = lateralIdx; % xMin is in absolute value, so this is good for both hemispheres
        
        % Recompute Cbest
        switch clusterMethod
            case 'kmedoids'
                Cbest(medialIdx) = nanmedian(xMin(idxbest==medialIdx));
                Cbest(lateralIdx) = nanmedian(xMin(idxbest==lateralIdx));
                disp(['Corrected distance between bundles: ' num2str(abs(Cbest(1)-Cbest(2)))])
            case 'kmeans'
                Cbest(medialIdx) = nanmean(xMin(idxbest==medialIdx));
                Cbest(lateralIdx) = nanmean(xMin(idxbest==lateralIdx));
                disp(['Corrected distance between bundles: ' num2str(abs(Cbest(1)-Cbest(2)))])
            case 'histogramClustering'
                Cbest(medialIdx) = nanmedian(xMin(idxbest==medialIdx));
                Cbest(lateralIdx) = nanmedian(xMin(idxbest==lateralIdx));
                disp(['Corrected distance between bundles: ' num2str(abs(Cbest(1)-Cbest(2)))])
   
            otherwise
                disp([clusterMethod ' not implemented for fixing clustering by using the medialToLgnFlag']);
        end
    end
    
    lateralIndices = 1:length(fg.fibers);
    lateralIndices (idxbest == medialIdx) = [];
    
    
    fg_medial = fgRetainIndices(fg, idxbest == medialIdx);
    fg_lateral = fgRetainIndices(fg, lateralIndices);
    
    %% Plot
    if plotFlag
        % find the y of each xmin, to see where the fibers were
        % separated, and plot
        yOfXMin = zeros(size(xMin));
        for t = 1:length(yOfXMin)
            yOfXMin(t) = y(t,xMinIndices(t));
        end
        AFQ_RenderFibers(fg_medial ,'color',[1 1 0], 'camera', [0 90], 'numfibers',length(fg_medial.fibers), 'tubes', [0], 'jittershading', 0.5);
        AFQ_RenderFibers(fg_lateral,'color',[0 1 1], 'camera', [0 90], 'numfibers',length(fg_lateral.fibers), 'tubes', [0], 'jittershading', 0.5,'newfig',0);
        cmap = lines(3);
        AFQ_RenderFibers(fg_medial ,'color',cmap(3,:), 'camera', [180 -90], 'numfibers',length(fg_medial.fibers), 'tubes', [0], 'jittershading', 0.2,'newfig',0);
        AFQ_RenderFibers(fg_lateral,'color',cmap(3,:), 'camera', [180 -90], 'numfibers',length(fg_lateral.fibers), 'tubes', [0], 'jittershading', 0.2,'newfig',0);
        
        view([180 -90])
        zlim = get(gca,'zlim');
        set(gca,'zlim',[-40,zlim(2)])
        if rightHemisphere
            scatter3(xMin(find(idxbest==1)), yOfXMin(find(idxbest==1)),-40*ones(size(yOfXMin(find(idxbest==1)))))
            scatter3(xMin(find(idxbest==2)), yOfXMin(find(idxbest==2)),-40*ones(size(yOfXMin(find(idxbest==2)))))
        else
            scatter3(-xMin(find(idxbest==1)), yOfXMin(find(idxbest==1)),-40*ones(size(yOfXMin(find(idxbest==1)))))
            scatter3(-xMin(find(idxbest==2)), yOfXMin(find(idxbest==2)),-40*ones(size(yOfXMin(find(idxbest==2)))))
        end
        title( ['Iteration: ' num2str(k), ', Distance: ' num2str(abs(Cbest(1)-Cbest(2))), ' ' clusterMethod] );
    end
    
    %% Check if bundles are well separated in space
    
    if abs(Cbest(1)-Cbest(2))>minDist % if the two bundles are really separated in space, by at least minDist mm in their centroids
        
        
        % Remove the coordinates of medial fibers
        medialIndices = find(idxbest == medialIdx);
        for m = 1:length(medialIndices )
            fg.fibers{medialIndices(m)} = [nan nan nan]';
        end
     
%%%%%%%% This was a fix used to make some medial fibers to lateral, under
%%%%%%%% some conditions like going anteriorily, but that's just a bad
%%%%%%%% patch.
% % % % %         for fI = medialIndices'
% % % % %             % If a fiber has y coordinates that extend anterioraly to the LGN more than 5 mm, consider
% % % % %             % it a lateral fiber, and not a medial fiber
% % % % %             if any(fg.fibers{fI}(2,:) > y1Median + 5)
% % % % %                 continue
% % % % %             end
% % % % %             
% % % % %             % % % %% I decided against this condition, which gave bad
% % % % %             % results, for example in subject 4
% % % % %             % % %             % If a fiber has z coordinates that start by going up, and not down, consider it a lateral fiber, and not a medial fiber
% % % % %             % % %             if sum(diff((fg.fibers{fI}(3,1:5)))) > 0
% % % % %             % % %                 continue
% % % % %             % % %             end
% % % % %             
% % % % %             fg.fibers{fI} = [nan nan nan]';
% % % % %         
% % % % %         
% % % % %         end
        
        continue
        
    else
        flag = 0;
        continue
    end
    
end

% AND THIS IS THE PART USED WHEN THE PATCH ABOVE WAS USED
% % % % lateralIndices = 1:length(fg.fibers);
% % % % excludeIndices = [];
% % % % for fI = 1:length(fg.fibers)
% % % %     if isnan(fg.fibers{fI}(1))
% % % %         excludeIndices(end+1) = fI;
% % % %     end
% % % % end
% % % % lateralIndices(excludeIndices) = [];

lateralIndices = 1:length(fg.fibers); % by default, all fibers are considered lateral
medialIndices = find(cellfun(@(x) all(isnan(x(:))), fg.fibers));
lateralIndices(medialIndices) = [];
