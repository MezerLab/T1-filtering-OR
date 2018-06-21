function fgOut = getFgBetweenCoords(fg, percentLength)

fgOut = fg;
minYFibers = zeros(1,size(fg.fibers,2));
maxYFibers = zeros(1,size(fg.fibers,2));

for fI = 1:size(fg.fibers,1)
    minYFibers(fI) = min(fg.fibers{fI}(2,:));
    maxYFibers(fI) = max(fg.fibers{fI}(2,:));
end

medianYmin = median(minYFibers);
medianYmax = median(maxYFibers);


thresholdY = medianYmax - percentLength*(medianYmax-medianYmin);


for fI = 1:size(fg.fibers,1)
    if fg.fibers{fI}(2,1) < fg.fibers{fI}(2,end) % if the first coordinate is more posterior
        fg.fibers{fI} = fg.fibers{fI}(:,end:-1:1); % make all fibers start at the thalamus
    end
    
    index = find(fg.fibers{fI}(2,:)<thresholdY);
    if isempty(index)
        continue
    end
    index = index(1);
    fgOut.fibers{fI} = fg.fibers{fI}(:,1:index);
end



% nodesNum = size(coords,2);
% rgb = vals2colormap(1:nodesNum);
% 
% scatter3(coords(1,:), coords(2,:), coords(3,:),[],rgb,'filled');
% pause(0.5)
% drawnow
