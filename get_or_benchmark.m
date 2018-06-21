function fg = get_or_benchmark(fg,dtDir,lgnRoi)
% or_get_benchmark takes a candidate fiber group of streamline
% connecting the thalamus and area V1, and returns a cleaned subset of
% streamlines representing the optic radiation.
%
% INPUT
% =====
% fg: a fiber group in vistasoft format (usually saved in a .mat file),
%     including all optic radiation candidates (see fgCreate)
% dtDir: the directory in which the subject's dt6.mat file is saved
% lgnRoi: ROI representing the LGN in vistasoft format (see roiCreate)
% 
% OUTPUT
% ======
% The same fg as the input, with added cells in the params field. These 
% cells contain the binary vectors that allow filtering the streamlines.
% For example, to get the benchmark fg from the output use:
% fgBenchmark = fgRetainIndices(fg, fgGetParams(fg,'lgn_not_cc_highz_lateral'));

%% (1) Find all fibers crossing the LGN
minDist = 1;
fgStartEnd = fg; % Keep only start and end nodes indices
fgStartEnd.fibers = cellfun(@(x) x(:,[1,end]), fgStartEnd.fibers,'UniformOutput',false);
[~,~, keepIndices] = dtiIntersectFibersWithRoi([], 'not', minDist, lgnRoi, fgStartEnd); % keepIndices are really those that intersect with the ROI, even when the 'not' option is used
% Create binary indices vector in a params field
fg.params{end+1}.name = 'is_in_lgn';
fg.params{end}.stat = keepIndices;

%% (2) Eliminate fascicles passing through the Corpus Callosum
ccRoiFile = fullfile(dtDir,'ROIs','callosum_rough.mat');
if ~exist(ccRoiFile,'file')
    % This piece of code was adapted from AFQ's AFQ_SegmentCallosum.m
    % Load Dt6
    dt = dtiLoadDt6(fullfile(dtDir,'dt6.mat'));
    % Create an ROI of the corpus callosum
    [~,~,ccRoi]=dtiCreateRoiFromMniNifti(dt.dataFile, fullfile(AFQ_directories,'templates','callosum','callosum_rough.nii.gz'));
    % Load wholebrain fiber group
    wholebrainFg = AFQ_get(afq,'wholebrain fg', ii);
    % Create callosum fiber group
    ccFg = dtiIntersectFibersWithRoi([],'and',2,ccRoi,wholebrainFg);
    ccFg.name = 'callosumFG';
    % Save callosum ROI and fiber group
    fprintf('\nSaving %s',fgPath)
    dtiWriteFiberGroup(ccFg,fgPath);
    dtiWriteRoi(ccRoi,ccRoiFile);
end

minDist = 2;

ccRoi = dtiReadRoi(ccRoiFile);
% Expand the ROI in the Z direction to ensure the CC is fully covered
ccRoi.coords = unique([ccRoi.coords; ccRoi.coords+repmat([0 0 1],size(ccRoi.coords,1),1); ccRoi.coords+repmat([0 0 2],size(ccRoi.coords,1),1); ccRoi.coords+repmat([0 0 3],size(ccRoi.coords,1),1); ccRoi.coords+repmat([0 0 4],size(ccRoi.coords,1),1); ccRoi.coords+repmat([0 0 -1],size(ccRoi.coords,1),1); ccRoi.coords+repmat([0 0 -2],size(ccRoi.coords,1),1); ccRoi.coords+repmat([0 0 -3],size(ccRoi.coords,1),1); ccRoi.coords+repmat([0 0 -4],size(ccRoi.coords,1),1)],'rows');

[~,~, notkeep] = dtiIntersectFibersWithRoi([], 'not', minDist, ccRoi, fg); % 'not' means return only fascicles that don't pass through the ROI
keepIndices = ~notkeep;
% Create binary indices vector in a params field
fg.params{end+1}.name = 'is_not_in_cc';
fg.params{end}.stat = keepIndices;

%% (3) Eliminate fascicles going down through the brainstem
cfg = [];
cfg.method = 'minz'; % or 'meanz'
cfg.minDist = 15;
cfg.percentLength = 0.5;

keepIndices = find(fgGetParams(fg, 'is_in_lgn') & fgGetParams(fg, 'is_not_in_cc')); % "find" is used to get the indices themselves, not a binary vector
fgTemp = fgetainIndices(fg, keepIndices);
[indicesLowz, zDist] = fgGetLowZIndices(cfg, fgTemp); % bad fibers (too low, go through Cerebellum/brainstem)
tempIndices = 1:length(fgTemp.fibers);
tempIndices  = tempIndices (~ismember(tempIndices , indicesLowz)); % Take only good fibers, not those that go through Cerebellum
% Create binary indices vector in a params field
highzIndices = zeros(1,length(fg.fibers));
highzIndices(keepIndices(tempIndices)) = 1;
fg.params{end+1}.name = 'lgn_not_cc_highz';
fg.params{end}.stat = highzIndices;

%% (4) Cluster by x coordinate and separate lateral bundle (the OR) from medial bundles (probably Splenium, not pulvinar bundles)
% Get a new fg with only some of the path (between LGN and V1, with
% going anteriorily)
cfg = [];
cfg.method = 'xmin_highprctile';            %'xmin_highprctile'; % xmin_highstd OR xmin OR xmin_highinterquartilerange OR xmin_highprctile
cfg.clusterMethod = 'histogramClustering'; % kmedoids
cfg.iterNum = 3;          % number of iterations (default 1)
cfg.minDist = 5;         % the minimal distance to separate medial and lateral bundles
cfg.percentLength = 0.5; % percent of fiber length to use, between its most anterior and most posterior coordinates (default 0.5)
cfg.plotFlag = 1;         % 1 ('yes') or 0 ('no') (default 0)

keepIndices = find(fgGetParams(fg, 'lgn_not_cc_highz')); % "find" is used to get the indices themselves, not a binary vector
fgTemp = fgRetainIndices(fg, keepIndices);
indicesLateral = fgGetLateralIndices(cfg, fgTemp);
indicesMedial = 1:length(fgTemp.fibers);
indicesMedial = indicesMedial(~ismember(indicesMedial, indicesLateral));
keepIndices(indicesMedial) = [];
% Create binary indices vector in a params field
lateralIndices = zeros(1,length(fg.fibers));
lateralIndices(keepIndices) = 1;
fg.params{end+1}.name = 'lgn_not_cc_highz_lateral';
fg.params{end}.stat = lateralIndices;
