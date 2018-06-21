%% Before you start, your subject should have the following:
% * The subject's output folder of dtiInit.m (e.g., dt64trilin)
% * The subject's output folder of FreeSurfer
% * The subject's output folder of mrQ

%% (1) Create the ROIs ndecessary for tractography and warp them to b0
fsDir = ''; % Subject's FreeSurfer directory
T1wFile = ''; % Subject's T1w file (preferably the one used to run FreeSurfer)
outDir = ''; % Output directory where the nifti ROIs will be written
extract_rois(fsDir, T1wFile, outDir);

roisFiles = {'','','',''}; % The ROIs files from the previous step
dtDir = ''; % The directory in which the subject's dt6.mat file is saved
dilationSteps = 2;
warp_rois_to_b0(roisFiles, dtDir, T1wFile, dilationSteps);


%% (2) Run ConTrack tractography
% Note: When running the ConTrack tractography, it is advised not to use
% the default wmProb.nii.gz (output of dtiInit.m) as the propagation mask,
% but rather use a binary WM mask. This can be the wmProb.nii.gz
% thresholded at some value (e.g., 0.65), a WM mask from FreeSurfer, warped
% and resampled to diffusion space, etc.


%% (3) Create an OR benchmark
fgFile = ''; % fgFile can be the .pdb output file of ConTrack, or a converted version saved in .mat format
fg = fgRead(fgFile);
extract_lgn_from_thalamus(thalamusFile,leftOrRight,lgnOutputFile);
fg = get_or_benchmark(fg,dtDir,lgnRoi); % This will add the relevant params fields to extract the OR benchmark
fgBenchmark = fgRetainIndices(fg, fgGetParams(fg,'lgn_not_cc_highz_lateral'));

%% (4) Add the T1 statistics
T1File = ''; % The quantitative T1 file
t1_img = readFileNifti(T1File);

% Mask the images
useWmMask = 0;
if useWmMask
    wmMaskFile = ''; % The WM mask file in nifti format, in b0 space
    wmMask = readFileNifti(wmMaskFile);
    t1_img.data(~wmMask.data) = nan;
end
perPointFlag = 0;
fg = dtiCreateQuenchStats(fg, 'T1_std', 'T1', perPointFlag, t1_img, 'nanstd', 1);
fg = dtiCreateQuenchStats(fg, 'T1_median', 'T1', perPointFlag, t1_img, 'nanmedian', 1);

% NOTE ==== if dtiCreateQuenchStats doesn't include a 'nanstd' or
% 'nanmedian' option, please add the following lines under the
% "switch(lower(fiberSummary))" in dtiCreateQuenchStats:
%     case 'nanstd'
%         zeroIndices = cellfun(@find,cellfun(@(x)
%         istrue(x==0),vals,'UniformOutput', false),'UniformOutput',
%         false); for i = 1:length(zeroIndices)
%             if isempty(zeroIndices{i})
%                 continue
%             end vals{i}(zeroIndices{i}) = nan;
%         end param.stat = cellfun(@nanstd,vals);
%     case 'nanmedian'
%         zeroIndices = cellfun(@find,cellfun(@(x)
%         istrue(x==0),vals,'UniformOutput', false),'UniformOutput',
%         false); for i = 1:length(zeroIndices)
%             if isempty(zeroIndices{i})
%                 continue
%             end vals{i}(zeroIndices{i}) = nan;
%         end param.stat = cellfun(@nanmedian,vals);

%% (5) Filter by T1-STD and T1_Mdn
% disp('Writing fiber density maps for T1wmThresh_std & T1wmThresh_median')
% fdIndicesWrite2(fgOrL, 'T1wmThresh_std', prctileVec1,
% 'T1wmThresh_median', prctileVec2, size(t1w.data), t1w.qto_xyz,
% fullfile(fdImgsDir,'fdIndices_t1wmThreshstdt1wmThreshmedian_L_moreCombinedThresholds'),false,false,'left')
% fdIndicesWrite2(fgOrR, 'T1wmThresh_std', prctileVec1,
% 'T1wmThresh_median', prctileVec2, size(t1w.data), t1w.qto_xyz,
% fullfile(fdImgsDir,'fdIndices_t1wmThreshstdt1wmThreshmedian_R_moreCombinedThresholds'),false,false,'right')

t1StdPrctileThresh = 35; % Or any other threshold
t1MdnPrctileThresh = 67; % Or any other threshold

vals = fgGetParams(fgR, 'T1_std');% T1wmThresh_std    T1_std
thrsh = prctile(vals, t1StdPrctileThresh);
indices1 = find(vals<thrsh);
vals = fgGetParams(fgR, 'T1_median');% T1wmThresh_median    T1_median
thrsh = prctile(vals, t1MdnPrctileThresh);
indices2 = find(vals<thrsh);
indices = intersect(indices1, indices2);

fgT1filtered = fgRetainIndices(fg,indices);

%% (6) Calculate the TPR and FPR of the filtered fiber group
leftOrRight = ''; % 'left' or 'right' hemisphere
xform = t1_img.qto_xyz;
sz = size(t1_img.data);
normalizeMap = false;

% (6.1) Find indices of voxels in the T1-filtered fg
fdImg = dtiComputeFiberDensityNoGUI(fgTmp, xform, sz, normalizeMap);
if strcmp(leftOrRight,'left') % Delete voxels in the other hemisphere, in case any survived
    fdImg(floor(size(fdImg,1)/2):end,:,:) = 0;
elseif strcmp(leftOrRight,'right')
    fdImg(1:floor(size(fdImg,1)/2),:,:) = 0;
else
    error('"leftOrRight" must be "left" or "right"')
end
fiberIndices_t1Filtered = find(fdImg);

% (6.2) Find indices of voxels in the candidate fg
fdImg = dtiComputeFiberDensityNoGUI(fg, xform, sz, normalizeMap);
if strcmp(leftOrRight,'left') % Delete voxels in the other hemisphere, in case any survived
    fdImg(floor(size(fdImg,1)/2):end,:,:) = 0;
elseif strcmp(leftOrRight,'right')
    fdImg(1:floor(size(fdImg,1)/2),:,:) = 0;
else
    error('"leftOrRight" must be "left" or "right"')
end
fiberIndices_candidates = find(fdImg);

% (6.3) Find indices of voxels in the benchmark fg
fdImg = dtiComputeFiberDensityNoGUI(fgBenchmark, xform, sz, normalizeMap);
if strcmp(leftOrRight,'left') % Delete voxels in the other hemisphere, in case any survived
    fdImg(floor(size(fdImg,1)/2):end,:,:) = 0;
elseif strcmp(leftOrRight,'right')
    fdImg(1:floor(size(fdImg,1)/2),:,:) = 0;
else
    error('"leftOrRight" must be "left" or "right"')
end
fiberIndices_benchmark = find(fdImg);


% Sensitivity (TPR, true positive rate)
true_positive =  sum(ismember(fiberIndices_t1Filtered,fiberIndices_benchmark));
tpr = true_positive/length(fiberIndices_benchmark); % proportion of positives correctly identified as such

% Specificity (FPR, false positive rate)
false_positive =  sum(~ismember(fiberIndices_t1Filtered,fiberIndices_benchmark));
true_negative = length(fiberIndices_candidates) - length(fiberIndices_benchmark);
fpr = false_positive/true_negative;

% Youden's Index
sensitivity = tpr(cI,dI);
specificity = 1 - fpr(cI,dI);
w = 0;
wYouden = (1-w)*sensitivity + (1+w)*specificity - 1; % weighted Youden's index