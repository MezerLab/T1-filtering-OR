function warp_rois_to_b0(roisFiles, dtDir, T1wFile, dilationSteps)
% warp_rois_to_b0 warps the ROIs extracted by extract_rois.m to the
% diffusion space of the b0 image.
%
% INPUT
% =====
% roisFiles: Cell array of ROIs to warp
% dtDir: the directory in which the subject's dt6.mat file is saved
% T1wFile: A strcutural MRI file at the same space as that used to run
%          FreeSurferThe (T1w or a quantitative T1 map)
% dilationSteps: Number of dilation steps to perform for the ROIs before
%                warping to b0 (recommended: 2).
%
% OUTPUT
% ======
% The Nifti of each ROI will be saved at the same directory as roisFiles,
% but with a _2DTI suffix.
%
% Based on the mrQ_registerMap2DTI function of the mrQ toolbox.


%%
b0File = fullfile(dtDir,'bin','b0.nii.gz');

%% Calculate the transform from T1w to b0
% Set the name for the output of the function
out = fullfile(dtDir,'bin','map2B0');

warpFile = [out 'Warp.nii.gz'];
affineFile = [out 'Affine.nii.gz'];

if ~exist(warpFile,'file') || ~exist(affineFile,'file')
    % make a trasformation from T1wFile (moving image) to Bofile (fixed image)
    cmANTS=['xterm -e ANTS 3 -m MI[' b0File ',' T1wFile ',1,2] -o ' out '.nii.gz --rigid-affine true'];
    
    % Run the command in unix and get back status and results:
    [status, ~] = system(cmANTS);
    
    % If the above command fails...
    if status ~= 0
        cmANTS=['ANTS 3 -m MI[' b0File ',' T1wFile ',1,2] -o ' out '.nii.gz --rigid-affine true'];
        % Run the command in unix and get back status and results:
        [~, ~] = system(cmANTS,'-echo');
    end
end


%% Apply the trasformation to the images
[~, name] = fileparts(T1wFile);
[~, name] = fileparts(name);

refFile = b0File; % Use b0 as reference, to move everything to diffusion space
T1wWarped = fullfile(outDir,[name '_2DTI_resamp.nii.gz']);
interpMethod = '--use-NN'; % Use nearest-neighbor interpolation

cmWarp=['xterm -e WarpImageMultiTransform  3 ' T1wFile  ' ' T1wWarped ' -R ' refFile ' ' out 'Warp.nii.gz ' out 'Affine.txt ' interpMethod];
% Run the command in unix and get back status and results:
[status, ~] = system(cmWarp);

% If the command fails, try without opening a new terminal (i.e., without xterm)
if status ~= 0
    cmWarp=['WarpImageMultiTransform  3 ' T1wFile  ' ' T1wWarped ' -R ' refFile ' ' out 'Warp.nii.gz ' out 'Affine.txt ' interpMethod];
    [~, ~] = system(cmWarp,'-echo');
end

%% Dilate the ROIs before warping, if wanted
if dilationSteps>0
    for rI = 1:length(roisFiles)
        inFile = roisFiles{rI};
        roi = readFileNifi(inFile);
        
        for voxStep = 1:dilationSteps
            roi.data = roi.data + circshift(roi.data,1,1) + circshift(roi.data,1,2) + circshift(roi.data,1,3) + circshift(roi.data,-1,1) + circshift(roi.data,-1,2) + circshift(roi.data,-1,3);
        end
        roi.data = double(logical(roi.data));
        
        roisFiles{rI} = [roisFiles{rI}(1:end-7) '_dilated.nii.gz'];
        dtiWriteNiftiWrapper(roi.data, roi.qto_xyz, roisFiles{rI});
    end
end

%% Warp the ROIs
for rI = 1:length(roisFiles)
    inFile = roisFiles{rI};
    outFile = inFile;
    outFile = [outFile(1:end-7), '_2b0.nii.gz.'];
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' inFile  ' ' outFile ' -R ' b0File  ' ' warpFile ' ' affineFile  ' ' interpMethod];
    [~, ~] = system(cmWarp);
end

%% If temporary dilated files were created, delete them
if dilationSteps>0
    delete(roisFiles{rI});
end