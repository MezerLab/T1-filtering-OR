function extract_rois(fsDir, refFile, outDir)
% extract_rois extracts V1 and Thalamus from the aparc+aseg parcellation.
% 
% INPUT
% =====
% fsDir: The full path of subject's FreeSurfer directory
% refImg: A Nifti of the same space and resolution as the one used for
%         running FreeSurfer (this is used for converting the
%         aparc+aseg.mgz file into Nifti format);
% outDir: The full path of the output directory where the ROI Nifti file
%         will be saved
%
% OUTPUT
% ======
% The Nifti of each ROI will be saved under outDir, for both hemispheres:
% L_V1.nii.gz and R_V1.nii.gz, L_Thalamus.nii.gz and R_Thalamus.nii.gz.
%

%% Convert the aparc+aseg.mgz file to .nii.gz format if necessary
segFile = fullfile(fsDir, 'mri/aparc+aseg.nii.gz');
if ~exist(segFile,'file')
    mgzIn  = fullfile(fsDir, 'mri/aparc+aseg.mgz');
    niiOut = segFile;
    if exist(refFile,'file')
        orientation = 'RAS';
        fs_mgzSegToNifti(mgzIn, refFile, niiOut, orientation);
    else
        error('Cannot convert aparc+aseg.mgz to nii.gz format. Missing reference structural image (refImg) for conversion.');
    end
end
seg = readFileNifti(segFile);

%% Extract V1 as calcarine fissure
fs_aparcAsegLabelToNiftiRoi(seg,1021,fullfile(outDir,'L_V1.nii.gz'));
fs_aparcAsegLabelToNiftiRoi(seg,2021,fullfile(outDir,'R_V1.nii.gz'));

%% Extract thalamus for tractography and/or for finding the LGN
fs_aparcAsegLabelToNiftiRoi(seg,10,fullfile(outDir,'L_Thalamus.nii.gz'));
fs_aparcAsegLabelToNiftiRoi(seg,49,fullfile(outDir,'R_Thalamus.nii.gz'));
end