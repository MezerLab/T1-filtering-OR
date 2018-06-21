function extract_lgn_from_thalamus(thalamusFile,leftOrRight,lgnOutputFile)
% extract_lgn_from_thalamus extract an estimated LGN ROI based on a binary
% mask of the thalamus in nifti format.
%
% INPUT
% =====
% thalamusFile: A nifti file of the thalamus mask (See extract_rois.m)
% leftOrRight: a string ('left' or 'right') indicating if this is the left
%              or right thalamus
% OUTPUT
% ======
% lgnOutputFile: The LGN ROI (Vistasoft's ROI format) output file. This
%                should have a .mat suffix

roiThalamusNii = readFileNifti(thalamusFile);
roiThalamus = roiThalamusNii.data;

%% Set sphere center at the 20% most posterior, 10% most lateral voxels (of the 20% most posterior), 10% z voxel
% X,Y,Z coordinates are in *voxels* here, and later transformed to mm
% Get center Y coordinate
posteriorThreshold = 20;
lateralThreshold = 10;
zThreshold = 10;

indices = find(roiThalamus);
[x,y,z] = ind2sub(size(roiThalamus),indices);
thresholds = prctile(unique(y),[posteriorThreshold]);
[~,idxOfPrctileY] = min(abs(y-thresholds));
idxOfPrctileY = idxOfPrctileY(1);
yCenter = y(idxOfPrctileY);
% Cut ROI anterior to Y coordinate
thresholds = prctile(unique(y),[0 posteriorThreshold]);
idxOfLowY = find(y>thresholds(1) & y<thresholds(2)); % posterior
idxOfHighY = find(y>=thresholds(2)); % anterior
roiThalamus(indices(idxOfHighY)) = 0;
% Get center X coordinate
indices = find(roiThalamus);
[x,y,z] = ind2sub(size(roiThalamus),indices);
if strcmpi(leftOrRight,'L') % Left thalamus
    thresholds = prctile(unique(x),[lateralThreshold]);
    [~,idxOfPrctileX] = min(abs(x-thresholds));
    idxOfPrctileX = idxOfPrctileX(1);
    xCenter = x(idxOfPrctileX);
    
    thresholds = prctile(unique(x),[0 lateralThreshold]);
    idxOfLowX = find(x>thresholds(1) & x<thresholds(2)); % medial
    idxOfHighX = find(x>=thresholds(2)); % lateral
    % Set density to zero in lateral points
    roiThalamus(indices(idxOfHighX)) = 0; % NOTE: why am I doing this? perhaps lateral indices are too high z?
elseif strcmpi(leftOrRight,'R') % Right thalamus
    thresholds = prctile(unique(x),100-[lateralThreshold]);
    [~,idxOfPrctileX] = min(abs(x-thresholds));
    idxOfPrctileX = idxOfPrctileX(1);
    xCenter = x(idxOfPrctileX);
    
    thresholds = prctile(unique(x),100-[0 lateralThreshold]);
    idxOfLowX = find(x<thresholds(2)); % medial
    idxOfHighX = find(x>=thresholds(2)); % lateral
    % Set density to zero in lateral points
    roiThalamus(indices(idxOfLowX)) = 0;
else
    error('The leftOrRight argument must be ''L'' or ''R');
end
[x,y,z] = ind2sub(size(roiThalamus),indices);
thresholds = prctile(unique(z),[zThreshold]);
[~,idxOfPrctileZ] = min(abs(z-thresholds));
idxOfPrctileZ = idxOfPrctileZ(1);
zCenter = z(idxOfPrctileZ);

if strcmpi(leftOrRight,'L') % Left thalamus
    roiName = 'L_LGN';    
else
    roiName = 'R_LGN';
end

centerCoords = mrAnatXformCoords(roiThalamusNii.qto_xyz, [xCenter,yCenter,zCenter], true);
radiusInVoxels = 4; % 4 mm (although I think inside dtiBuildSphereCoords it says it's voxels. It's not).
roi = dtiNewRoi(['sphere_', num2str(radiusInVoxels)]);
roi.coords = dtiBuildSphereCoords(centerCoords, radiusInVoxels);

% Save LGN ROI in .mat format
dtiWriteRoi(roi, lgnOutputFile);

% Convert LGN ROI to .nii.gz format
lgnOutputFileNii = lgnOutputFile;
lgnOutputFileNii = [lgnOutputFileNii(1:end-4),'.nii.gz'];
% NOTE: I STOPPED WORKING
% WITH dtiRoiNiftiFromMat, BECAUSE IT HAS PROBLEM WITH THE
% IMAGE SIZE. E.G. IT MAKES THE bb BIGGER IN 1 IN SOME OF THE
% DIMENSIONS (THOSE WITH AN EVEN NUMBER OF VOXELS, I THINK).
% cd(orDir); dtiRoiNiftiFromMat(fullfile(orDir,[roiName
% '.mat']),b0File,roiName,1); % created with fg of
% opts.stepSizeMm = 0.2; opts.faThresh = 0.1;, after mlooking
% only at the 0-30 precentile of lowest y coordinated in the
% trackDensityMasked map
coordsImg = mrAnatXformCoords(inv(roiThalamusNii.qto_xyz), roi.coords, false);
coordsVoxels = floor(coordsImg)+1; % Note: I should check if this is the right way here.
coordsIndices = sub2ind(size(roiThalamusNii.data), coordsVoxels(1,:),coordsVoxels(2,:),coordsVoxels(3,:));
roiMask = zeros(size(roiThalamusNii.data));
roiMask(coordsIndices) = 1;
nii
dtiWriteNiftiWrapper(roiMask,roiThalamusNii.qto_xyz,lgnOutputFileNii);

end