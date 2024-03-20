function DICOM2MHA(dicomList)
%% DICOM2MHA(dicomList)
% ------------------------------------------
% FILE   : DICOM2MHA.m
% AUTHOR : Nicholas Hindley, ACRF Image X Institute, The University of Sydney
% DATE   : 2024-03-21  Created.
% ------------------------------------------
% PURPOSE
%   Convert a series of DICOM files to a .mha file in IEC geometry.
% ------------------------------------------
% INPUT
%   dicomList:      A list of dicom files in cell format (use "lscell")

%% Input check
if nargin < 1
    dicomList = lscell('CT*');
end
outputFile = 'CT.mha';
waterAtt = 0.013;

%% Read DICOM files
for k = 1:length(dicomList)
    dcmHeader{k} = dicominfo(dicomList{k});
    M(:,:,k) = dicomread(dicomList{k});
    sl(k) = dcmHeader{k}.SliceLocation;
end

%% Convert to IEC geometry and scale intensity value
M = permute(single(M),[1 3 2]);
signTransformMat = dcmHeader{1}.ImageOrientationPatient;
signTransformMat = signTransformMat(signTransformMat~=0);

% Inverting LR and AP if needed
if dcmHeader{1}.ImageOrientationPatient(1) ~= 0
    if signTransformMat(1) < 0
        M = M(:,:,end:-1:1);
    end
    if signTransformMat(2) > 0
        M = M(end:-1:1,:,:);
    end
else
    if signTransformMat(1) < 0
        M = M(end:-1:1,:,:);
    end
    if signTransformMat(2) > 0
        M = M(:,:,end:-1:1);
    end
end

% Axial slice ordering
[~,indSlice] = sort(sl);
[~,indSlice] = unique(sl);
M = M(:,indSlice(end:-1:1),:);

% Swap LR and AP if needed
if dcmHeader{1}.ImageOrientationPatient(1) ~= 0
    M = permute(M,[3 2 1]);
end

if ~isnan(waterAtt)
    if min(M(:)) >= 0
        M = M * waterAtt / 1000;
    else
        M = (M + 1000) * waterAtt / 1000;
    end
end

%% Construct header
load(which('mhaHeaderTemplate.mat'));
mhaHeader = mhaHeaderTemplate;
mhaHeader.Dimensions = size(M);
mhaHeader.PixelDimensions = [dcmHeader{1}.PixelSpacing(1),dcmHeader{1}.SliceThickness,dcmHeader{1}.PixelSpacing(2)];
mhaHeader.Offset = -0.5 * (mhaHeader.Dimensions - 1) .* mhaHeader.PixelDimensions;

%% Write to file
MhaWrite(mhaHeader,M,outputFile);

end
