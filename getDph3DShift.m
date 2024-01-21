function [shift3D,shiftedMask] = getDph3DShift(dphMask,img,header,rLR,rSI,rAP)
%% [shift3D,shiftedMask] = getDph3DShift(dphMask,img,header,rLR,rSI,rAP)
% ------------------------------------------
% FILE   : getDph3DShift.m
% AUTHOR : Andy Shieh, ACRF Image X Institute, The University of Sydney
% DATE   : 2018-07-12 Created.
% ------------------------------------------
% PURPOSE
%  Calculate 3D shift of the diaphragm model in order to match with the
%  input image.
% ------------------------------------------
% INPUT
%   dphMask:        The binary mask the 3D diaphragm model. Suggest to use
%                   the joint mask of the left and right diaphragm segments.
%   img:            The 3D image to be matched to. Must have the same
%                   dimension and occupy the same physical space as the
%                   dphMask.
%   header:         MHA header struct of either dphMask or img.
%   rLR:            (Optional) Search vector for LR shift (mm). Default: 0
%   rSI:            (Optional) Search vector for SI shift (mm). Default: -20:20 mm
%   rAP:            (Optional) Search vector for AP shift (mm). Default: -20:20 mm
% OUTPUT
%   shift3D:        The calculated 3D shift (1x3 vector).

%% Input check and hard-coded paramters
if ~isequal(size(dphMask),size(img))
    error('ERROR: dphMask and img must have the same dimensions.');
end

if nargin < 4
    rLR = 0;
end

if nargin < 5
    rSI = -20:20;
end

if nargin < 6
    rAP = -20:20;
end

% Margins used to calculate mean pixel intensity above and below the
% diaphragm (5 mm)
w = round(5 / header.PixelDimensions(2));

%% Fitting the diaphragm
% Convert search vector to pixel index
rpixLR = unique(round(rLR / header.PixelDimensions(1)));
rpixSI = unique(round(rSI / header.PixelDimensions(2)));
rpixAP = unique(round(rAP / header.PixelDimensions(3)));

% Pixel index list of the diaphragm mask
[dphIdx(:,1),dphIdx(:,2),dphIdx(:,3)] = ind2sub(size(dphMask),find(dphMask));

% Precalculate the w-pixel mean image
imgAvgAbv = zeros(size(img));
imgAvgBlw = zeros(size(img));
for kw = 1:w
    imgAvgAbv(:,w:end,:) = imgAvgAbv(:,w:end,:) + img(:,kw:end-w+kw,:) / w;
    imgAvgBlw(:,1:end-w+1,:) = imgAvgBlw(:,1:end-w+1,:) + img(:,kw:end-w+kw,:) / w;
    % Boundary pixels
    if kw < w
        for kb = 1:kw
            imgAvgAbv(:,kw,:) = imgAvgAbv(:,kw,:) + img(:,kw-kb+1,:) / kw;
            imgAvgBlw(:,end-kw+1,:) = imgAvgAbv(:,end-kw+1,:) + img(:,end-kw+kb,:) / kw;
        end
    end
end

% Loop through search vectors
metricVals = zeros([length(rpixLR),length(rpixSI),length(rpixAP)]);
for kLR = 1:length(rpixLR)
    for kSI = 1:length(rpixSI)
        for kAP = 1:length(rpixAP)
            x = dphIdx(:,1) + rpixLR(kLR);
            y = dphIdx(:,2) + rpixSI(kSI);
            z = dphIdx(:,3) + rpixAP(kAP);
            % y > 1 instead of >=1 because we need at least a pixel above
            % to calculate difference
            idxValid = x >= 1 & x <= size(img,1) & y > 1 & y <= size(img,2) & z >= 1 & z <= size(img,3);
            x = x(idxValid);
            y = y(idxValid);
            z = z(idxValid);
            idxBlw = sub2ind(size(img),x,y,z);
            idxAbv = sub2ind(size(img),x,y-1,z);
            metricVals(kLR,kSI,kAP) = mean(imgAvgBlw(idxBlw) - imgAvgAbv(idxAbv));
        end
    end
end
[idxBest(:,1),idxBest(:,2),idxBest(:,3)] = ind2sub(size(metricVals),find(metricVals == max(metricVals(:))));

shift3D = [mean(rpixLR(idxBest(:,1))),mean(rpixSI(idxBest(:,2))),mean(rpixAP(idxBest(:,3)))] .* header.PixelDimensions;

% Generate the new shifted mask
shiftedMask = imtranslate(dphMask,[round(mean(rpixSI(idxBest(:,2)))),round(mean(rpixLR(idxBest(:,1)))),round(mean(rpixAP(idxBest(:,3))))]);