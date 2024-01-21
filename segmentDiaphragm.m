function [dphL,dphR] = segmentDiaphragm(img,header,outputDir)
%% [dphL,dphR] = segmentDiaphragm(img,header,outputDir)
% ------------------------------------------
% FILE   : segmentDiaphragm.m
% AUTHOR : Andy Shieh, ACRF Image X Institute, The University of Sydney
% DATE   : 2018-06-16  Created.
% ------------------------------------------
% PURPOSE
%   Segment diaphragm from an input 3D image. The function assumes the
%   input image is in attenuation values, i.e. minimum valu in the image
%   should be 0.
% ------------------------------------------
% INPUT
%   img:            The input 3D image matrix in IEC geometry (LR, SI, AP)
%   header:         The MHA header struct of the img to be segmented.
%   outputDir:      If this is input, the program will write the left and
%                   right diaphragm masks into "dphL.mha" and "dphR.mha" in
%                   the outputDir folder. If left empty, no file is
%                   written.
%                   Default: empty

% OUTPUT
%   dphL:           Segmentatoin of the left diaphragm (binary mask).
%   dphR:           Segmentatoin of the right diaphragm (binary mask).

%% Intialize some parameters
spacing = header.PixelDimensions;

% The upper and lower bound on the concavity of the parabola as the diaphragm is
% unlikely to have very large or small curvature
PMax = 0.1 * spacing(1)^2 / spacing(2);
PMin = 0.0025 * spacing(1)^2 / spacing(2);

% This is the tolerance in pixel distance for whether to use points
% predicted by the parabola or not. 1 is for right, 2 is for left
tol(1) = ceil(3 * spacing(2) / spacing(1));
tol(2) = ceil(3 * spacing(2) / spacing(1));

% Mininum points to be considered a segment
NMINSEG = 10;

%% Extract body and lung as a reference
fprintf('Segmenting the body and the lungs ......');
tic;

% First, threshold negative values in the image
img(img<0) = 0;

% Find the best tissue cutoff value
cutoffs = linspace(min(img(:)),max(img(:)),100);
for k = 1:length(cutoffs)
    numOfPix(k) = sum(img(:) > cutoffs(k)); % Number of pixels with intensities larger than a sliding cutoff value
end
dN = numOfPix(3:end) - numOfPix(1:end-2); % dN = derivative of numOfPix
dN = [max(abs(dN)),dN,max(abs(dN))]; % Pad the first and last entry
% Cutoff that resulted in less than 10% of pixels being soft tissue is
% very unlikely to be correct
dN(numOfPix/numel(img) < 0.1) = max(abs(dN));
% Then within the likely values, pick the one that is most stable to
% segmentation
indBestCutoff = round(mean(find(abs(dN) == min(abs(dN)))));
softCutoff = cutoffs(indBestCutoff);

% Segment the body and use imfill to fill out holes within and outside the
% body. Also enforce connectivity constraint
body = ~imfill(~imfill(img > softCutoff,'holes'),'holes');
body = findConnected(body);

% Lung is calculated as region that is not part of the body.
% An additional constraint is that it should not be connected to the
% boundary pixels in the LR and AP directions
lung = ~body;
% Do in along the AP direction
for z = 1:size(lung,3)
    slice = lung(:,:,z);
    lconn = bwconncomp(slice);
    for nconn = 1:length(lconn.PixelIdxList)
        clear ind2D;
        [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),lconn.PixelIdxList{nconn});
        if any(ind2D(:,1) == 1) || any(ind2D(:,1) == size(lung,1))
            slice(lconn.PixelIdxList{nconn}) = false;
        end
    end
    lung(:,:,z) = slice;
end
% Also do it along the LR direction
for x = 1:size(lung,1)
    slice = permute(lung(x,:,:),[3,2,1]);
    lconn = bwconncomp(slice);
    for nconn = 1:length(lconn.PixelIdxList)
        clear ind2D;
        [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),lconn.PixelIdxList{nconn});
        if any(ind2D(:,1) == 1) || any(ind2D(:,1) == size(lung,3))
            slice(lconn.PixelIdxList{nconn}) = false;
        end
    end
    lung(x,:,:) = permute(slice,[3,2,1]);
end

% We also want to exlude small cavities. We cannot simply pick the first
% one or two largest connected regions, as the two lungs are sometimes
% connected and sometimes not. We look for the largest connected region
% first, and exclude anything <50% of its volume.
lconn = bwconncomp(lung);
indLargest = find(cellfun(@length,lconn.PixelIdxList) == max(cellfun(@length,lconn.PixelIdxList)));
for nconn = 1:length(lconn.PixelIdxList)
    if length(lconn.PixelIdxList{nconn}) < 0.5 * length(lconn.PixelIdxList{indLargest})
        lung(lconn.PixelIdxList{nconn}) = false;
    end
end

% We are not interested in holes in the lung, so we use imfill to remove
% them
for z = 1:size(lung,3)
    lung(:,:,z) = imfill(lung(:,:,z),'holes');
end

fprintf('COMPLETED using %s\n',seconds2human(toc));
%% Generate estimated diaphragm mask
fprintf('Getting the initial estimated mask ......');
tic;

% First, find the pixels that go from lung to body
estMask = false(size(img));
dph{1} = false(size(img));
dph{2} = false(size(img));
for y = 2:size(estMask,2)
    estMask(:,y,:) = body(:,y,:) & lung(:,y-1,:);
end
% If any pixel is a 8-neighbor of the dphMask pixel, it should also be
% counted as a dphMask pixel. This is to make the algorithm more robust
% against blurring, artifacts, and steep diaphragm gradient
estMaskAug = estMask;
for z = 1:size(estMaskAug,3)
    estMaskAug(2:end-1,2:end-1,z) = estMaskAug(2:end-1,2:end-1,z)   | ...
        estMaskAug(1:end-2,2:end-1,z) | estMaskAug(3:end,2:end-1,z) | ...
        estMaskAug(1:end-2,1:end-2,z) | estMaskAug(3:end,1:end-2,z) | ...
        estMaskAug(1:end-2,3:end,z)   | estMaskAug(3:end,3:end,z)   | ...
        estMaskAug(2:end-1,1:end-2,z) | estMaskAug(2:end-1,3:end,z);
end

fprintf('COMPLETED using %s\n',seconds2human(toc));
%% Finding the optimal starting slice and parabolas
fprintf('Finding best diaphragm candidate in the coronal view ......');
tic;

% We find the coronal slice with the best parabolic segment in the following senses:
% 1. largest LR width. 2. The far end is furthest from the center of the
% image (diaphragm should not be at the center of the image).
% We do this separately for the left and right lung.
% This gives us a starting point of segments that must be the diaphragm
warning('off','all'); % The polynomial fitting will throw warnings that are not of concern
% We do this for the right (1) and left side (2) separately
for nside = 1:2
    optScore = zeros(size(estMask,3),1);
    for z = 1:size(estMask,3)
        slice = estMaskAug(:,:,z);
        sconn = bwconncomp(slice);
        if isempty(sconn.PixelIdxList)
            continue;
        end
        scoreVec = [];
        % Although we use the augmented mask to pick out the connected
        % segment (to avoid disrupted segments due to artifacts), we only
        % want to fit the parabola to points in the original esimated mask
        estMaskIdxList = find(estMask(:,:,z));
        for nconn = 1:length(sconn.PixelIdxList)
            validIdxList = sconn.PixelIdxList{nconn}(ismember(sconn.PixelIdxList{nconn},estMaskIdxList));
            clear ind2D;
            [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),validIdxList);
            % Identify if this is the left and right side
            if (nside == 1 && mean(ind2D(:,1)) > size(estMask,1)/2) || (nside == 2 && mean(ind2D(:,1)) <= size(estMask,1)/2)
                continue;
            end
            % Fit a parabola
            [p,s,mu] = polyfit(ind2D(:,1),ind2D(:,2),2);
            if p(1) / var(ind2D(:,1)) >= PMin && p(1) / var(ind2D(:,1)) <= PMax && size(ind2D,1) >= NMINSEG
                % Find the left most, center, and right most points here
                indLeft = find(ind2D(:,1) == min(ind2D(:,1)));
                indRight = find(ind2D(:,1) == max(ind2D(:,1)));
                xCenter = 0.5 * max(ind2D(:,1)) + 0.5 * min(ind2D(:,1));
                indCenter = find(abs(ind2D(:,1) - xCenter) == min(abs(ind2D(:,1) - xCenter)));
                scoreVec = [scoreVec, ...
                    (max(ind2D(:,1)) - min(ind2D(:,1))) ...         % Span in the LR direction
                    * abs(mean(ind2D(:,1)) - size(img,1)/2) ...    % Segments that are not around the center of the image
                    * sign(max(ind2D(indRight,2)) - min(ind2D(indCenter,2))) * sign(max(ind2D(indLeft,2)) - min(ind2D(indCenter,2)))]; % The center of the diaphragm should be higher than the sides
            end
        end
        if isempty(scoreVec)
            optScore(z) = 0;
        else
            optScore(z) = max(scoreVec);
        end
    end
    z_best = find(optScore == max(optScore)); z_best = z_best(1);
    z_ref(nside) = z_best;
end

% For that coronal slice, we find the parabola that best fit the
% diaphragm curve
for nside = 1:2
    slice = estMaskAug(:,:,z_ref(nside));
    sconn = bwconncomp(slice);
    scoreVec = [];
    nconnVec = [];
    estMaskIdxList = find(estMask(:,:,z_ref(nside)));
    clear p s mu nValid;
    for nconn = 1:length(sconn.PixelIdxList)
        validIdxList = sconn.PixelIdxList{nconn}(ismember(sconn.PixelIdxList{nconn},estMaskIdxList));
        clear ind2D;
        [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),validIdxList);
        if (nside == 1 && mean(ind2D(:,1)) > size(estMask,1)/2) || (nside == 2 && mean(ind2D(:,1)) <= size(estMask,1)/2)
            continue;
        end
        % Fit a parabola
        [p{nconn},s,mu{nconn}] = polyfit(ind2D(:,1),ind2D(:,2),2);
        if p{nconn}(1) / var(ind2D(:,1)) >= PMin && p{nconn}(1) / var(ind2D(:,1)) <= PMax && size(ind2D,1) >= NMINSEG
            % Find the left most, center, and right most points here
            indLeft = find(ind2D(:,1) == min(ind2D(:,1)));
            indRight = find(ind2D(:,1) == max(ind2D(:,1)));
            xCenter = 0.5 * max(ind2D(:,1)) + 0.5 * min(ind2D(:,1));
            indCenter = find(abs(ind2D(:,1) - xCenter) == min(abs(ind2D(:,1) - xCenter)));
            scoreVec = [scoreVec, ...
                (max(ind2D(:,1)) - min(ind2D(:,1))) ...         % Span in the LR direction
                * abs(mean(ind2D(:,1)) - size(img,1)/2) ...    % Segments that are not around the center of the image
                * sign(max(ind2D(indRight,2)) - min(ind2D(indCenter,2))) * sign(max(ind2D(indLeft,2)) - min(ind2D(indCenter,2)))]; % The center of the diaphragm should be higher than the sides
            nconnVec = [nconnVec, nconn];
        end
    end
    [~,sortid] = sort(scoreVec);
    p_side{nside} = p{nconnVec(sortid(end))}; mu_side{nside} = mu{nconnVec(sortid(end))};
    idxList{nside} = sconn.PixelIdxList{nconnVec(sortid(end))};
    idxList{nside} = idxList{nside}(ismember(idxList{nside},estMaskIdxList));
end

% Now for this slice, for every x value, we only keep the point that best matches
% with the parabola of the left or right diaphragm segment. If the matching
% error is larger than the tolerance, we use the average of the point predicted
% by the parabola and the point in the estimated mask. If there is no point
% for that x value, we use the point predicted by the parabola instead.
for nside = 1:2
    clear ind2D;
    [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),idxList{nside});
    slice = false(size(estMask(:,:,1)));
    % We put a 1-pixel margin here because the 8-neighbor augmentation will
    % give us an additional pixel on each side
    for x = min(ind2D(:,1))+1:max(ind2D(:,1))-1
        idx = find(ind2D(:,1) == x);
        err = abs(ind2D(idx,2) - polyval(p_side{nside},ind2D(idx,1),[],mu_side{nside}));
        yofx = polyval(p_side{nside},x,[],mu_side{nside});
        if isempty(err)
            y = round(yofx);
        else
            if min(err) > tol(nside)
                y = round(0.5 * yofx + 0.5 * mean(ind2D(idx(err == min(err)),2)));
            else
                y = round(mean(ind2D(idx(err == min(err)),2)));
            end
        end
        if x >= 1 && x <= size(slice,1) && y >= 1 && y <= size(slice,2)
            slice(x,y) = true;
        end
    end
    
    % Use the original transition mask to filter out incorrection
    % points due to fitting
    pixIdxEst = find(estMaskAug(:,:,z_ref(nside)));
    pixIdxFit = sub2ind(size(slice),ind2D(:,1),ind2D(:,2));
    indRemove = find(~ismember(pixIdxFit,pixIdxEst));
    ind2D(indRemove,:) = [];
    
    % If we end up with less than 2 points, we cannot perform
    % spline interpolation. We probably failed to find the diaphragm at all
    % anyway
    if size(ind2D,1) < 2
        dph{nside}(:,:,z_ref(nside)) = false(size(slice));
        break;
    end
    
    % Finall, we use smoothing spline to smooth out the segment
    clear ind2D;
    [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),find(slice));
    smoothfit = fit(ind2D(:,1),ind2D(:,2),'smoothingspline','SmoothingParam',0.005);
    yy = feval(smoothfit,ind2D(:,1));
    slice = false(size(slice));
    for k = 1:size(ind2D,1)
        xidx = ind2D(k,1);
        yidx = round(yy(k));
        if yidx >= 1 && yidx <= size(slice,2)
            slice(xidx,yidx) = true;
        end
    end
    dph{nside}(:,:,z_ref(nside)) = slice;
end

fprintf('COMPLETED using %s\n',seconds2human(toc));
%% Working towards the first and the last slices
fprintf('Propagating the results to other coronal slices ......');
tic;

% Now we work from this slice towards the first slice and exclude all
% segments have both an mean fitting error larger than the tolerance with
% the parabolas in its consecutive frame, and <25% overlap with the diaphragm in
% the previous slice within +-1 SI pixel. The parabolas are updated as we
% go through the slices
% We also work towards the end slice as well. To avoid repeating codes, we
% put this in a loop where ndir = 1 is working towards the first slice and
% ndir = 2 is working towards the last slice.
for ndir = 1:2
    for nside = 1:2
        pz = p_side{nside}; muz = mu_side{nside};
        if ndir == 1
            zvec = (z_ref(nside)-1):-1:1;
        else
            zvec = (z_ref(nside)+1):size(img,3);
        end
        for z = zvec
            % Indices for the diaphragm in the reference slice
            if ndir == 1
                indRef = find(dph{nside}(:,:,z+1));
            else
                indRef = find(dph{nside}(:,:,z-1));
            end
            
            % Get segment candidates for ths slice
            slice = estMask(:,:,z);
            sconn = bwconncomp(slice);
            if length(sconn.PixelIdxList) == 0
                break;
            end
            for nconn = 1:length(sconn.PixelIdxList)
                % Fit with reference parabola
                clear ind2D;
                [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),sconn.PixelIdxList{nconn});
                if (nside == 1 && mean(ind2D(:,1)) > size(estMask,1)/2) || (nside == 2 && mean(ind2D(:,1)) <= size(estMask,1)/2)
                    slice(sconn.PixelIdxList{nconn}) = false;
                    continue;
                end
                % Fitting error
                err = norm(ind2D(:,2) - polyval(pz,ind2D(:,1),[],muz)) / sqrt(length(sconn.PixelIdxList{nconn}));
                
                % Overlap
                overlap = sum(sum( ...
                    ismember(sconn.PixelIdxList{nconn}, indRef) | ...
                    ismember(sconn.PixelIdxList{nconn} - size(slice,1),indRef) | ...
                    ismember(sconn.PixelIdxList{nconn} + size(slice,1),indRef) )) / ...
                    length(sconn.PixelIdxList{nconn});
                
                if (err > tol(nside) && overlap < 0.25) || isnan(err)
                    slice(sconn.PixelIdxList{nconn}) = false;
                end
            end
            
            % If this slice has less than NMINSEG, we should stop the loop here
            if sum(sum(slice)) < NMINSEG
                if ndir == 1
                    dph{nside}(:,:,1:z) = false;
                else
                    dph{nside}(:,:,z:end) = false;
                end
                break;
            end
            
            % If the remaining points do not form a parabola with positive
            % concavity, the loop should stop at this slice
            clear ind2D;
            [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),find(slice));
            [p,s,mu] = polyfit(ind2D(:,1),ind2D(:,2),2);
            if p(1) <= 0
                if ndir == 1
                    dph{nside}(:,:,1:z) = false;
                else
                    dph{nside}(:,:,z:end) = false;
                end
                break;
            end
            
            % Keep only one point that matches best with the parabola for each
            % x value. Here we also make it such that the first column in ind2D
            % is in ascending order
            ind2DTmp = [];
            yPara = polyval(p,min(ind2D(:,1)):max(ind2D(:,1)),[],mu);
            for x = min(ind2D(:,1)):max(ind2D(:,1))
                idx = find(ind2D(:,1) == x);
                if isempty(idx)
                    continue;
                end
                if length(idx) == 1
                    ind2DTmp = [ind2DTmp;x,ind2D(idx,2)];
                    continue;
                end
                fitError = abs(yPara(x - min(ind2D(:,1)) + 1) - ind2D(idx,2));
                indBest = find(fitError == min(fitError));
                if length(indBest) > 1
                    indBest = indBest(ind2D(idx(indBest),2) == max(ind2D(idx(indBest),2)));
                end
                ind2DTmp = [ind2DTmp;x,ind2D(idx(indBest),2)];
            end
            ind2D = ind2DTmp;
            % And then update the parabola
            [p,~,mu] = polyfit(ind2D(:,1),ind2D(:,2),2);
            
            % Sometimes the segments curve up on the inner side of the lung
            % usually because we accidentally included the lateral lung
            % boundary. We remove those points here
            xCenter = - mu(2) * p(2) / p(1) / 2 + mu(1);
            if nside == 1 % Right
                idxSide = find(ind2D(:,1) > xCenter + mu(2) * sqrt(tol(nside) / p(1)));
                dy = diff(ind2D(:,2));
                idxdyPos = find(dy < 0);
                indCurveUp = idxdyPos(ismember(idxdyPos,idxSide));
                if isempty(indCurveUp)
                    indRemove = [];
                else
                    indRemove = indCurveUp(1):size(ind2D,1);
                end
            else % Left
                idxSide = find(ind2D(:,1) < xCenter - mu(2) * sqrt(tol(nside) / p(1)));
                dy = ind2D(1:end-1,2) - ind2D(2:end,2);
                idxdyPos = find(dy < 0);
                indCurveUp = idxdyPos(ismember(idxdyPos,idxSide));
                if isempty(indCurveUp)
                    indRemove = [];
                else
                    indRemove = 1:indCurveUp(end);
                end
            end
            ind2D(indRemove,:) = [];
            
            % If we end up with less than 2 points, we cannot perform
            % spline interpolation. We should probably stop the loop here
            % anyway
            if size(ind2D,1) < 2
                if ndir == 1
                    dph{nside}(:,:,1:z) = false;
                else
                    dph{nside}(:,:,z:end) = false;
                end
                break;
            end
            
            % We use spline interpolation to smooth out the points
            smoothfit = fit(ind2D(:,1),ind2D(:,2),'smoothingspline','SmoothingParam',0.005);
            yy = feval(smoothfit,ind2D(1,1):ind2D(end,1));
            slice = false(size(slice));
            % We only interpolate for points not originally segmented if
            % they account for <=25% of the LR range
            missingTooMuch = sum(ismember(ind2D(:,1),ind2D(1,1):ind2D(end,1))) ...
                < 0.75 * (ind2D(end,1) - ind2D(1,1) + 1);
            for xidx = ind2D(1,1):ind2D(end,1)
                if missingTooMuch && ~ismember(xidx,ind2D(:,1))
                    continue;
                end
                yidx = round(yy(xidx - ind2D(1,1) + 1));
                if yidx >= 1 && yidx <= size(slice,2)
                    slice(xidx,yidx) = true;
                end
            end
            
            % Use the original transition mask to filter out incorrection
            % points due to fitting
            pixIdxEst = find(estMaskAug(:,:,z));
            pixIdxFit = find(slice);
            slice(pixIdxFit(~ismember(pixIdxFit,pixIdxEst))) = false;
            
            dph{nside}(:,:,z) = slice;
            
            % Update parabola
            clear ind2D;
            [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),find(slice));
            if nside == 1
                [pz,~,muz] = polyfit(ind2D(find(ind2D(:,1) <= size(estMask,1)/2),1),ind2D(find(ind2D(:,1) <= size(estMask,1)/2),2),2);
            else
                [pz,~,muz] = polyfit(ind2D(find(ind2D(:,1) > size(estMask,1)/2),1),ind2D(find(ind2D(:,1) > size(estMask,1)/2),2),2);
            end
        end
    end
end

fprintf('COMPLETED using %s\n',seconds2human(toc));
%% Finally, we smooth the segmentation in AP view as well
fprintf('Apply spline smoothing in AP view ......');
tic;

for nside = 1:2
    for x = 1:size(img,1)
        slice = permute(dph{nside}(x,:,:),[3 2 1]);
        
        if sum(sum(slice)) < NMINSEG
            slice(:) = false;
        else
            clear ind2D;
            [ind2D(:,1),ind2D(:,2)] = ind2sub(size(slice),find(slice));
            % Spline interpolation to smooth out the points
            smoothfit = fit(ind2D(:,1),ind2D(:,2),'smoothingspline','SmoothingParam',0.005);
            xstart = min(ind2D(:,1));
            xend = max(ind2D(:,1));
            yy = feval(smoothfit,xstart:xend);
            slice = false(size(slice));
            % We only interpolate for points not originally segmented if
            % they account for <=25% of the LR range
            missingTooMuch = sum(ismember(ind2D(:,1),xstart:xend)) ...
                < 0.75 * (xend - xstart + 1);
            for xidx = xstart:xend
                if missingTooMuch && ~ismember(xidx,ind2D(:,1))
                    continue;
                end
                yidx = round(yy(xidx - xstart + 1));
                if yidx >= 1 && yidx <= size(slice,2)
                    slice(xidx,yidx) = true;
                end
            end
            
            % Use the original transition mask to filter out incorrection
            % points due to fitting
            pixIdxEst = find(permute(estMaskAug(x,:,:),[3 2 1]));
            pixIdxFit = find(slice);
            slice(pixIdxFit(~ismember(pixIdxFit,pixIdxEst))) = false;
        end
        dph{nside}(x,:,:) = permute(slice,[3 2 1]);
    end
    
    % Usually diaphragm around the spine is hard to segment due to
    % contrast. We eliminate these points (within tol pixels of the most
    % posterior side)
    pixIdxAll = find(dph{nside});
    clear ind3D;
    [ind3D(:,1),ind3D(:,2),ind3D(:,3)] = ind2sub(size(img),find(dph{nside}));
    dph{nside}(pixIdxAll(ind3D(:,3) < min(ind3D(:,3)) + tol(nside))) = false;
end

fprintf('COMPLETED using %s\n',seconds2human(toc));
%% And then assign the results to output
dphR = dph{1};
dphL = dph{2};

%% Write to files if outputDir is provided
if nargin < 3 || isempty(outputDir)
    return;
end

fprintf('Processing 3D diaphragm model files ......');
tic;

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

header.DataType = 'uchar';
header.BitDepth = 8;
MhaWrite(header,dphL,fullfile(outputDir,'dphL.mha'));
MhaWrite(header,dphR,fullfile(outputDir,'dphR.mha'));

fprintf('COMPLETED using %s\n',seconds2human(toc));