function [trj3D,metricVal,trackFrame] = ...
    trackDiaphragm_corrcoef(projList,model,geometryFile,flipTag,invertTag,r,excMargin,frameInd)
%% Update to trackDiaphragm.m to score candidate positions based on Pearson correlation
% for model: suggest get2DDiaphragmModel.m
% for pcVec: suggest getDph3DShift.m
%% Input check
if nargin < 4
    flipTag = false;
end

if nargin < 5
    invertTag = false;
end

if nargin < 6
    r = 3;
end

if nargin < 7
    excMargin = [10,10,10,20];
end

if nargin < 8
    frameInd = 1:length(projList);
end

%% Hard coded parameters
initr = 10; % Initial search radius
w = 2;      % Margin (in 3D space mm) for calculating neighboring intensity

%% Initialization
fprintf('Initializing......');
tic;

ang = ReadAnglesFromRTKGeometry(geometryFile);
sid = 785;%ReadTagFromRTKGeometry('Geometry.xml','SourceToIsocenterDistance');
sdd = 1300;%ReadTagFromRTKGeometry('Geometry.xml','SourceToDetectorDistance');

% Pre-calculated v-sum-weighted map
% The V-sum-weighted map makes sure each u value is treated equally
map = model.Map2D_R + model.Map2D_L;
mapWeighted = map;
for k = 1:size(mapWeighted,3)
    mapWeighted(:,:,k) = mapWeighted(:,:,k) ./ (sum(mapWeighted(:,:,k),2) * ones(1,size(mapWeighted,2)));
end
mapWeighted(isnan(mapWeighted) | isinf(mapWeighted)) = 0;

% Pre-allocation
[~,P] = ProjRead(projList{frameInd(1)});
trackedMap = single(zeros([size(P),length(frameInd)]));

dv_prev = 0; % dv is SI shift in projcetion pixel coordinate
fh = figure('units','normalized','outerposition',[0.1,0.1,0.8,0.8]);
nCount = 0;

fprintf('%f seconds\n',toc);
%% Start tracking
for n = frameInd
    fprintf('Frame#%05d......',n);
    nCount = nCount + 1;
    % Read projection
    tic;
    [~,~,ext] = fileparts(projList{n});
    [projInfo,P] = ProjRead(projList{n});
    P = single(P);

    if flipTag
        P = P(:,end:-1:1);
    end
    if invertTag
        P = 65535 - P;
    end
    if strcmpi(ext,'.hnc')
        P = log(65536 ./ (P + 1));
        spacingX = projInfo.dIDUResolutionX;
        spacingY = projInfo.dIDUResolutionY;
    elseif strcmpi(ext,'.hnd')
        P = log(65536 ./ (P + 1));
        % Hnd files acquired from iTools often have the resolution fiels
        % written incorrectly. Enforcing 0.388 mm pixel size here
        spacingX = 0.388;
        spacingY = 0.388;
    elseif strcmpi(ext,'.his')
        P = log(65536 ./ (65536 - P));
        spacingX = projInfo.PixelSpacingX;
        spacingY = projInfo.PixelSpacingY;
    end
    fprintf('%f seconds...... ',toc);
    
    fprintf('Tracking......');
    tic;
    
    % Search radius in terms of projection pixel
    if n == frameInd(1)
        rPix = round(initr * sdd / sid / spacingY);
    else
        rPix = round(r * sdd / sid / spacingY);
    end
    
    % Margin for calculating neighboring intensity
    wPix = round(w * sdd / sid / spacingY);
    
    % Exclude margins
    P(1:excMargin(1),:) = nan;
    P(end-excMargin(3)+1:end,:) = nan;
    P(:,1:excMargin(2)) = nan;
    P(:,end-excMargin(4)+1:end) = nan;
    
   % Pre-calculate w-Pixel averaged difference map
    imgAvgAbv = zeros(size(P));
    imgAvgBlw = zeros(size(P));
    for kw = 1:wPix
        imgAvgAbv(excMargin(1)+1:end-excMargin(3),wPix+1+excMargin(2):end-excMargin(4)) = imgAvgAbv(excMargin(1)+1:end-excMargin(3),wPix+1+excMargin(2):end-excMargin(4)) + P(excMargin(1)+1:end-excMargin(3),kw+excMargin(2):end-wPix+kw-excMargin(4)-1) / wPix;
        imgAvgBlw(excMargin(1)+1:end-excMargin(3),1+excMargin(2):end-wPix+1-excMargin(4)) = imgAvgBlw(excMargin(1)+1:end-excMargin(3),1+excMargin(2):end-wPix+1-excMargin(4)) + P(excMargin(1)+1:end-excMargin(3),kw+excMargin(2):end-wPix+kw-excMargin(4)) / wPix;
        % Boundary pixels
        if kw < wPix
            for kb = 1:kw
                imgAvgAbv(excMargin(1)+1:end-excMargin(3),kw+1+excMargin(2)) = imgAvgAbv(excMargin(1)+1:end-excMargin(3),kw+1+excMargin(2)) + P(excMargin(1)+1:end-excMargin(3),kw+1-kb+excMargin(2)) / kw;
                imgAvgBlw(excMargin(1)+1:end-excMargin(3),end-kw+1-excMargin(4)) = imgAvgBlw(excMargin(1)+1:end-excMargin(3),end-kw+1-excMargin(4)) + P(excMargin(1)+1:end-excMargin(3),end-kw+kb-excMargin(4)) / kw;
            end
        end
    end
    diffImg = imgAvgBlw - imgAvgAbv;
    diffImg(:,1+excMargin(2)) = 0;
    
    % Select which frame of the 2D model to use
    %ang = projInfo.dCBCTPositiveAngle;
    indModel = find(abs(mod(ang(n),360) - mod(model.Angles,360)) == min(abs(mod(ang(n),360) - mod(model.Angles,360))));
    %indModel = find(abs(mod(ang,360) - mod(model.Angles,360)) == min(abs(mod(ang,360) - mod(model.Angles,360))));
    indModel = indModel(1);
    map_this = map(:,:,indModel);
    mapWeighted_this = mapWeighted(:,:,indModel);
    
    % Go through search window
    dvCand = (dv_prev-rPix):(dv_prev+rPix);
    % Calculate the lateral shift for this particular SI shift by
    % using the principle component vector
     duCand = round((cosd(ang(n)) * model.PCVec(1) - sind(ang(n)) * model.PCVec(3)) ...
         * dvCand / model.PCVec(2) * spacingY / spacingX);
%duCand = round((cosd(ang) * model.PCVec(1) - sind(ang) * model.PCVec(3)) ...
        %* dvCand / model.PCVec(2) * spacingY / spacingX);
    metricVec = zeros(length(dvCand),1);
    % Find 2D index of non-zero map pixels
    mapIdx = find(mapWeighted_this>0);
    [mapX,mapY] = ind2sub(size(mapWeighted_this),mapIdx);
    for k = 1:length(dvCand)
        % Calculate shifted pixel index
        shiftedX = mapX + duCand(k);
        shiftedY = mapY + dvCand(k);
        % Find pixels within the valid FOV
        indIncl = find(shiftedX >= 1 + excMargin(1) & shiftedX <= size(P,1) - excMargin(3) & shiftedY >= 1 + excMargin(2) & shiftedY <= size(P,2) - excMargin(4));
        % Convert x y to 1D index
        shiftedIdx = sub2ind(size(mapWeighted_this),shiftedX(indIncl),shiftedY(indIncl));
        
        % Calculate tracking metric
        metricVec(k) = sum(diffImg(shiftedIdx).* mapWeighted_this(mapIdx(indIncl)));
    end
    
    % Find best match
    indBest = find(metricVec == max(metricVec));
    indBest = round(mean(indBest));
    metricVal(n) = metricVec(indBest);
    dv_prev = dvCand(indBest);
    du_prev = duCand(indBest);
    
    % Convert to patient coordinate
    trj3D(n,2) = dv_prev * spacingY * sid / sdd;
    trj3D(n,[1,3]) = trj3D(n,2) * [model.PCVec(1),model.PCVec(3)] / model.PCVec(2);

    fprintf('%f seconds...... ',toc);
    
    % Create results for visualisation
    fprintf('Visualizing......');
    tic;
    % We record and use the non-weighted map for visualization
    mapShifted = imtranslate(map_this,[dv_prev,du_prev]);
    trackedMap(:,:,nCount) = mapShifted;
    
    winMin = prctile(P(:),1); winMax = prctile(P(:),99);
    pVis = (P - winMin) / (winMax - winMin) * 255;
    pVis(isnan(pVis)) = 0;
    pVis(pVis>255) = 255; pVis(pVis<0) = 0;
    pVis = uint8(round(pVis'));
    maxMapVal = max(mapShifted(:));
    pB = pVis;
    pVis = pVis - uint8(150 * mapShifted' / maxMapVal);
    pVis(:,:,2) = pVis; pVis(:,:,3) = pB;
    imagesc(pVis); axis image; axis off;
    [~,fileName] = fileparts(projList{n});
    text(20,30,['Frame: ', num2str(n),'   SI shift: ',num2str(trj3D(n,2),'%.1f'),' mm'],'color','blue','FontSize',14);
    trackFrame(nCount) = getframe(fh);
    
    fprintf('%f seconds......\n',toc);
end

close(fh);

%% Write tracking movie to .avi
vw = VideoWriter(fullfile(pwd,'trackDph.avi'));
vw.FrameRate = 10;
vw.open;
vw.writeVideo(trackFrame);
vw.close;

end
