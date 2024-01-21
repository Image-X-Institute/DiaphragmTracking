function model = get2DDiaphragmModel(dphLFile,dphRFile,pcVec,volOffset,sid,sdd,detOffset,projSize,projSpacing,outputDir)
%% model = get2DDiaphragmModel(dphLFile,dphRFile,pcVec,volOffset,sid,sdd,detOffset,projSize,projSpacing,outputDir)
% ------------------------------------------
% FILE   : get2DDiaphragmModel.m
% AUTHOR : Andy Shieh, ACRF Image X Institute, The University of Sydney
% DATE   : 2018-07-09  Created.
% ------------------------------------------
% PURPOSE
%  Generate 2D diaphragm model in the 2D image space.
%  The model is generated for 0-360 degree gantry angle with 0.5 degree
%  increment.
% ------------------------------------------
% INPUT
%   dphLFile:       MHA file of the 3D model of the left diaphragm.
%   dphRFile:       MHA file of the 3D model of the right diaphragm.
%   pcVec:          The principle component vector (3x1) of 3D diaphragm
%                   movement. This can be done by using getDph3DShift to
%                   match the exhale diaphragm mask with the inhale CT. The
%                   magnitude of the vector is irrelevant.
%                   A good vector to try is [0,0.7880,0.6156], but the
%                   optimal strategy is to get this using getDph3DShift for
%                   each individual patient.
%   volOffset:      1x3 vector of the offset needed to shift the diaphragm
%                   model to align with patient position.
%   sid:            Source-to-isocenter distance (mm).
%   sdd:            Source-to-detector distance (mm).
%   detOffset:      1x2 vector of the [lateral,vertical] detector offset
%                   (mm).
%   projSize:       1x2 vector of the projection image dimension (int).
%   projSpacing:    1x2 vector of the projection pixel sizes (double).
%   outputDir:      (Optional) If input, the program will save the model as
%                   Dph2DModel.mat in outputDir.
%                   Default: does not save

% OUTPUT
%   model:          The model in a struct. The struct records the geometric
%                   specificis, geometry matrices, and the 2D diaphragm
%                   points (in pixel index) for each angular view.

%% Initialize parameters and output sturct
model.SID = sid;
model.SDD = sdd;
model.VolOffset = volOffset;
model.DetOffset = detOffset;
model.ProjSize = projSize;
model.ProjSpacing = projSpacing;
model.Angles = 0:0.5:(360-0.5);
model.PCVec = pcVec;

%% Generate the RTK geometry matrices
fprintf('Generating RTK geometry matrices ......');
tic;

rtksimulatedgeometry = which('rtksimulatedgeometry.exe');
if isempty(rtksimulatedgeometry)
    rtksimulatedgeometry = 'rtksimulatedgeometry';
end
geoFile = [tempname,'.xml'];
system(['"',rtksimulatedgeometry,'" ','-f 0 -n 720 -a 360 ',...
    '--sid ',num2str(sid,'%f'),' ',...
    '--sdd ',num2str(sdd,'%f'),' ',...
    '--proj_iso_x ',num2str(detOffset(1),'%f'),' ',...
    '--proj_iso_y ',num2str(detOffset(2),'%f'),' ',...
    '-o "',geoFile,'"']);
model.G = ReadRTKGeometryMatrices(geoFile);

if nargin < 8
    system(['del "',geoFile,'"']);
end

fprintf('COMPLETED using %s\n',seconds2human(toc));

%% Reading the 3D diaphragm model
fprintf('Reading the 3D diaphragm model ......');
tic;

% From here on we assign 1 to the right diaphragm and 2 to the left.
% This allows us to process the two diaphragms in a loop without repeating
% codes

[header{1},dph{1}] = MhaRead(dphRFile);
[header{2},dph{2}] = MhaRead(dphLFile);

% Convert binary mask into 3D points (in physical coordinate mm)
for nside = 1:2
    ptsIdx = find(dph{nside});
    [pts3D{nside}(:,1),pts3D{nside}(:,2),pts3D{nside}(:,3)] = ...
        ind2sub(size(dph{nside}),ptsIdx);
    pts3D{nside} = (pts3D{nside} - 1) .* (ones(length(ptsIdx),1) * header{nside}.PixelDimensions) + ...
        ones(length(ptsIdx),1) * (header{nside}.Offset + volOffset);
end

fprintf('COMPLETED using %s\n',seconds2human(toc));

%% Projecting 3D points to 2D projection space
fprintf('Projecting 3D points to 2D projection space ......');
tic;

for nside = 1:2
    % Projection space pixel margin to account for finite 3D voxel size
    wx = 0.5 * max(header{nside}.PixelDimensions(1),header{nside}.PixelDimensions(3)) * sdd / sid / projSpacing(1);
    wy = 0.5 * header{nside}.PixelDimensions(2) * sdd / sid / projSpacing(2);
    
    map2D = single(zeros([projSize(1),projSize(2),size(model.G,3)]));
    for k = 1:size(model.G,3)
        pixIdx2D{k} = [];
        slice = zeros(projSize);
        pts2D = rtk3DTo2D(model.G(:,:,k),pts3D{nside});
        pts2D = pts2D ./ (ones(size(pts2D,1),1) * projSpacing) + ones(size(pts2D,1),1) * (projSize + 1) / 2;
        
        % Sum points up in the projection space
        for np = 1:size(pts2D)
            startU = max(1,floor(pts2D(np,1) - wx));
            endU = min(projSize(1),ceil(pts2D(np,1) + wx));
            startV = max(1,floor(pts2D(np,2) - wy));
            endV = min(projSize(2),ceil(pts2D(np,2) + wy));
            slice(startU:endU,startV:endV) = slice(startU:endU,startV:endV) + 1;
        end
        
        map2D(:,:,k) = slice;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Older implementation of only keeping a diaphragm curve. Not
        %%% used now.
        %         % For each u value, we only keep the point with the maximum value
        %         for u = 1:projSize(1)
        %             if sum(slice(u,:)) == 0
        %                 continue;
        %             end
        %             indV = find(slice(u,:) == max(slice(u,:)));
        %             pixIdx2D{k} = [pixIdx2D{k};u,indV(1)];
        %         end
        %         % Spline smoothing
        %         if size(pixIdx2D{k},1) > 1
        %             smoothfit = fit(pixIdx2D{k}(:,1),pixIdx2D{k}(:,2),'smoothingspline','SmoothingParam',0.005);
        %             pixIdx2D{k}(:,2) = round(feval(smoothfit,pixIdx2D{k}(:,1)));
        %             % Remove invalid points
        %             pixIdx2D{k}(pixIdx2D{k}(:,2) < 1 | pixIdx2D{k}(:,2) > projSize(2),:) = [];
        %         end
        %%% Older implementation of only keeping a diaphragm curve. Not
        %%% used now.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Older implementation of only keeping a diaphragm curve. Not
    %%% used now.
    % Spline smoothing in the k direction as well
    %     for x = 1:projSize(1)
    %         kprof = [];
    %         for k = 1:size(model.G,3)
    %             if isempty(pixIdx2D{k})
    %                 continue;
    %             end
    %             idx = find(pixIdx2D{k}(:,1) == x);
    %             if ~isempty(idx)
    %                 y = pixIdx2D{k}(idx,2);
    %                 kprof = [kprof; k, y, idx];
    %             end
    %         end
    %         if size(kprof,1) > 1
    %             smoothfit = fit(kprof(:,1),kprof(:,2),'smoothingspline','SmoothingParam',0.005);
    %             kprof(:,2) = round(feval(smoothfit,kprof(:,1)));
    %             for nk = 1:size(kprof,1)
    %                 if kprof(nk,2) < 1 || kprof(nk,2) > projSize(2)
    %                     pixIdx2D{kprof(nk,1)}(kprof(nk,3),:) = [];
    %                 else
    %                     pixIdx2D{kprof(nk,1)}(kprof(nk,3),2) = kprof(nk,2);
    %                 end
    %             end
    %         end
    %     end
    %     Put results to mask
    %     for k = 1:size(model.G,3)
    %         maskTmp = false(projSize);
    %         maskTmp(sub2ind(projSize,pixIdx2D{k}(:,1),pixIdx2D{k}(:,2))) = true;
    %         mask2D(:,:,k) = maskTmp;
    %     end
    %%% Older implementation of only keeping a diaphragm curve. Not
    %%% used now.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nside == 1
        model.Map2D_R = map2D;
    else
        model.Map2D_L = map2D;
    end
end

fprintf('COMPLETED using %s\n',seconds2human(toc));

%% Saving
if nargin >= 10
    fprintf('Saving the 2D model ......');
    tic;
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir,'Dph2DModel.mat'),'model','-v7.3');
    
    copyfile(geoFile,fullfile(outputDir,'ModelGeometry.xml'));
    system(['del "',geoFile,'"']);
    
    fprintf('COMPLETED using %s\n',seconds2human(toc));
end