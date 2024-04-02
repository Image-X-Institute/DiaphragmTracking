function [info M] = AttRead(filename)
%% [info M] = AttRead(filename)
% ------------------------------------------
% FILE   : AttRead.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2013-07-05  Created.
% ------------------------------------------
% PURPOSE
%   Read the header and image body of a customary image format which
%   contains a typical Varian Hnc header and attenuation value image in
%   single format
% ------------------------------------------
% INPUT
%   filename : The full path to the hnc file.
% ------------------------------------------
% OUTPUT
%   info  : Header information in a struct.
%   M     : The image stored in a 2D matrix (single).
% ------------------------------------------

%% Checking input arguments & Opening the file

M = [];

% If no input filename => Open file-open-dialog
if nargin < 1
    
    % go into a default directory
    DefaultDir = pwd;
    
    % get the input file (hnc) & extract path & base name
    [FileName,PathName] = uigetfile( {'*.att','Attenuation Image Files (*.att)'}, ...
        'Select an image file', ...
        DefaultDir);
    
    % catch error if no file selected
    if isnumeric(FileName)
        error('ERROR: No file selected. \n');
        info = struct([]);
        return
    end
    
    % make same format as input
    filename = fullfile(PathName, FileName); 
end

filename = strtrim(filename);

% Open the file
fid = fopen(filename,'r');

    % catch error if failure to open the file
    if fid == -1
        fprintf('ERROR : Failure in opening the file. \n');
        info = struct([]);
        return
    end

%% Reading the header (same for hnc and hnd)

info.('bFileType') = fread(fid,32,'uint8=>char')';
info.('uiFileLength') = fread(fid,1,'uint32')';
info.('bChecksumSpec') = fread(fid,4,'uint8=>char')';
info.('uiCheckSum') = fread(fid,1,'uint32')';
info.('bCreationDate') = fread(fid,8,'uint8=>char')';
info.('bCreationTime') = fread(fid,8,'uint8=>char')';
info.('bPatientID') = fread(fid,16,'uint8=>char')';
info.('uiPatientSer') = fread(fid,1,'uint32')';
info.('bSeriesID') = fread(fid,16,'uint8=>char')';
info.('uiSeriesSer') = fread(fid,1,'uint32')';
info.('bSliceID') = fread(fid,16,'uint8=>char')';
info.('uiSliceSer') = fread(fid,1,'uint32')';
info.('uiSizeX') = fread(fid,1,'uint32')';
info.('uiSizeY') = fread(fid,1,'uint32')';
info.('dSliceZPos') = fread(fid,1,'double')';
info.('bModality') = fread(fid,16,'uint8=>char')';
info.('uiWindow') = fread(fid,1,'uint32')';
info.('uiLevel') = fread(fid,1,'uint32')';
info.('uiPixelOffset') = fread(fid,1,'uint32')';
info.('bImageType') = fread(fid,4,'uint8=>char')';
info.('dGantryRtn') = fread(fid,1,'double')';
info.('dSAD') = fread(fid,1,'double')';
info.('dSFD') = fread(fid,1,'double')';
info.('dCollX1') = fread(fid,1,'double')';
info.('dCollX2') = fread(fid,1,'double')';
info.('dCollY1') = fread(fid,1,'double')';
info.('dCollY2') = fread(fid,1,'double')';
info.('dCollRtn') = fread(fid,1,'double')';
info.('dFieldX') = fread(fid,1,'double')';
info.('dFieldY') = fread(fid,1,'double')';
info.('dBladeX1') = fread(fid,1,'double')';
info.('dBladeX2') = fread(fid,1,'double')';
info.('dBladeY1') = fread(fid,1,'double')';
info.('dBladeY2') = fread(fid,1,'double')';
info.('dIDUPosLng') = fread(fid,1,'double')';
info.('dIDUPosLat') = fread(fid,1,'double')';
info.('dIDUPosVrt') = fread(fid,1,'double')';
info.('dIDUPosRtn') = fread(fid,1,'double')';
info.('dPatientSupportAngle') = fread(fid,1,'double')';
info.('dTableTopEccentricAngle') = fread(fid,1,'double')';
info.('dCouchVrt') = fread(fid,1,'double')';
info.('dCouchLng') = fread(fid,1,'double')';
info.('dCouchLat') = fread(fid,1,'double')';
info.('dIDUResolutionX') = fread(fid,1,'double')';
info.('dIDUResolutionY') = fread(fid,1,'double')';
info.('dImageResolutionX') = fread(fid,1,'double')';
info.('dImageResolutionY') = fread(fid,1,'double')';
info.('dEnergy') = fread(fid,1,'double')';
info.('dDoseRate') = fread(fid,1,'double')';
info.('dXRayKV') = fread(fid,1,'double')';
info.('dXRayMA') = fread(fid,1,'double')';
info.('dMetersetExposure') = fread(fid,1,'double')';
info.('dAcqAdjustment') = fread(fid,1,'double')';
info.('dCTProjectionAngle') = fread(fid,1,'double')';
info.('dCBCTPositiveAngle') = info.('dCTProjectionAngle') + 270.0;
info.('dCTNormChamber') = fread(fid,1,'double')';
info.('dGatingTimeTag') = fread(fid,1,'double')';
info.('dGating4DInfoX') = fread(fid,1,'double')';
info.('dGating4DInfoY') = fread(fid,1,'double')';
info.('dGating4DInfoZ') = fread(fid,1,'double')';
info.('dGating4DInfoTime') = fread(fid,1,'double')';
info.('dOffsetX') = fread(fid,1,'double');
info.('dOffsetY') = fread(fid,1,'double');
info.('dUnusedField') = fread(fid,1,'double');

if nargout == 1
    fclose(fid);
    return
end

%% Reading the projection

% Calculate bytes per pixels
nPixels = info.('uiSizeX') * info.('uiSizeY');

    % Check if the number of elements is correct
    if mod( info.('uiFileLength') - ftell(fid), nPixels) ~= 0
        fclose(fid);
        M = [];
        error('ERROR: Incompatible file format. \n');
        return
    end
    
bytesPerPixel = (info.('uiFileLength') - ftell(fid)) / nPixels;

% Read the projection according to the bytes per pixel information
switch bytesPerPixel
    case (4)
        M = single(fread(fid,nPixels,'single'));
    otherwise
        fclose(fid);
        M = [];
        error('ERROR : Incompatible file format. \n');
        return
end

% Reshape the projection data matrix
M = reshape(M,info.('uiSizeX'),info.('uiSizeY'));

% Close the file
fclose(fid);

return
