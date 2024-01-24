function [info M] = HndRead(filename)
%% [info M] = HndRead(filename)
% ------------------------------------------
% FILE   : HndRead.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2016-01-10  Created.
% ------------------------------------------
% PURPOSE
%   Read Varian hnd header and image.
% ------------------------------------------
% INPUT
%   filename : The full path to the hnd file.
% ------------------------------------------
% OUTPUT
%   info  : Header information in a struct.
%   M     : The image stored in a 2D matrix.
% ------------------------------------------

%% Checking input arguments & Opening the file

M = [];

% If no input filename => Open file-open-dialog
if nargin < 1
    
    % go into a default directory
    DefaultDir = pwd;
    
    % get the input file (hnc) & extract path & base name
    [FileName,PathName] = uigetfile( {'*.hnd;*.hnc;','Varian Image Files (*.hnd,*.hnc,*.mat,*.mdl)';
        '*.hnc',  'Portal Images (*.hnc)'; ...
        '*.hnd',  'OBI kv Images (*.hnd)'}, ...
        'Select an image file', ...
        DefaultDir);
    
    % catch error if no file selected
    if isnumeric(FileName)
        M = []; info = struct([]);
        error('ERROR: No file selected. \n');
    end
    
    % make same format as input
    filename = fullfile(PathName, FileName); 
end

filename = strtrim(filename);

% Open the file
fid = fopen(filename,'r');

    % catch error if failure to open the file
    if fid == -1
        info = struct([]);
        error('ERROR : Failure in opening the file. \n');
    end

%% Reading the header (same for hnc and hnd)

info.('bFileType') = fread(fid,32,'uint8=>char')';
info.('uiFileLength') = fread(fid,1,'uint32')';

% Sometimes the actual file length is different from the one specified in
% the header
fileinfo = dir(filename);
info.('uiActualFileLength') = fileinfo.('bytes');

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
info.('dCBCTPositiveAngle') = mod(info.('dCTProjectionAngle') + 270.0, 360.0);
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

fclose(fid);

%% Use Varian Image toolbox to read the image body
javaaddpath(which('VarianReader.jar'));

readerObj = CPS_Reader();
M = readerObj.getHNC_HND(filename);
clear readerObj;
javarmpath(which('VarianReader.jar'));

% Reshape the projection data matrix
M = reshape(M,info.('uiSizeX'),info.('uiSizeY'));

return
