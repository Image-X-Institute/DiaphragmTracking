function [info,M] = HisRead(filename)
%% [info M] = HisRead(filename)
% ------------------------------------------
% FILE   : HisRead.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2014-06-02  Created.
% ------------------------------------------
% PURPOSE
%   Read his header and image.
% ------------------------------------------
% INPUT
%   filename : The full path to the hnc file.
% ------------------------------------------
% OUTPUT
%   info  : Header information (in uint16 un-decoded)
%   M     : The image stored in a 2D matrix (uint16).
% ------------------------------------------

%% Checking input arguments & Opening the file

M = [];

HISHeadLength = 68;
ElektaDetSizeX = 409.6;
ElektaDetSizeY = 409.6;

% If no input filename => Open file-open-dialog
if nargin < 1
    
    % go into a default directory
    DefaultDir = pwd;
    
    % get the input file (hnc) & extract path & base name
    [FileName,PathName] = uigetfile( {'*.his;','Elekta Image Files (*.his)'}, ...
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

%% Reading the header

fileinfo = dir(filename);
info.('uiActualFileLength') = fileinfo.('bytes');

header = fread(fid,HISHeadLength,'int8');
if (header(1) ~=0 || header(2) ~= 112 || header(3) ~= 68 || header(4) ~=0)
    error(['ERROR: File:',filename,' is not in Heimann HIS format version 100']);
end

info.('HeaderSize') = header(11) + bitshift(header(12),8) + HISHeadLength;
ulx = header(13) + bitshift(header(14),8);
uly = header(15) + bitshift(header(16),8);
brx = header(17) + bitshift(header(18),8);
bry = header(19) + bitshift(header(20),8);
info.('NFrames') = header(21) + bitshift(header(22),8);
info.('Type') = header(33) + bitshift(header(35),8);

info.('SizeX') = bry-uly+1;
info.('SizeY') = brx-ulx+1;
info.('PixelSpacingX') = ElektaDetSizeX / info.SizeX;
info.('PixelSpacingY') = ElektaDetSizeY / info.SizeY;
info.('OriginX') = -0.5 * (info.('SizeX') - 1) * info.('PixelSpacingX');
info.('OriginY') = -0.5 * (info.('SizeY') - 1) * info.('PixelSpacingY');

if nargout == 1
    fclose(fid);
    return
end

%% Reading the projection

% Calculate bytes per pixels
nPixels = info.('SizeX') * info.('SizeY') * info.('NFrames');

% Read the projection according to the bytes per pixel information
fseek(fid,info.HeaderSize,-1);
switch info.('Type')
    case (4)
        M = uint16(fread(fid,nPixels,'uint16'));
    otherwise
        M = uint16(fread(fid,nPixels,'uint16'));
end

% Reshape the projection data matrix
M = reshape(M,info.('SizeX'),info.('SizeY'),info.('NFrames'));

% Close the file
fclose(fid);

return
