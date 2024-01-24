function [info M] = MhaRead(filename)
%% [info M] = MhaRead(filename)
% ------------------------------------------
% FILE   : MhaRead.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2012-11-08  Created.
% ------------------------------------------
% PURPOSE
%   Read metaimage mha header and image.
% ------------------------------------------
% INPUT
%   filename : The full path to the file.
% ------------------------------------------
% OUTPUT
%   info  : header information in a struct.
%   M     : the image stored in a 3D matrix.
% ------------------------------------------

%% Checking input arguments & Opening the file

M = [];

% If no input filename => Open file-open-dialog
if nargin < 1
    
    % go into a default directory
    DefaultDir = pwd;
    
    % get the input file (hnc) & extract path & base name
    [FileName,PathName] = uigetfile( {'*.mha','MetaImage (*.mha)';}, ...
        'Select an image file', ...
        DefaultDir);
    
    % catch error if no file selected
    if isnumeric(FileName)
        info = struct([]);
        error('ERROR: No file selected.');
    end
    
    % make same format as input
    filename = fullfile(PathName, FileName);
    
end
filename = strtrim(filename);

%% Reading the header and body using external module
info = mha_read_header(filename);

if nargout == 1
    return;
else
    M = mha_read_volume(info);
end

return
