function [info M] = ProjRead(filename)
%% [info M] = ProjRead(filename)
% ------------------------------------------
% FILE   : ProjRead.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2014-06-10  Created.
% ------------------------------------------
% PURPOSE
%   Read the header and image body of hnc, his, or att projection format.
% ------------------------------------------
% INPUT
%   filename : The full path to the file.
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
    [FileName,PathName] = uigetfile( {'*.att','Attenuation Image Files (*.att)'}, ...
        'Select an image file', ...
        DefaultDir);
    
    % catch error if no file selected
    if isnumeric(FileName)
        error('ERROR: No file selected. \n');
        return
    end
    
    % make same format as input
    filename = fullfile(PathName, FileName);
    
end
    
%% Reading

[~,~,ext] = fileparts(filename);

if strcmpi(ext(2:end),'hnc')
    [info,M] = HncRead(filename);
elseif strcmpi(ext(2:end),'hnd')
    [info,M] = HndRead(filename);
elseif strcmpi(ext(2:end),'att')
    [info,M] = AttRead(filename);
elseif strcmpi(ext(2:end),'his')
    [info,M] = HisRead(filename);
else
    error('ERROR: Unsupported projection format.');
end

return;
