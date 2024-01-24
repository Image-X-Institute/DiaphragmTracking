function angles = ReadAnglesFromRTKGeometry(rtkFile)
%% angles = ReadAnglesFromRTKGeometry(rtkFile)
% ------------------------------------------
% FILE   : ReadAnglesFromRTKGeometry.m
% AUTHOR : Andy Shieh, The University of Sydney
% DATE   : 2016-03-22 Created.
% ------------------------------------------
% PURPOSE
%   Get a list of angles (within 0 to 360 degrees) from RTK geometry file.
%
% ------------------------------------------
% INPUT
%   rtkFile:       The RTK geometry file path.
% ------------------------------------------
% OUTPUT
%   angles:         The list of angles in degrees.

%% 
% If no input filename => Open file-open-dialog
if nargin < 1
    
    % go into a default directory
    DefaultDir = pwd;
    
    % get the input file (hnc) & extract path & base name
    [FileName,PathName] = uigetfile( {'*.xml;*.XML;','RTK geometry file (*.xml,*XML)'}, ...
        'Select an RTK geometry file', ...
        DefaultDir);
    
    % catch error if no file selected
    if isnumeric(FileName)
        result = [];
        error('ERROR: No file selected. \n');
    end
    
    % make same format as input
    rtkFile = fullfile(PathName, FileName);
    rtkFile = strtrim(rtkFile);
end


fid = fopen(rtkFile,'r');
line = fgetl(fid); k = 0;
while(line~=-1);
    if(strfind(line,'<GantryAngle>'));
        k = k + 1;
        line = regexp(line,'<GantryAngle>','split');
        line = regexp(line{2},'</GantryAngle>','split');
        angles(k) = str2double(line{1});
    end;
    line = fgetl(fid);
end;
fclose(fid);

angles = mod(angles,360);

return;