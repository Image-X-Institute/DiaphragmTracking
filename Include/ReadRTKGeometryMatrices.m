function result = ReadRTKGeometryMatrices(filename)
%% result = ReadRTKGeometryMatrices(filename)
% ------------------------------------------
% FILE   : ReadRTKGeometryMatrices.m
% AUTHOR : Andy Shieh, The University of Sydney
% DATE   : 2016-04-05  Created.
% ------------------------------------------
% PURPOSE
%   Read RTK geometry matrices from RTK geometry file.
% ------------------------------------------
% INPUT
%   filename : The path to RTK geometry .xml file.
% ------------------------------------------
% OUTPUT
%   result:    A 3 by 4 by N matrix storing N geometry matrices.
% ------------------------------------------

%% Checking input arguments & Opening the file

result = [];

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
    filename = fullfile(PathName, FileName);
    filename = strtrim(filename);
end

%% Read

fid = fopen(filename,'r');
k = 0;
line = fgetl(fid);
while(line~=-1)
    if strcmpi(strtrim(line),'<Matrix>');
        k = k+1;
        for n = 1:3;
            line = fgetl(fid);
            parts = regexp(strtrim(line),' ','split');
            m = 0;
            for p = 1:length(parts);
                if ~isempty(parts{p})
                    m = m +1;
                    result(n,m,k) = str2double(parts{p});
                end;
            end;
        end;
    end;
    line = fgetl(fid);
end;
fclose(fid);