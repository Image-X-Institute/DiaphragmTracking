function MhaWrite(info,M,outputfp)
%% MhaWrite(info,M,outputfp)
% ------------------------------------------
% FILE   : MhaWrite.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2013-01-09  Created.
% ------------------------------------------
% PURPOSE
%   Write the input header and 3D image to a Metaimage file. Support only
%   3D data.
% ------------------------------------------
% INPUT
%   info:       The header MATLAB struct. The template of the info struct:
%                                  Filename: 'SideADisBin01.mha'
%                                    Format: 'MHA'
%                            CompressedData: 'false'
%                                ObjectType: 'image'
%                        NumberOfDimensions: 3
%                                BinaryData: 'true'
%                                 ByteOrder: 'false'
%                           TransformMatrix: [1 0 0 0 1 0 0 0 1]
%                                    Offset: [-224.8400 -99 -224.8400]
%                          CenterOfRotation: [0 0 0]
%                     AnatomicalOrientation: 'RAI'
%                           PixelDimensions: [0.8800 2 0.8800]
%                                Dimensions: [512 100 512]
%                                  DataType: 'float'
%                                  DataFile: 'LOCAL'
%                                  BitDepth: 32
%                                HeaderSize: 318
%   M:          The image body, 3D.
%   outputfp:   The path of the output file.
% ------------------------------------------

%% Input argument check

% if (info.NumberOfDimensions ~= 3)
%     error('ERROR: Data is not 3-dimensional.');
%     return;
% end

if ~(strcmp(info.DataType,'float') || strcmp(info.DataType,'double') || ...
        strcmp(info.DataType,'int8') || strcmp(info.DataType,'char') || ...
        strcmp(info.DataType,'uint8') || strcmp(info.DataType,'uchar') || ...
        strcmp(info.DataType,'int16') || strcmp(info.DataType,'short') || ...
        strcmp(info.DataType,'uint16') || strcmp(info.DataType,'ushort') || ...
        strcmp(info.DataType,'int'))
    error('ERROR: Unsupported data type.');
    return;
end

if ~strcmp(info.CompressedData,'false')
    error('ERROR: MHA is an uncompressed format.');
end

if ~strcmp(info.DataFile,'LOCAL')
    error('ERROR: MHA is a single file format.');
end

fid = fopen(strtrim(outputfp),'w');
if(fid<=0) 
    error('ERROR: Invalid file path %s\n', outputfp);
    return;
end

%% Writing the header
fprintf(fid, 'ObjectType = Image\n');
fprintf(fid, num2str(info.NumberOfDimensions,'NDims = %d\\n'));
if isfield(info,'BinaryData')
    fprintf(fid, ['BinaryData = ',info.BinaryData,'\n']);
end;
fprintf(fid, ['BinaryDataByteOrderMSB = ',info.ByteOrder,'\n']);
fprintf(fid, ['CompressedData = ',info.CompressedData,'\n']);
if isfield(info,'TransformMatrix')
    fprintf(fid, ['TransformMatrix =',num2str(info.TransformMatrix,' %g'),'\n']);
end;
fprintf(fid, ['Offset =',num2str(info.Offset,' %g'),'\n']);
if isfield(info,'CenterOfRotation')
    fprintf(fid, ['CenterOfRotation =',num2str(info.CenterOfRotation,' %g'),'\n']);
end;
if isfield(info,'AnatomicalOrientation')
    fprintf(fid, ['AnatomicalOrientation = ',info.AnatomicalOrientation,'\n']);
end;
fprintf(fid, ['ElementSpacing =',num2str(info.PixelDimensions,' %g'),'\n']);
fprintf(fid, ['DimSize =',num2str(info.Dimensions,' %d'),'\n']);
if isfield(info,'ElementNumberOfChannels')
    fprintf(fid, ['ElementNumberOfChannels =',num2str(info.ElementNumberOfChannels,' %d'),'\n']);
end;
switch (info.DataType)
    case 'float'
        fprintf(fid, ['ElementType = MET_FLOAT\n']);
    case 'double'
        fprintf(fid, ['ElementType = MET_DOUBLE\n']);
    case {'uchar','uint8'}
        fprintf(fid, ['ElementType = MET_UCHAR\n']);
    case {'char','int8'}
        fprintf(fid, ['ElementType = MET_CHAR\n']);   
    case {'short','int16'}
        fprintf(fid, ['ElementType = MET_SHORT\n']);
    case {'ushort','uint16'}
        fprintf(fid, ['ElementType = MET_USHORT\n']);
    case 'int'
        fprintf(fid, ['ElementType = MET_INT\n']);
end
fprintf(fid, ['ElementDataFile = ',info.DataFile,'\n']);

%% Writing the body

% Need to reverse the dimension order if there are multi-entries
if isfield(info,'ElementNumberOfChannels') & info.ElementNumberOfChannels > 1
    M = reshape(M,[prod(info.Dimensions),str2double(info.ElementNumberOfChannels)]);
    M = M';
end;

switch(info.DataType)
    case {'char','int8'}
        fwrite(fid,M,'char');
    case {'uchar','uint8'}
        fwrite(fid,M,'uchar');
    case {'short','int16'}
        fwrite(fid,M,'short');
    case {'ushort','uint16'}
        fwrite(fid,M,'ushort');
    case 'int'
        fwrite(fid,M,'int');
    case 'uint'
        fwrite(fid,M,'uint');
    case 'float'
        fwrite(fid,M,'float');
    case 'double'
        fwrite(fid,M,'double');
end

fclose(fid);

return
