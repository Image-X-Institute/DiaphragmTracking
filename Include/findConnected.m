function result = findConnected(M,minN)
%% result = findConnected(M,minN)
% ------------------------------------------
% FILE   : findConnected.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2013-12-10   Created.
% ------------------------------------------
% PURPOSE
%   Find the maximum connected region of the input image, or excluded
%   regions with less than minN connected pixels.
%   Support only logical input for now.
%
% ------------------------------------------
% INPUT
%
%  M
%   The input logical map.
%
%  minN
%   (Optional) Regions with less than minN connected pixels will be
%   excluded.
%   This field is optional. If it is not input, or input as "0", the
%   function will look for the largest connected region.
%   (default: 0, where findConnected looks for the largest connected
%   region)
%
% OUTPUT
%
%  result
%   The resultant logical map with the unconnected region eliminated.

%% Using MATLAB built-in function bwconncomp

if ~islogical(M)
    error('ERROR: The input M must be in logical format.');
end;

if nargin < 2
    minN = 0;
else
    if ~isnumeric(minN) || ~isscalar(minN) || minN < 0 || mod(minN,1) ~= 0
        error('ERROR: The input minN must be zero or a positive integer.');
    end
end;

conn = bwconncomp(M);
nPixels = cellfun(@numel,conn.PixelIdxList);

if minN == 0
    result = false(size(M));
    [~,idx] = max(nPixels);
    result(conn.PixelIdxList{idx}) = true;
else
    result = M;
    noise_index = find(nPixels < minN);
    for k = 1:length(noise_index);
        result(conn.PixelIdxList{noise_index(k)}) = false;
    end;
end;

return;