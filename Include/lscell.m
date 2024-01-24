function fl_out = lscell(rexp,dirmode)
%% fl_out = lscell(rexp,dirmode)
% ------------------------------------------
% FILE   : lscell.m
% AUTHOR : Andy Shieh, The University of Sydney
% DATE   : 2016-04-08   Creation
% ------------------------------------------
% PURPOSE
%   Same as lscp.m but output file list in cells.
%
% ------------------------------------------
% INPUT
%   rexp:         The regular expression. rexp can contain file path as
%                 well. e.g. '\Projection\bin01\SPI*'
%   dirmode:      If the user input 'dir' for the second input, then only
%                 directories will be considered.
% ------------------------------------------
% OUTPUT
%   fl:             The file list in cells.

%%
if nargin < 2;  dirmode = '';   end;

[dirname,rexp,ext] = fileparts(rexp);
rexp = [rexp,ext];

currentDir = pwd;
if ~isempty(dirname);
    try
        cd (dirname);
    catch
        fl_out = {};
        cd (currentDir);
        return;
    end
end;

fl = '';
if ispc
    fl = ls(rexp);
elseif isunix || ismac
    fl_s = dir(rexp);
    n = length(fl_s);
    for k = 1:n
        l = length(fl_s(k).name);
        fl(k,1:l) = fl_s(k).name;
    end
else
    fl = ls(rexp);
end

if isempty(fl)
    fl_out = {};
    cd (currentDir);
    return;
end;

% Exclude . and ..
ind = find(fl(:,1)=='.');
fl(ind,:) = [];

if isempty(fl)
    cd (currentDir);
    fl_out = {};
    return;
end;

if strcmp(dirmode,'dir')
    count = 0;
    for k = 1:size(fl,1)
        if ~isdir(fl(k-count,:))
            fl(k-count,:) = [];
            count = count + 1;
        end
    end;
end;    

% Add in the directory path
% fl = [ char(ones(size(fl,1),1) * [pwd,fsepcp]) , fl];
if size(fl,1) > 0
    for k = 1:size(fl,1)
        fl_out{k} = fullfile(pwd,strtrim(fl(k,:)));
    end;
else
    fl_out = {};
end;

if ~isempty(dirname);
    cd (currentDir);
end;

return;
