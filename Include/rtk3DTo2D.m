function [twoDPos, rayScaleFactors] = rtk3DTo2D(matrices,threeDPos)
%% [twoDPos, rayScaleFactors] = rtk3DTo2D(matrices,threeDPos)
% ------------------------------------------
% FILE   : rtk3DTo2D.m
% AUTHOR : Andy Shieh, Sydney Medical Schoo, The University of Sydney
% DATE   : 2016-04-05   Created.
% ------------------------------------------
% PURPOSE
%   Convert one or multiple 3D positions to detector 2D positions (in mm) 
%   using RTK geomery matrices.
%
% ------------------------------------------
% INPUT
%   matrix :    The RTK geometry matrix or matrices (3 by 4 by N).
%   threeDPos:  A N by 3 vector representing the 3D position to be
%               converted.
% ------------------------------------------
% OUTPUT
%   twoDPos  :  The corresponding 2D positions in detector space (in mm,
%               0,0 corresponds to detector center)
%   rayScaleFactors: The distance dependent ray scaling factors used in RTK.
% ------------------------------------------

%% 
rayScaleFactors = ...
    squeeze(matrices(3,1,:)) .* threeDPos(:,1) + ...
    squeeze(matrices(3,2,:)) .* threeDPos(:,2) + ...
    squeeze(matrices(3,3,:)) .* threeDPos(:,3) + ...
    squeeze(matrices(3,4,:));
twoDPos(:,1) = ...
    ( squeeze(matrices(1,1,:)) .* threeDPos(:,1) + ...
    squeeze(matrices(1,2,:)) .* threeDPos(:,2) + ...
    squeeze(matrices(1,3,:)) .* threeDPos(:,3) + ...
    squeeze(matrices(1,4,:)) ) ./ rayScaleFactors;
twoDPos(:,2) = ...
    ( squeeze(matrices(2,1,:)) .* threeDPos(:,1) + ...
    squeeze(matrices(2,2,:)) .* threeDPos(:,2) + ...
    squeeze(matrices(2,3,:)) .* threeDPos(:,3) + ...
    squeeze(matrices(2,4,:)) ) ./ rayScaleFactors;

return;
