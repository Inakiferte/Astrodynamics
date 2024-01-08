function [jul_day,jul_cent] = Julian(y,m,d,utch,utcm,utcs)
%=======================================
% This function computes the julian day
% and the julian century.
%
% INPUTS:
%       - y    : year.
%       - m    : month.
%       - d    : day.
%       - utch : UTC hour.
%       - utcm : UTC minute.
%       - utcs : UTC seconds.
% OUTPUTS:
%       - jul_day  : julian day.
%       - jul_cent : julian century.
%=======================================

    % Compute j0
    j0 = 367 * y - fix(7 * (y + fix((m + 9) / 12)) / 4) + fix(275 * m / 9) + d + 1721013.5;

    % Compute UT
    UT = utch + (utcm / 60.0) + (utcs / 3600.0);

    % Compute julian day
    jul_day = j0 + (UT / 24.0);
   
    % Compute julian century
    jul_cent = (jul_day - 2451545) / 36525.0;
end