function [jul_day,jul_cent] = Julian(y,m,d,utch,utcm,utcs)
    % Compute j0
    j0 = 367 * y - round(7 * (y + round((m + 9) / 12)) / 4) + round(275 * m / 9) + d + 1721013.5;

    % Compute UT
    UT = utch + (utcm / 60.0) + (utcs / 3600.0);

    % Compute julian day
    jul_day = j0 + (UT / 24.0);

    jul_cent = (jul_day - 2451545) / 36525;
end