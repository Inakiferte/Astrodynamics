function [d,m,s] = radians_to_dms(radians)
degrees = rad2deg(radians);
d = round(degrees);
m = abs(round((degrees - d) * 60));
s = abs((degrees - d - m / 60) * 3600);
end