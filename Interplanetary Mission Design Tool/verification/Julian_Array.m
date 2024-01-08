function [jul_cent_D, jul_cent_A, tof] = Julian_Array(jul_day_D, tof_Days)
%=======================================
% We compute the required time arrays
% for PCP
%=======================================

tof = zeros(1,tof_Days);
jul_cent_A = zeros(1,tof_Days);

% Compute tof days
for i=1:tof_Days
    tof(i) = i * 86400;
end

% Julian centuries
jul_cent_D = (jul_day_D - 2451545) / 36525;

for i=1:tof_Days
    jul_cent_A(i) = ((jul_day_D + i) - 2451545) / 36525;
end

end