function [hours,minutes,seconds] = degrees_to_time(degrees)
    total_hours = degrees * 24 / 360;
    hours = round(total_hours);
    remaining_minutes = (total_hours - hours) * 60;
    minutes = round(remaining_minutes);
    seconds = abs((remaining_minutes - minutes) * 60);
end