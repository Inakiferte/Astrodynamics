function loading_bar(i,d_Days,departure_planet,arrival_planet)
%=======================================
% Displays the loading text in the 
% command window.
% 
% INPUTS: 
%       - i : iteration
%       - d_Days           : Departure 
%                            maximun days.
%       - departure_planet : Departure 
%                            planet index.
%       - arrival_planet   : Arrival 
%                            planet index.
%=======================================

planet_names = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto (extra)'};
departure_name = planet_names{departure_planet};
arrival_name = planet_names{arrival_planet};

if i==1
disp( "=======================================================");
disp( "Interplanetary Mission Analysis Design Tool.");
disp( "Master in Space and Aeronautical Engineering.");
disp( "Astrodynamics Course, Project.");
disp( "By: Jorge Simón, Javier Sanchez & Iñaki Fernandez.");
disp( "=======================================================");
disp(" ");
disp("Code has started runing. It might take a few minutes ;)")
disp(['Departure planet is : ' , departure_name]);
disp(['Arrival planet is: ', arrival_name]);

elseif i == fix(d_Days / 4)
disp(" ")
disp("25 % of the code has been done")

elseif i == fix(d_Days / 2)
disp(" ")
disp("50 % of the code has been done")

elseif i == fix((d_Days / 4) + (d_Days / 2))
disp(" ")
disp("75 % of the code has been done")

elseif i == d_Days
disp(" ")
disp("Code has finished succesfully!")
disp(" ")
end

end