function PCP_plot(departure_planet, arrival_planet, d_Days, tof_Days, jul_day_D, deltaV, Mode,tof)
%=======================================
% Function that creates the PCP and 
% saves deltaV matrix.
%
% INPUTS:
%       - departure_planet : Departure 
%                            planet index.
%       - arrival_planet   : Arrival 
%                            planet index.
%       - d_Days           : Departure 
%                            maximun days.
%       - tof_Days         : Time of flight 
%                            maximun days.
%       - juk_day_D        : Departure
%                            julian day.
%       - deltaV           : Computed
%                            deltaV matrix.
%       - Mode             : PCP mode.
%       - tof              : Time of flight 
%                            array.
% OUTPUT: It creates and saves the PCP
%=======================================
    % Convert axis
    date = zeros(1, d_Days);
    for i = 1:d_Days
        date(i) = jul_day_D + i - jul_day_D;
    end
    
    for i = 1:tof_Days
        tof(i) = tof(i) / 86400;
    end

    % Set values exceeding 50 to NaN
    deltaV(deltaV > 50) = NaN;

    % Create colormap plot using pcolor
    [X, Y] = meshgrid(tof, date);

    figure;
    pcolor(X, Y, deltaV);
    shading interp; % Optional: interpolate colors for a smoother appearance
    colormap('jet');
    c = colorbar;
    ylabel(c, '$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Optionally, you can set axis labels
    xlabel('TOF [days]', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('Departure Days From 01-01-2020', 'Interpreter', 'latex', 'FontSize', 20);

    % Title for the plot
    planet_names = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto (extra)'};
    departure_name = planet_names{departure_planet};
    arrival_name = planet_names{arrival_planet};
    title_str = sprintf('PCP: From %s to %s', departure_name, arrival_name);
    title(title_str, 'Interpreter', 'latex', 'FontSize', 20);

    axis([2500 tof_Days 0 d_Days]);

    %Add grid lines if needed
    grid on;

    % Folder to save all the plots
    folder_name = 'PCP_Plots';
    folder_name2 = 'deltaV';

    % Save the plot in the specified folder with a dynamic filename
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end
    if ~exist(folder_name2, 'dir')
        mkdir(folder_name2);
    end

    % Convert Mode to a string
    Mode_str = num2str(Mode);

    % Generate a filename based on the title
    filename = sprintf('PCP_%s_to_%s_Mode_%s.png', departure_name, arrival_name, Mode_str);
    saveas(gcf, fullfile(folder_name, filename));

    % Save deltaV
    filename2 = sprintf('deltaV_%s_to_%s_Mode_%s.mat', departure_name, arrival_name, Mode_str);
    save(fullfile(folder_name2, filename2),'deltaV')

    % Optionally, you can close the figure after saving
    close(gcf);
    
    disp(" ");
    disp('Plot has been successfully saved!');
    disp('You can find the plot and deltaV matrix in the following path: ')
    disp([folder_name '/' filename])
    disp([folder_name2 '/' filename2])
end
