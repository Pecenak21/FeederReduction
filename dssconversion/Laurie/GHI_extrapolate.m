%% Program for reading in the interpolated irradiance curve generated from 
%% the GHI files and extrapolating around the data to the zero crossings 
%% to extend the curve into morning and evening times when there is no data
%% 
%% Input file is a .csv file with the interpolated irradiance curve in 
%% the sixth column and 931 rows of data (first three columns are hours, 
%% minutes, seconds of the timestamp, fourth column is the whole timestamp,
%% fifth column is the uninterpolated irradiance data)
%% 
%% Output file is a .csv file containing the information from the input 
%% file and with the extrapolated irradiance curve in the seventh column

strDirectory='C:\Users\lmiller\Documents\exmainpush\SDGE\GHI';
cd(strDirectory); 
formatspec = '%02i';

fileToRead1 = ['GHI_20091004_irradiance_interpolated.csv'];

% Import the file
if (exist (fileToRead1)) 
    irradiance_interpolated = importdata(fileToRead1);

    %% create a matrix for the early morning hours
    % hour 7
    early_morning = zeros(121,6);
    n_hours = 8;
    n_minutes = -1;
    for i = 1:120
        early_morning(i,1)= n_hours;
        if mod(i,2) == 1
            n_seconds = 0;
            n_minutes = n_minutes + 1;
        else
            n_seconds = 30;
        end
        early_morning(i,2)= n_minutes;
        early_morning(i,3)= n_seconds;
        hours = num2str(n_hours, formatspec);
        minutes = num2str(n_minutes, formatspec);
        seconds = num2str(n_seconds, formatspec);
        early_morning(i,4) = str2num([hours minutes seconds]);
    end
    % the first half of the first minute of hour nine:
    n_hours = 9;
    n_minutes = 0;
    n_seconds = 0;
    early_morning(121,1) = n_hours;
    early_morning(121,2)= n_minutes;
    early_morning(121,3)= n_seconds;
    hours = num2str(n_hours, formatspec);
    minutes = num2str(n_minutes, formatspec);
    seconds = num2str(n_seconds, formatspec);
    early_morning(121,4) = str2num([hours minutes seconds]);

    %% create a matrix for the evening hours
    % 16:46:00-16:59:30
    evening_16 = zeros(28,6);
    n_hours = 16;
    n_minutes = 45;
    for i = 1:28
        evening_16(i,1)= n_hours;
        if mod(i,2) == 1
            n_seconds = 0;
            n_minutes = n_minutes + 1;
        else
            n_seconds = 30;
        end
        evening_16(i,2)= n_minutes;
        evening_16(i,3)= n_seconds;
        hours = num2str(n_hours, formatspec);
        minutes = num2str(n_minutes, formatspec);
        seconds = num2str(n_seconds, formatspec);
        evening_16(i,4) = str2num([hours minutes seconds]);
    end
    % 17:00:00 - 19:59:30
    evening_17_20 = zeros(360,6);
        for h = 1:3
        if h == 1
            n_hours = 17;
        elseif h == 2
            n_hours = 18;
        else
            n_hours = 19;
        end
        n_minutes = -1;
        for i = 1:120
            evening_17_20(i+(h-1)*120,1)= n_hours;
            if mod(i,2) == 1
                n_seconds = 0;
                n_minutes = n_minutes + 1;
            else
                n_seconds = 30;
            end
            evening_17_20(i+(h-1)*120,2)= n_minutes;
            evening_17_20(i+(h-1)*120,3)= n_seconds;
            hours = num2str(n_hours, formatspec);
            minutes = num2str(n_minutes, formatspec);
            seconds = num2str(n_seconds, formatspec);
            evening_17_20(i+(h-1)*120,4) = str2num([hours minutes seconds]);
        end
    end %for h = 1:3
    evening = [evening_16; evening_17_20];
    
    % concatenate the early_morning and evening matrices onto the 
    % irradiance_interpolated matrix
    irradiance_extrapolated = [early_morning; irradiance_interpolated; evening]; 
    
    % copy the existing data into the 7th column:
    for i = 1:1440
        irradiance_extrapolated(i,7) = irradiance_extrapolated(i,6);
    end
    
    % extrapolation for the morning
    % put data the sixth column of irradiance_extrapolated (containing the interpolated 
    % irradiance values) into a vector:
    morning_extrapolation = irradiance_extrapolated(1:300, 6);
    % set the extrapolation point for the morning 
    morning_extrapolation(6) = 0.1;
    % create an index vector from 1..300:
    for i = 1:300
        x_morning(i) = i;
    end
    % create list of indices of timesteps that have data
    % timesteps with no data are skipped
    xi_morning=x_morning(find(morning_extrapolation));
    % create list of irradiance values of timesteps that have data
    % timesteps with no data are skipped
    yi_morning=morning_extrapolation(find(morning_extrapolation));   
    % interpolate the missing values:
    irradiance_morning=interp1(xi_morning,yi_morning,x_morning,'linear','extrap');
    irradiance_morning(irradiance_morning<0)=0;
    % write the new irradiance vector to the 7th column of
    % irradiance_extrapolated:
    for i = 1:121
        irradiance_extrapolated(i,7) = irradiance_morning(i);
    end

    % extrapolation for the evening
    % put data the sixth column of irradiance_extrapolated (containing the interpolated 
    % irradiance values) into a vector:
    evening_extrapolation = irradiance_extrapolated(801:1440, 6);
    % set the extrapolation point for the evening 
    %evening_extrapolation(620) = 0.1;
    % create an index vector from 1..640:
    for i = 1:640
        x_evening(i) = i;
    end
    % create list of indices of timesteps that have data
    % timesteps with no data are skipped
    xi_evening=x_evening(find(evening_extrapolation));
    % create list of irradiance values of timesteps that have data
    % timesteps with no data are skipped
    yi_evening=evening_extrapolation(find(evening_extrapolation));   
    % interpolate the missing values:
    irradiance_evening=interp1(xi_evening,yi_evening,x_evening,'linear', 'extrap');
    irradiance_evening(irradiance_evening<0)=0;
    % write the new irradiance vector to the 7th column of
    % irradiance_extrapolated:
    for i = 1:388
        irradiance_extrapolated(i+1052,7) = irradiance_evening(i+252);
    end
   
end % if file exists    

% clear files from workspace
%clear rawData1;

%% write peak irradiance values to .csv file
csvwrite('GHI_20091004_irradiance_extrapolated.csv',irradiance_extrapolated);


