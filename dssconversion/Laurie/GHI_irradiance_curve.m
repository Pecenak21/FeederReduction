%% Program for reading in GHI files and outputing an irradiance curve in
%% .csv
%%
%% Output file is a .csv file with the irradiance curve in the fifth column
%% and 931 rows of data (first three columns are hours, minutes, seconds of
%% the timestamp, fourth column is the whole timestamp)


strDirectory='C:\Users\lmiller\Documents\exmainpush\SDGE\GHI';
cd(strDirectory); 

i = 1;
irradiance = zeros(931,5);
formatspec = '%02i';

% construct the filename:
n_hours = 09;
n_minutes = 00;
seconds = '00';
while ((n_hours < 16)||(n_hours == 16 && n_minutes < 46))
%while ((n_minutes < 4))
    if seconds == '00'
        seconds = '30';
    elseif n_minutes == 59
        n_hours = n_hours + 1;
        n_minutes = 0;
        seconds = '00';
    else
        n_minutes = n_minutes + 1;
        seconds = '00';
    end
    hours = num2str(n_hours, formatspec);
    minutes = num2str(n_minutes, formatspec);
    %time = str2num([hours minutes seconds]);

    fileToRead1 = ['GHI_20091004' hours minutes seconds '.csv'];
  
    % Import the file
    if (exist (fileToRead1)) 
        rawData1 = importdata(fileToRead1);
        rawData1(isnan(rawData1))= 0;

    %% find maximum value in this timeslice
        irradiance(i,5) = max(max(rawData1));
    end % if file exists
    
    % record timestamp information to .csv file
    irradiance(i,1) = n_hours;
    irradiance(i,2) = n_minutes;
    if seconds == '30'
        irradiance(i,3) = 30;
    end
    irradiance(i,4) = str2num([hours minutes seconds]);
    i = i + 1  
    
    % clear files from workspace
    clear rawData1;

end % while ((n_hours < 16)||(n_hours = 16 && n_minutes < 46))

irradiance_trimmed = irradiance( [1:931] , : )

%% write peak irradiance values to .csv file
csvwrite('GHI_20091004_irradiance_incomplete.csv',irradiance_trimmed);


