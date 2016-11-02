%% Program for reading in GHI files and outputing a set of irradiance 
%% curves in .csv
%%
%% Output file is a .csv file with the irradiance curve in the fifth column
%% and 931 rows of data (first three columns are hours, minutes, seconds of
%% the timestamp, fourth column is the whole timestamp)


strDirectory='C:\Users\lmiller\Documents\exmainpush\SDGE\GHI';
cd(strDirectory); 

i = 1;
irradiance = zeros(931,49);
formatspec = '%02i';
x_offset = 297;
y_offset = 297;
number_PVs = 45;

fileToRead1 = ['520_PV.csv'];

% Import the file
if (exist (fileToRead1)) 
    PV_data = importdata(fileToRead1);
end
PV_loc = PV_data.data(:,2:3);
% 1 kft = 305 meters, 1 data point = 60 meters
PV_loc = PV_loc.*5.083;
% find minimum y value
y_min = min(PV_loc(:,1));
% decrement each y by y_min
PV_loc(:,1)=PV_loc(:,1)-y_min;
% find minimum x value
x_min = min(PV_loc(:,2));
% decrement each x by x_min
PV_loc(:,2)=PV_loc(:,2)-x_min;

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
        %rawData1(isnan(rawData1))= 0;
        
        % get the irradiance values
        for j = 1:number_PVs
            x = x_offset + round(PV_loc(j,2));
            y = y_offset + round(PV_loc(j,1));
            irradiance(i,j+4) = rawData1(x,y);
        end
    end % if file exists
    
    % record timestamp information 
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

irradiance_trimmed = irradiance( [1:931] , : );

%% write peak irradiance values to .csv file
csvwrite('GHI_20091004_irradiance_set_incomplete.csv',irradiance_trimmed);


