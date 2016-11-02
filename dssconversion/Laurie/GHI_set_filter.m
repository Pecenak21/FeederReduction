%% Program for reading in an irradiance curve generated from the GHI files
%% and interpolating the missing data
%% 
%% Input file is a .csv file with the irradiance curve in the fifth column
%% and 931 rows of data (first three columns are hours, minutes, seconds of
%% the timestamp, fourth column is the whole timestamp)
%% 
%% Output file is a .csv file containing the information from the input 
%% file and with the interpolated irradiance curve in the sixth column

strDirectory='C:\Users\lmiller\Documents\exmainpush\SDGE\GHI';
cd(strDirectory); 
number_PVs = 45;
filter_paramter = 0.3; % 0.4 for a little filtering, 0.1 for a lot

fileToRead1 = ['GHI_20091004_irradiance_set_interpolated.csv'];
fileToRead2 = ['520_PV.csv'];

if (exist (fileToRead2)) 
    PV_data = importdata(fileToRead2);
end
PV_loc = PV_data.data(:,2:3);

% Import the file
if (exist (fileToRead1)) 
    irradiance_interpolated = importdata(fileToRead1);
    
    for j = 1:number_PVs
        % put the (j+4)th column of rawData1 (containing the irradiance
        % values) into a vector:
        irradiance_partial = irradiance_interpolated(1:931, j+4);
        [b,a] = butter(1, filter_paramter, 'low');
        irradiance = filtfilt(b, a, irradiance_partial);
        % write the now complete irradiance vector to the (j+4)th column of
        % irradiance_interpolated:
        for i = 1:931
            irradiance_interpolated(i,j+4) = irradiance(i);
        end 
        %fig_title = PV_data.textdata{j+1,1};
        index = num2str(j);
        fig_title = [index '  ' PV_data.textdata{j+1,1}];
        figureFileName = fig_title;
        drawIrradianceCurves(irradiance, fig_title, figureFileName);
    end % for j = 1:number_PVs

end % if file exists    

% close figures
close all;

%% write peak irradiance values to .csv file
csvwrite('GHI_20091004_irradiance_set_filtered.csv',irradiance_interpolated);


