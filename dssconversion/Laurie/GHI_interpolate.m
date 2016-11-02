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


    fileToRead1 = ['GHI_20091004_irradiance_incomplete.csv'];
  
    % Import the file
    if (exist (fileToRead1)) 
        irradiance_interpolated = importdata(fileToRead1);
        
        % put the fourth column of rawData1 (containing the irradiance
        % values) into a vector:
        irradiance_partial = irradiance_interpolated(1:931, 5);
        % create an index vector from 1..931:
        for i = 1:931
            x(i) = i;
        end
        % create list of indices of timesteps that have data
        % timesteps with no data are skipped
        xi=x(find(irradiance_partial));
        % create list of irradiance values of timesteps that have data
        % timesteps with no data are skipped
        yi=irradiance_partial(find(irradiance_partial));
        % interpolate the missing values:
        irradiance=interp1(xi,yi,x,'spline');
        % write the now complete irradiance vector to the 5th column of
        % irradiance_interpolated:
        for i = 1:931
            irradiance_interpolated(i,6) = irradiance(i);
        end

    end % if file exists    

    % clear files from workspace
    %clear rawData1;

    %% write peak irradiance values to .csv file
    csvwrite('GHI_20091004_irradiance_interpolated.csv',irradiance_interpolated);


