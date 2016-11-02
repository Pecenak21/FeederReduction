%% Program for computing offsets for each timestep of the GHI data
%% Offsets are to be used when employing the "chase the data around the
%% square" option
%%
%% Output: csv file with a row per timestep, each row containing 
%% timestamp
%% first_row  = first row that has data
%% last_row  = last row that has data
%% first_column  = first column that has data
%% last_column = last column that has data
%% x_radius = difference between first and last columns
%% y_radius = difference between first and last rows
%% x_offset = midpoint between first and last columns
%% y_offset = midpoint between first and last rows

% not done yet -- if any of the first/last rows/columns are at an edge, assume data field
% is not round

strDirectory='C:\Users\lmiller\Documents\exmainpush\SDGE\GHI';
cd(strDirectory); 

GHI_offset_matrix = zeros(931,9);
i = 1;

% construct the filename:
n_hours = 09;
n_minutes = 00;
seconds = '00';
while ((n_hours < 16)||(n_hours == 16 && n_minutes < 46))
%while ((n_minutes < 4))
    if seconds == '00';
        seconds = '30';
    elseif n_minutes == 59;
        n_hours = n_hours + 1;
        n_minutes = 0;
        seconds = '00';
    else
        n_minutes = n_minutes + 1;
        seconds = '00';
    end
    formatspec = '%02i';
    hours = num2str(n_hours, formatspec);
    minutes = num2str(n_minutes, formatspec);

    fileToRead1 = ['GHI_20091004' hours minutes seconds '.csv'];
  
    % Import the file
    if (exist (fileToRead1)) 
        rawData1 = importdata(fileToRead1);
        rawData1(isnan(rawData1))= 0;

        % create a true_false matrix of rawData1
        true_false = rawData1==0;
        % create a true_false vector of which rows have data
        true_false_rows = ~all(true_false,2);
        % find the first and last rows that have data
        first_row = find(true_false_rows, 1, 'first');
        last_row = find(true_false_rows, 1, 'last');
        % calculate the height of the data area
        y_radius = last_row - first_row;
        % calculate the middle row that has data to get the y-offset
        y_offset = floor(y_radius/2);
        
        % create a true_false vector of which columns have data
        true_false_columns = ~all(true_false,1);
        % find the first and last columns that have data
        first_column = find(true_false_columns, 1, 'first');
        last_column = find(true_false_columns, 1, 'last');
        % calculate the width of the data area
        x_radius = last_column - first_column;
        % calculate the middle column that has data to get the x-offset
        x_offset = floor(x_radius/2);
        
        % record:
        % timestamp
        GHI_offset_matrix(i,1) =  str2num([hours minutes seconds]);
        % first_row  = first row that has data
        GHI_offset_matrix(i,2) = first_row;
        % last_row  = last row that has data
        GHI_offset_matrix(i,3) = last_row;
        % first_column  = first column that has data
        GHI_offset_matrix(i,4) = first_column;
        % last_column = last column that has data
        GHI_offset_matrix(i,5) = last_column;
        % x_radius = difference between first and last columns
        GHI_offset_matrix(i,6) = x_radius;
        % y_radius = difference between first and last rows
        GHI_offset_matrix(i,7) = y_radius;
        % x_offset = midpoint between first and last columns
        GHI_offset_matrix(i,8) = x_offset;
        % y_offset = midpoint between first and last rows
        GHI_offset_matrix(i,9) = y_offset;

    end % if file exists
  
    i = i+1;
        
    % clear data from workspace
    clear rawData1;
    clear true_false;
    % close figure
    %close all;

end % while ((n_hours < 16)||(n_hours = 16 && n_minutes < 46))



%% write  values to .csv file
csvwrite('GHI_20091004_offsets.csv', GHI_offset_matrix);





