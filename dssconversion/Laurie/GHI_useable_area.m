%% Program for computing usable area of GHI data set

strDirectory='C:\Users\lmiller\Documents\exmainpush\SDGE\GHI';
cd(strDirectory); 

% Import the first file
fileToRead1 = ['GHI_20091004090030.csv'];
%fileToRead1 = ['GHI_20091004093000.csv'];
%fileToRead1 = ['GHI_20091004103000.csv'];
% fileToRead1 = ['GHI_20091004113000.csv'];
%fileToRead1 = ['GHI_20091004100000.csv'];
if (exist (fileToRead1)) 
    rawData1 = importdata(fileToRead1);
end % if file exists

% Import the last file
fileToRead1 = ['GHI_20091004164530.csv'];
%fileToRead1 = ['GHI_20091004160000.csv'];
%fileToRead1 = ['GHI_20091004143000.csv'];
% fileToRead1 = ['GHI_20091004133000.csv'];
% fileToRead1 = ['GHI_20091004140000.csv'];
if (exist (fileToRead1)) 
    rawData2 = importdata(fileToRead1);
end % if file exists

% zero out the NaNs
rawData1(isnan(rawData1))= 0;
rawData1(isnan(rawData2))= 0;

% if we want to use the full time period and stretch the data by a factor
% of six, need to find offsets in x and y that will work 
% range in y is 8201.25m/10 = 822 data points / 6 = 137 data points
% range in x is 11127.9m/10 = 1113 data points / 6 = 186 data points
range_y = 137;
range_x = 186; 
x_offset = 297;
y_offset = 297;
rawData1(y_offset,:)=0;
rawData1(y_offset+range_y,:)=0;
rawData1(:,x_offset)=0;
rawData1(:,x_offset+range_x)=0;

    
% % Create figure
fig_title = ['Usable area, 09:00:30 to 16:45:30'];
figureFileName = ['Usable_Area_090030_164530'];
%fig_title = ['Usable area, 09:30:00 to 16:00:00'];
%figureFileName = ['Usable_Area_093000_160000'];
%fig_title = ['Usable area, 10:30:00 to 14:30:00'];
%figureFileName = ['Usable_Area_103000_143000'];
% fig_title = ['Usable area, 11:30:00 to 13:30:00'];
% figureFileName = ['Usable_Area_113000_133000'];
% fig_title = ['Usable area, 010:00:00 to 16:45:30'];
% figureFileName = ['Usable_Area_090030_164530'];
% fig_title = ['Usable area, 10:00:00 to 14:00:00'];
% figureFileName = ['Usable_Area_100000_140000']
drawGHIFigures(rawData1, fig_title, figureFileName);






