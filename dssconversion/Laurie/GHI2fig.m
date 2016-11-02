%% Program for reading in GHI files and outputing Matlab images in .eps

strDirectory='C:\Users\lmiller\Documents\exmainpush\SDGE\GHI';
cd(strDirectory); 

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
    GHI_data = ['GHI_20091004' hours minutes seconds];
    figureFileName = GHI_data;
  
    % Import the file
    if (exist (fileToRead1)) 
        rawData1 = importdata(fileToRead1);
        rawData1(isnan(rawData1))= 0;
    % % Create figure
        fig_title = ['GHI 2009 10 04   ' hours ' ' minutes ' ' seconds];
        drawGHIFigures(rawData1, fig_title, figureFileName);

    end % if file exists
  

    % clear input from workspace
    clear rawData1;
    % close figure
    close all;

end % while ((n_hours < 16)||(n_hours = 16 && n_minutes < 46))


