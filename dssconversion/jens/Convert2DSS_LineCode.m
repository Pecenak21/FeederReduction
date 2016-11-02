function lc = Convert2DSS_LineCode(o)
% lc = Convert2DSS_LineCode(o)
%
% PURPOSE : Converts line codes from object to OpenDSS
%
%
% INPUT :   o: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Line Code object
%           see page 106 in OpenDSS manual V7.4.3(March 2012)

lc(length(o)) = dsslinecode;
w=2*pi()*60;
% if ischar(o(i).R1)
    for i = 1:length(o)
        lc(i).Name = o(i).ID; 
        lc(i).R1 =  str2num(o(i).R1)*o(i).ConversionFactor;
        lc(i).X1 =  str2num(o(i).X1)*o(i).ConversionFactor;
        lc(i).R0 =  str2num(o(i).R0)*o(i).ConversionFactor;
        lc(i).X0 =  str2num(o(i).X0)*o(i).ConversionFactor;
    %     lc(i).C1 = str2num(o(i).B1)/w*1e9*o(i).ConversionFactor; % OpenDSS unit is nF
    %     lc(i).C0 = str2num(o(i).B0)/w*1e9*o(i).ConversionFactor;
        lc(i).C1 = str2num(o(i).B1)/w*o(i).ConversionFactor; % OpenDSS unit is nF
        lc(i).C0 = str2num(o(i).B0)/w*o(i).ConversionFactor;
    %     lc(i).Rg = 0;
    %     lc(i).Xg = 0;
    % 	lc(i).Units = 'kft';
        lc(i).Units = o(i).Unit_length;
    end
% elseif iscell(o(i).R1)
% 	for i = 1:length(o)
%         lc(i).Name = o(i).ID; 
%         lc(i).R1 =  str2num(o(i).R1)*o(i).ConversionFactor;
%         lc(i).X1 =  str2num(o(i).X1)*o(i).ConversionFactor;
%         lc(i).R0 =  str2num(o(i).R0)*o(i).ConversionFactor;
%         lc(i).X0 =  str2num(o(i).X0)*o(i).ConversionFactor;
%     %     lc(i).C1 = str2num(o(i).B1)/w*1e9*o(i).ConversionFactor; % OpenDSS unit is nF
%     %     lc(i).C0 = str2num(o(i).B0)/w*1e9*o(i).ConversionFactor;
%         lc(i).C1 = str2num(o(i).B1)/w*o(i).ConversionFactor; % OpenDSS unit is nF
%         lc(i).C0 = str2num(o(i).B0)/w*o(i).ConversionFactor;
%     %     lc(i).Rg = 0;
%     %     lc(i).Xg = 0;
%     % 	lc(i).Units = 'kft';
%         lc(i).Units = o(i).Unit_length;
%     end  
% end

end

% 	obj.defaults = struct( ...
%     'Name'      ,{''},  ...
%     'Nphases'   ,3      ,...
%     'R1'        ,{[]}     ,...   % Ohms per unit length
%     'X1'        ,{[]}     ,... 
%     'R0'        ,{[]}     ,...
%     'X0'        ,{[]}     ,...
%     'C1'        ,{[]}     ,...    % nF per unit length
%     'C0'        ,{[]}     ,...
%     'Units'     ,{'none'}	,...	% 'mi', 'km', 'kft', 'm', 'ft', 'in', 'cm'; see note below
%     'Rmatrix'   ,{[]}     ,...
%     'Xmatrix'   ,{[]}     ,...
%     'Cmatrix'   ,{[]}     ,...
%     'BaseFreq'  ,60.0     ,...
%     'Normamps'  ,{[]}     ,...
%     'Emergamps' ,{[]}     ,...
%     'Faultrate' ,{0.0005}     ,... % Number of faults per year per unit length
%     'Pctperm'   ,100     ,... 
%     'Kron'      ,'N'     ,... 
%     'Rg'        ,0.155081    ,... % ohms per 1000 ft at 60 Hz, see note below
%     'Xg'        ,0.01805      ,... % ohms per 1000 ft at 60 Hz
%     'Rho'       ,100     ,...% meter ohms.
%     'Like'      ,'',...
% 	'repair'	,[] ...
%     );

	% Units: The default assumption ('None') for units implies that the
	% impedances are specified per unit length in the same units as the
	% line length.  OpenDSS allows you to specify any of a variety of units
	% (miles, feet, kfeet, cm, m, km), which are matched based on the first
	% two letters (or one letter for 'm').  Unit conversion only works well
	% in the case where you are applying a linecode with one set of units
	% (e.g. kft; so impedances in ohms/kft) to a line with different units
	% (e.g. ft).  In this case, the length is converted to the same units
	% as impedance when calculating the Y matrix, or the impedances are
	% converted to the same units as the length (and the 'units' param of
	% the line) for display (by OpenDSS; we're not that sophisticated here
	% yet).  Attempting to specify impedance and line length with different
	% units for just one line object is tricky.  Conversions are applied as
	% before (i.e. to the length for calculation purposes or the impedance
	% for display purposes) however _setting_ the unit is complicated:
	%	* The first time a unit is set serves only to adjust the 'current'
	%	unit, and does not effect the scale factor, because opendss
	%	refuses to convert from 'none' units to anything else
	%	* setting any of the impedance values or a linecode resets the
	%	units to none
	%	* by applying the units property a SECOND time, you can achieve the
	%	same effect as specifying a linecode with the first set of units
	%	and a line with the second set.  e.g. to interpret impedances in
	%	ohms/kft and lengths in ft, one might do:
	%	New line.myline R=0.05 length=534 units=kft units=ft
	%	Unfortunately, the current matlab class architecture we've setup is
	%	incapable of specifying the same parameter twice, so for creating
	%	OpenDSS files that behave like this, you either need to use
	%	linecodes or else just convert the units yourself.
    
%     Rg: Carson earth return resistance per unit length used to compute impedance values
%           at base frequency. See description above. For making better adjustments of line
%           impedance values for frequency for harmonics studies. Default= 0.01805 ohms per
%           1000 ft at 60 Hz. If you do not wish to adjust the earth return impedance for
%           frequency, set both Rg and Xg to zero. Generally avoid Kron reduction if you will be
%           solving at frequencies other than the base frequency and wish to adjust the earth
%           return impedance.
