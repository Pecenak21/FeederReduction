function d = faultread(fn)
% read OpenDSS's fault-study CSV files

% read the data into a cell array
fid = fopen(fn);
tline = fgetl(fid);
if(~strcmp(tline,'FAULT STUDY REPORT'))
	y = textscan(fid,'%s %f %f %f','Delimiter',',');

	fclose(fid);

	% convert to useful struct format
	d.bus = y{1};
	d.I3Phase = num2cell(y{2});
	d.I1Phase = num2cell(y{3});
	d.LL = num2cell(y{4});
	d = structconv(d);
else
	%ALL-Node Fault Currents
	while(isempty(strfind(fgets(fid),'Bus'))); end
	fgets(fid);
	tline = fgetl(fid);
	i = 1;
	o = struct();
	while(ischar(tline) && ~isempty(tline))
		tline = regexp(tline,'[\s,]+','split');
		o(i).bus = regexprep(tline{1},'"','');
		o(i).I = str2double(tline(2:2:end));
		o(i).xr = str2double(tline(3:2:end));
		tline = fgetl(fid);
		i = i+1;
	end
	d.allnode = o';
	
	%ONE-Node to ground Faults
	while(isempty(strfind(fgets(fid),'Bus'))); end
	fgets(fid);
	tline = fgetl(fid);
	i = 1;
	o = struct();
	while(ischar(tline) && ~isempty(tline))
		tline = regexp(tline,'[\s,]+','split');
		o(i).bus = regexprep(tline{1},'"','');
		o(i).node = str2double(tline{2});
		o(i).I = str2double(tline{3});
		o(i).puV = str2double(tline(4:end));
		tline = fgetl(fid);
		i = i+1;
	end
	d.onenode = o';

	%Adjacent Node-Node Faults
	while(isempty(strfind(fgets(fid),'Bus'))); end
	fgets(fid);
	tline = fgetl(fid);
	i = 1;
	o = struct();
	while(ischar(tline) && ~isempty(tline))
		tline = regexp(tline,'[\s,]+','split');
		o(i).bus = regexprep(tline{1},'"','');
		o(i).nodes = str2double(tline(2:3));
		o(i).I = str2double(tline{4});
		o(i).puV = str2double(tline(5:end));
		tline = fgetl(fid);
		i = i+1;
	end
	d.adjNN = o';

	% cleanup:
	fclose(fid);
end

end
