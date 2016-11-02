function [ ob ] = excel2obj( filename, rowToStartFrom )
% Convert excel file with multiple sheets to object oriented struct with
% each component to be the name of each sheet.
% 
% inputs:
%		filename
%		rowToStartFrom :if there are some additional data on top that you might want to skip

if( ~strcmp(filename(1:2),'\\') && filename(2)~=':' && ~exist(filename,'file') )
	filename = [pwd '/' filename];
end

d = excelRead(filename);
%% Get sheets' names
sns = {d(:).name};

%%
% Loop through each datasheet
for i = 1:length(d)

	% if no data in the sheet, skip
	if(isempty( d(i).text )), continue; end;

	% get raw data for that sheet
	r = d(i).raw;

	% get struct of data organized by sheet's columns
	ob.(sanitize(sns{i})) = structconv(sheetRead( r ));
end

function sh = sheetRead( r )

if ~exist('rowToStartFrom','var')
	rowToStartFrom = 1;
end

% Remove blank rows
rows = all(cellfun('length',r)==1,2);
rows2 = all(cellfun(@(c)(isnan(c)),r(rows,:)),2);
rows = find(rows); rows = rows(rows2);
r(rows,:) = [];

% Get columns' names
ns = r(rowToStartFrom,:);

% If no data, return
if size(r,1) < rowToStartFrom + 1; sh = []; return; end;

% loop through each column and record the data starting from 2nd row
for j = 1:length(ns)
	try
		if isnan( ns{j} )
			continue;
		end;
		sh.(ns{j}) = r(rowToStartFrom + 1:end,j);
	catch err
		if strcmp( err.identifier, 'MATLAB:AddField:InvalidFieldName')
			% when an invalid field name is specified, fix it:
			sh.(sanitize(ns{j})) = r(rowToStartFrom + 1:end,j);
		elseif strcmp( err.identifier, 'MATLAB:mustBeFieldName')
			% when a non-char field name is specified, consider fixing it
			% if it's a number, otherwise issue a warning and skip it
			if(isnumeric(ns{j}) && ~isempty(ns{j}))
				sh.(sanitize(num2str(ns{j}))) = r(rowToStartFrom + 1:end,j);
			else
				warning('excel2obj:badColumnName','couldn''t convert column %d header to field name.',j);
			end
		else
			sh.(ns{j}) = [];
		end
	end
end

if ~exist('sh','var') || isempty(sh)
	sh = [];
end
end

end

function [ o ] = excelRead( fn )
% Read excel file and put it into a struct format with each item is a sheet
% with name, num (number data), text (text data), raw (raw data).
[fi, nms] = xlsfinfo(fn);
for i = 1:length(nms)
	o(i).name = nms{i};
	[o(i).num, o(i).text, o(i).raw] = xlsread(fn,i);
end
end

function s = sanitize(s)
s = fnSanitize(s,'x');
end
