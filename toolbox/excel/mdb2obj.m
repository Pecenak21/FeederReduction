function [obj moddate] = mdb2obj(dbPath)
% mdb2obj loads data from a Microsoft Access Database (.mdb) file
% Pass it a path to the database file, and it returns a struct containing
% the tables of the database.  Each table is represented as a struct array
% with fields named after the columns.

%% Setup Database Connection
% Instatiate a new ActiveX connection
access = actxserver('ADODB.Connection');

% Open connection to database
set(access,'CursorLocation',3); % magic?
Open(access, ['Provider=Microsoft.Jet.OLEDB.4.0;Data Source=' dbPath]);
% Open(access, ['Provider=Microsoft.ACE.OLEDB.12.0;Data Source=' dbPath]);
% Open(access, ['Provider=MSDASQL; Driver={Microsoft Access Driver (*.mdb)}; DBQ=' dbPath]);
% Provider=MSDASQL; Driver={Microsoft Access Driver (*.mdb)}; DBQ=C:\path\filename.mdb;

%% List Tables
% list the db schema
tabs = OpenSchema(access,'adSchemaTables');
% locate the columns we want
nam_i = 0; typ_i = 0; mod_i = 0;
for i=1:tabs.Fields.Count
	switch(tabs.Fields.Item(i-1).Name)
		case 'TABLE_NAME'
			nam_i = i;
		case 'TABLE_TYPE'
			typ_i = i;
		case 'DATE_MODIFIED'
			mod_i = i;
	end
end
if(~all([nam_i,typ_i,mod_i])), error('mdb2obj:badSchema','MDB file schema doesn''t include all the required parameters'); end
% mask to only include tables of type 'TABLE' (ignore system/access tables)
tabs = tabs.GetRows();
typ_i = strcmp(tabs(typ_i,:),'TABLE');
% extract date and name info
moddate = tabs(mod_i,typ_i)';
tabs = tabs(nam_i,typ_i)';

%% Extract Data
rmdate = false(size(moddate));
for j = 1:length(tabs)
	tn = tabs{j};
	% query the db
	try
		% This code is implicitly written for importing data from Synergee,
		% which has many tables with a SectionID column, so we'll try
		% ordering by that.  The catch block means we still get the data
		% even if this column is missing.
		try
			d1 = Execute(access,['SELECT * FROM ' tn ' order by SectionId,YearNumber;']);
		catch
			d1 = Execute(access,['SELECT * FROM ' tn ' order by SectionId;']);
		end
	catch
		try
			d1 = Execute(access,['SELECT * FROM ' tn ' order by NodeId;']);
		catch
			d1 = Execute(access,['SELECT * FROM ' tn ';']);
		end
	end
	% skip empty tables
	if(d1.RecordCount == 0); rmdate(j) = true; continue; end
	% extract fieldnames
	fns = {}';
	for i = 1:d1.Fields.Count;
		fns{i} = d1.Fields.Item(i-1).Name;
	end
	% make sure they're suitable to use as struct names
	fns = sanitize(fns');
	d1 = d1.GetRows();
	d1 = cell2struct(d1,fns);
	% and save this table to the output
	obj.(sanitize(tn)) = d1;
end
moddate(rmdate) = [];

%% Clean Up
% Clean up ActiveX resources
Close(access);
delete(access);


end

function s = sanitize(s)
s = fnSanitize(s,'x');
end
