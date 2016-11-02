function[fn] = GetFieldNames(opt,varargin)
% function[fn] = GetFieldNames(opt,varargin)
%
% PURPOSE : Returns fieldnames of input structures
%
%
% INPUT : opt:      'intersect' returns only field names that are common
%                   'unique' returns unique field names
%                   'all' returns all field names
%         varargin: struct1,fieldname1,struct2,fieldname2,...
%
%
% OUTPUT : Field Names; intersection, unique, or all

fn=[];
for i=1:2:nargin-1
    if isfield(varargin{i},varargin(i+1))
        y=fieldnames(varargin{i}.(varargin{i+1}));
        fn=unique([fn;y]);
    end
end
if strcmp(opt,'intersect')
    for i=1:2:nargin-1
        if isfield(varargin{i},varargin(i+1))
            y=fieldnames(varargin{i}.(varargin{i+1}));
            fn=intersect(fn,y);
        end
    end
elseif strcmp(opt,'intersect')
    fn=unique(fn);
end
end

