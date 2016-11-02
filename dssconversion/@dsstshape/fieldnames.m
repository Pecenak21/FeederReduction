function n = fieldnames(s)
% List the properties of the object

% we override this function so that user have some idea what kind of data
% they can set.  Also matlab may use this for tab completion with '.'
% subscripting?

% just return the field names fo the data object
n = fieldnames(s(1).data);
end