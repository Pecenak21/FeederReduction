function[y] = RemoveFileExtension(str)

%RemoveFileExtension(str)
%
% PURPOSE : Removes the file extension from a file name
%
%
% INPUT : a string (e.g. 'FileName.ext')

n=length(str);

for i=1:n
    if strcmp(str(n-i+1),'.')
        k=n-i;
        break
    end
    k=n;
end

y=str(1:k);
    