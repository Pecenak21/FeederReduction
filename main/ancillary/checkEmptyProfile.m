function id = checkEmptyProfile(fc)
% check for empty profiles for forecast struct from 

id=find(~(sum(fc.profile>0)>0));
if ~isempty(id)
    disp('Empty profile index:');
    disp(id);
    warning('There are %d empty forecast profiles. Please check forecast data again!',length(id));
end
end