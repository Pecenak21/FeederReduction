function dssstructcmp(cc,tt)
% compare two opendss structs objects (created from dssconversion or
% dssparse). Will report difference in numeric values if possible in form of fractional difference.
% Usage:
%			dssstructcmp(circuit1, circuit2)

% Brief check using structcmp function
checkf = [];
f = structcmp(cc,tt);
for i = f
	i = i{1};
	fprintf('in field %s\n',i);
	try
		checkf.(i) = structcmp(cc.(i),tt.(i));
	catch
		fprintf('\tonly on one side\n');
	end
end;

% calculate fractional difference for numeric values or report number of
% differences for other types
for i = f; i=i{1};
	for j=checkf.(i); j=j{1};
		if(isnumeric(cc.(i).(j)))
			m = cc.(i).(j)~=tt.(i).(j);
			fprintf('numeric field %s.%s differs in %d places by at most %g\n',i,j,sum(m(:)),max((cc.(i).(j)(m)-tt.(i).(j)(m))./cc.(i).(j)(m)));
		elseif(iscell(cc.(i).(j)))
			if(isnumeric(cc.(i).(j){1}))
				ma = cell2mat(cc.(i).(j));
				mb = cell2mat(tt.(i).(j));
				try 
					m = ma~=mb;
					fprintf('numeric field %s.%s differs in %d places by at most %g\n',i,j,sum(m(:)),max((ma(m)-mb(m))./ma(m)));
				catch err
                    if exist('m','var')
                        fprintf('numeric field %s.%s differs in %d places in matrix dimension\n',i,j,sum(m(:)));
                    else
                        fprintf('numeric field %s.%s differs in matrix dimension\n',i,j);
                    end
				end
			elseif(ischar(cc.(i).(j){1}))
				try
					m = ~strcmp(cc.(i).(j),tt.(i).(j));
					fprintf('string field %s.%s differs in %d places\n',i,j,sum(m(:)));
				catch
					fprintf('some fields are strings and some are not\n');
				end
			else
				fprintf('non numeric field %s.%s\n',i,j);
			end
		else
			fprintf('non numeric field %s.%s\n',i,j);
		end
	end
end

end

function m = cell2mat(c)
	m = [];
	for i=1:numel(c)
		m = [m; c{i}(:)];
	end
end
