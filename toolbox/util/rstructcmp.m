function rstructcmp(cc,tt)
% compare two structures recursively using structcmp

% Brief check using structcmp function
checkf = [];
f = structcmp(cc,tt);
for i = f; i = i{1};
	fprintf('in field %s\n',i);
	if(isnumeric(cc.(i)))
		m = cc.(i)~=tt.(i);
		fprintf('numeric field %s differs in %d places by at most %g\n',i,sum(m(:)),max((cc.(i)(m)-tt.(i)(m))./cc.(i)(m)));
	elseif(iscell(cc.(i)))
		if(isnumeric(cc.(i){1}))
			ma = cell2mat(cc.(i));
			mb = cell2mat(tt.(i));
			try
				m = ma~=mb;
				fprintf('numeric field %s differs in %d places by at most %g\n',i,sum(m(:)),max((ma(m)-mb(m))./ma(m)));
			catch err
				fprintf('numeric field %s differs in %d places in matrix dimension\n',i,sum(m(:)));
			end
		elseif(ischar(cc.(i){1}))
			try
				m = ~strcmp(cc.(i),tt.(i));
				fprintf('string field %s differs in %d places\n',i,sum(m(:)));
			catch
				fprintf('some fields are strings and some are not\n');
			end
		else
			fprintf('non numeric field %s\n',i);
		end
	elseif(isstruct(cc.(i)))
		fprintf('recursively comparing field %s\n',i);
		if(length(cc.(i))>1)
			for j = 1:length(cc.(i))
				rstructcmp(cc.(i)(j),tt.(i)(j));
			end
		else
			rstructcmp(cc.(i),tt.(i));
		end
	else
		fprintf('non numeric field %s\n',i);
	end

end

end

function m = cell2mat(c)
	m = [];
	for i=1:numel(c)
		m = [m; c{i}(:)];
	end
end
