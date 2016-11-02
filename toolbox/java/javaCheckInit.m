function javaCheckInit(  )
%JAVACHECKINIT Verifies that all needed java classes are loaded
%	These should generally all be done together to avoid issues with matlab
%	clearing java midway through your calculation because you loaded a new
%	class, so the best practice is just to always call this function before
%	you get started.

try
	bu.util.config.Configuration;
	com.mysql.jdbc.Driver;
	com.osisoft.jdbc.Driver;
catch err
	jpath = javaclasspath('-all');
	jfiles = {'libbu.jar'};
	for i = 1:length(jfiles)
		if( isempty( cell2mat(strfind(jpath,jfiles{i})) ) )
			javaaddpath(which(jfiles{i}));
		end
	end
end

end

