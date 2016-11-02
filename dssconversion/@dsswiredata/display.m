function disp(s)
% Display the contents of the object

% override the builtin display function to also show the object contents
% after we print its type
builtin('disp',s)
disp(get(s));
end
