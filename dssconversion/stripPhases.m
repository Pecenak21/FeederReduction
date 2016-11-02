function v = stripPhases(v)
% stripPhases removes the trailing .1.2.3 (or whatever) that indicate phase
% numbers on OpenDSS busnames

v = regexprep(v,'(\.\d)+$','');

end
