function cmpplot(myvar, myfield, lim, option, factor)
% plot comparison plot for myvar(n).myfield with 'option' field
% option: 'imagesc', 'image', 'rgb','hdr', 'rbr'

if ~exist('option','var') || isempty(option)
	option = 'imagesc';
end

if ~exist('factor','var') || isempty(factor)
	factor = 1;
end

figure;
c = length(myvar);
for i = 1:c
	h(i) = subplot(1,c,i);
	switch option
		case 'imagesc'
			imagesc(myvar(i).(myfield));
		case {'image','rgb'}
			image(myvar(i).(myfield)*factor);
		case 'hdr'
			image(previewHDR(myvar(i).(myfield)));
		case 'rbr'
			imagesc(rbr(myvar(i).(myfield),38.15));
		otherwise
			disp('Check your arguments');
	end
end

if exist('lim','var') && ~isempty(lim)
    set(h(1),'clim',lim);
    set(h(2),'clim',lim);
end

% linkprop(h,{'xlim','ylim','clim'});

end