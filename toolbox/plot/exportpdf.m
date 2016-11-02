function exportpdf(filename,savepath)
% To export current figure to a nice pdf format
% Format: no boundary, white background
% Usage: exportpdf('myfilename');
%
% NOTE: Only works for Linux and MAC OSX at the moment. Sorry Windows!!!

if ~exist('savepath','var') || isempty(savepath)
	savepath = pwd;
end

% print(gcf,'-dpdf',filename);
style = hgexport('factorystyle');
style.Format = 'pdf';
hgexport(gcf,['./' filename '.pdf'],style,'Format','pdf')
if ismac
	system(['/usr/texbin/pdfcrop --hires --margins 1 --gscmd /usr/local/bin/gs --pdftexcmd /usr/texbin/pdftex ' filename]);
else %Ubuntu
	system(['/usr/bin/pdfcrop ' filename]);
end

system(['mv ' filename '-crop.pdf' ' ' savepath '/' filename '.pdf']);
system(['rm ' filename '.pdf']);

% saveas(gcf,[path '/' filename '.png']);
% saveas(gcf,[path '/' filename '.fig']);
end
