function mailsend(recepient,subject,message,attachments)
% simple function to send notification mail from a common email user for
% our solar group (solarucsd@gmail.com)
%
% Usage example:
%			mailsend('andunguyen@ucsd.edu','Subject: Santa Claus coming to town','Message: Cheers');

% Process inputs

% Config gmail client
persistent mailclient;
if isempty(mailclient)
	mailclient.address = 'solarucsd@gmail.com';
	mailclient.password = 'solar@ucsd';
	
	setpref('Internet','E_mail',mailclient.address);
	setpref('Internet','SMTP_Server','smtp.gmail.com');
	setpref('Internet','SMTP_Username',mailclient.address);
	setpref('Internet','SMTP_Password',mailclient.password);
	
	props = java.lang.System.getProperties;
	props.setProperty('mail.smtp.auth','true');
	props.setProperty('mail.smtp.socketFactory.class', ...
		'javax.net.ssl.SSLSocketFactory');
	props.setProperty('mail.smtp.socketFactory.port','465');
end

% clean up email input
recepient = regexp(recepient,'[\w\.]*@[\.\w]*','match');

if isempty(recepient)
	warning('mailsend:Novalidrecepient','No valid emailing recepients!!!');
end

sendmail(recepient,subject,message);

end