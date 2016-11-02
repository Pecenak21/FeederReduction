function validateSC(c1,d1,p1)

[subdir subdir] = fileparts(fileparts(p1));

Vbase = 12470;
Mbase = 100000000;
Zbase = Vbase^2/Mbase;
for i=1:length(c1.linecode)% The SC validation is achieved when the line impedances are assumed to be PU
    c1.linecode(i).R1 = c1.linecode(i).R1/(Zbase);
    c1.linecode(i).X1 = c1.linecode(i).X1/(Zbase);
    c1.linecode(i).R0 = c1.linecode(i).R0/(Zbase);
    c1.linecode(i).X0 = c1.linecode(i).X0/(Zbase);
    c1.linecode(i).C1 = c1.linecode(i).C1/(Zbase);
    c1.linecode(i).C0 = c1.linecode(i).C0/(Zbase);
end

p1 = dsswrite(c1,subdir,1,subdir);

faultstudy(c1,d1,p1,1);

for i=1:length(c1.linecode)% The remaining validations (real and reactive powers) are achieved when the line impedances are assumed to be in Ohms  
    c1.linecode(i).R1 = c1.linecode(i).R1*(Zbase);
    c1.linecode(i).X1 = c1.linecode(i).X1*(Zbase);
    c1.linecode(i).R0 = c1.linecode(i).R0*(Zbase);
    c1.linecode(i).X0 = c1.linecode(i).X0*(Zbase);
    c1.linecode(i).C1 = c1.linecode(i).C1*(Zbase);
    c1.linecode(i).C0 = c1.linecode(i).C0*(Zbase);
end

p1 = dsswrite(c1,subdir,1,subdir);
