
strDirectory='D:\Sync Software\Other Code\SystemConversion\HydroOttawa';
d = excel2obj( [strDirectory '\custdata\HOL Short Circuit Results_sorted_distance.xlsx'] );

d=structconv(d.Short_Circuit_Detailed_1_);
d.Node_Id_=UnifyFormat(d,'Node_Id_','str');
d.Node_Id_=fnSanitize(d.Node_Id_);
d.Node_Id_=d.Node_Id_(1:20);
o=cymeread_dir([strDirectory '\input']);
Section=structconv(o.network.section.section);
[flg,index_FromNode]=ismember(Section.FromNodeID,d.Node_Id_);
index=find(index_FromNode);
Section=structreduce(Section,index);
[flg,index_ToNode]=ismember(Section.ToNodeID,d.Node_Id_);
index=find(index_ToNode);
Section=structreduce(Section,index);






