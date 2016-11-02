function[struc_out] = StandardizeFieldNames(struc,opt)
% StandardizeFieldNames(struc,opt)
%
% PURPOSE : translate field names of a structure so they follow a
% standardized format
%
%
% INPUT :   struc: as structure with field names
%           opt: option based on the source of the data
%                1: Data provided by HydroOttawa for CEATI Inrush project
%
% OUTPUT :  a structure with renamed field names

if opt==1
    FieldName_old='Feeder_Id_';
    FieldName_new='FeederID';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='Node_Id_';
    FieldName_new='NodeID';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='Equipment_Id_';
    FieldName_new='EquipmentID';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='Phase_';
    FieldName_new='Phase';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='kVLN_';
    FieldName_new='kVLN';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='C_Fact_';
    FieldName_new='CFact';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLL_Amps_';
    FieldName_new='Symmetrical_Amps_3Ph';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLL_Kmax_Amps_';
    FieldName_new='Symmetrical_Amps_3Ph_Kmax';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLL_Kmax_Z_Amps_';
    FieldName_new='Symmetrical_Amps_3Ph_Kmax_Z';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLL_Kmin_Amps_';
    FieldName_new='Symmetrical_Amps_3Ph_Kmin';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLL_Kmin_Z_Amps_';
    FieldName_new='Symmetrical_Amps_3Ph_Kmin_Z';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLG_Amps_';
    FieldName_new='Symmetrical_Amps_LLG';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLG_Kmax_Amps_';
    FieldName_new='Symmetrical_Amps_LLG_Kmax';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LLG_Kmin_Amps_';
    FieldName_new='Symmetrical_Amps_LLG_Kmin';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LL_Amps_';
    FieldName_new='Symmetrical_Amps_LL';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LL_Min_Amps_';
    FieldName_new='Symmetrical_Amps_LL_Min';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LL_Kmax_Amps_';
    FieldName_new='Symmetrical_Amps_LL_Kmax';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LL_Kmin_Amps_';
    FieldName_new='Symmetrical_Amps_LL_Kmin';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LL_Kmax_Z_Amps_';
    FieldName_new='Symmetrical_Amps_LL_Kmax_Z';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LL_Kmin_Z_Amps_';
    FieldName_new='Symmetrical_Amps_LL_Kmin_Z';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LG_Max_Amps_';
    FieldName_new='Symmetrical_Amps_LG_Max';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LG_Min_Amps_';
    FieldName_new='Symmetrical_Amps_LG_Min';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LG_Kmax_Amps_';
    FieldName_new='Symmetrical_Amps_LG_Kmax';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LG_Kmax_Z_Amps_';
    FieldName_new='Symmetrical_Amps_LG_Kmax_Z';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LG_Kmin_Amps_';
    FieldName_new='Symmetrical_Amps_LG_Kmin';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='LG_Kmin_Z_Amps_';
    FieldName_new='Symmetrical_Amps_LG_Kmin_Z';
    struc = RenameField(struc,FieldName_old,FieldName_new);
    FieldName_old='Total_distance_ft';
    FieldName_new='Distance_ft';
    struc = RenameField(struc,FieldName_old,FieldName_new);

end

struc_out=struc;