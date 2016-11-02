function regcontrol = resetRegcontrol(regcontrol)
% working setting that has been validated for intializing or reseting regcontrol
% DO NOT CHANGE if new config is not tested thoroughly
regcontrol.winding = 2;
% regcontrol.vreg = 120*.99;
regcontrol.PTPhase = 1;
regcontrol.EventLog = '';
regcontrol.vlimit = [];
regcontrol.ptratio = 57.75; % assume 12kv feeder
regcontrol.band = 2;
regcontrol.revNeutral = [];
regcontrol.reversible = '';
end