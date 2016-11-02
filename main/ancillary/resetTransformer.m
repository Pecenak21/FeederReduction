function transformer = resetTransformer(transformer)
% working setting that has been validated for intializing or reseting regcontrol
% DO NOT CHANGE if new config is not tested thoroughly
transformer.Conns = {'y' , 'y'};
transformer.Wdg = [];
transformer.Windings = 2;
transformer.Wdg = [];
transformer.kVAs = [20000 20000];
transformer.Phases = 3;
transformer.sub = '';
transformer.Rs = [];
transformer.XHL = [];
transformer.Taps = [];
end