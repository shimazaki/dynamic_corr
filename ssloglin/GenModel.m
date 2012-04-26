function model = GenModel(raw,struct,init)

model.struct = struct;

%%%%%%%%%%%%%%%%%%%% Generate Binary Data %%%%%%%%%%%%%%%%%%
model.binary = GenBinary(raw,struct);
		
%%%%%%%%%%%%%%%%%%%% Initalization  %%%%%%%%%%%%%%%%%%
model.init = GenInit(init,model.binary.d);

%%%%%%%%%%%%%%%%%%%% Estimation %%%%%%%%%%%%%%%%%%
model.param = ssloglin(model.binary,model.init);

