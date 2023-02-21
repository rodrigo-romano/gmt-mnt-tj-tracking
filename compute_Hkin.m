%% Function to compute Hkin based on a particular model 
%%
function [Hkin_m12, Hkin_pmtnodes, Hkin_hp, Hkin_fsm, Hkin_m2p] = ...
    compute_Hkin(ModelFolder)

fprintf('Wait... Computing Hkin from model\n%s\n',ModelFolder);
% Load model
load(fullfile(ModelFolder,"modal_state_space_model_2ndOrder.mat"),...
    'outputTable','modalDisp2Outputs');
Phi = modalDisp2Outputs; clear modalDisp2Outputs;

% M1/M2 RBM indices
m1_out = outputTable{'OSS_M1_lcl','indices'}{1};
m2_out = outputTable{'MC_M2_lcl_6D','indices'}{1};
pmt_out = outputTable{'PMT_3D','indices'}{1};
% MNT enc output indices
mntAZ_out = outputTable{'OSS_AzEncoder_Angle',"indices"}{1};
mntEL_out = outputTable{'OSS_ElEncoder_Angle',"indices"}{1};
mntGIR_out = outputTable{'OSS_RotEncoder_Angle',"indices"}{1};
% HP node displacement indices
hp_out = outputTable{'OSS_Hardpoint_D','indices'}{1};
% FSM node displacement indices
fsm_out = outputTable{'MC_M2_PZT_D','indices'}{1};
% M2 Positioner node displacement indices
m2p_out = outputTable{'MC_M2_SmHex_D','indices'}{1};


Phi_mnt_pinv = pinv([mean(Phi(mntAZ_out,1:3),1);...
    mean(Phi(mntEL_out,1:3),1);...
    mean(Phi(mntGIR_out,1:3),1)]);

% Kinematic comprensation matrices
Hkin_m12 = [Phi(m1_out,1:3); Phi(m2_out,1:3)] * Phi_mnt_pinv;
Hkin_pmtnodes = Phi(pmt_out,1:3) * Phi_mnt_pinv;
Hkin_hp = Phi(hp_out,1:3) * Phi_mnt_pinv;
Hkin_fsm = Phi(fsm_out,1:3) * Phi_mnt_pinv;
Hkin_m2p = Phi(m2p_out,1:3) * Phi_mnt_pinv;
end