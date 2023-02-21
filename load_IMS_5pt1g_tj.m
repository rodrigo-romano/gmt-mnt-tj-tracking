%
% load_IMS_5pt1g.m
% Load data required to run integrated model simulator (IMS).
%
%
% Versions
%
% - 5.1g: Some high-fidelity M1 FEM updates (Quad actuators). Mount
% controller sampling frequency matching.
% - 5p1f: high-fidelity M1 design; 2021CFD wind load inputs; DC-gain
% mismatch compensation
% - 5p1e: 2021 mount PDR design; M2 POS loop included. 
% - 5p1d: 1KHz.
% - 5p1c: AcO including 27 BMs.
% - 5p1b: AcO just with M1 and M2 rigid-body motions. Use with files of
% controllers & interfaces of version 5p1a.
% - 5p1a: Review fast-TT and FSM control loops according to P.Thompson's
% FSMS control study tech note.
%
% To do
% - Kinematic compensation in case mount SP is different from zero
% (tracking case).

% clear
simulink_fname = 'ims_Build5pt1g_tj';
fprintf("\n<-- Loading End-to-end simulation variables ->\n|-\tVersion:%s\t-|\n",simulink_fname);

%% General settings
%%
% Set telescope structural dynamics damping
sysDamp = 0.02;

% Set structural dynamics model sampling period
FEM_Ts = 1/1e3;


% Zenith angle string (affects the control model for now)
sZa = "00";

% - - - - - Simulation setting switches (0:disables) - - - - -
dc_mismatch_comp = 1; % [bool] DC mismatch compensation flag

mntFBc_en = 1;	% [bool] Mount feedback control loop switch
en_servo = 1;   % [bool] Enable/disable mount trajectory
az_ff_en = 1;	% [bool] Azimuth feedforward action switch      
el_ff_en = 1;	% [bool] Elevation feedforward action switch
gir_ff_en = 1;	% [bool] GIR feedforward action switch

m1olf_en = 0;   % [bool] M1 outer force loop switch

m2PZT_en = 0;   % [bool] M2 PZT control loop switch
m2_pos_en = 0;	% [bool] M2 Positioner control loop


%% Load telescope structural dynamics model
%%

% Model folder
switch sZa
    case '00'
        ModelFolder = fullfile(im.lfFolder,"20220611_1945_MT_mount_zen_00_m1HFN_FSM_");
        FileName = fullfile(ModelFolder,"modal_state_space_model_2ndOrder.py.mat");
    otherwise        
        ModelFolder = fullfile(im.lfFolder,"20220610_1023_MT_mount_zen_30_m1HFN_FSM_");
        FileName = fullfile(ModelFolder,"modal_state_space_model_2ndOrder.py.mat");
end

if(~exist('FEM_IO','var') || 1)
    load(FileName,'FEM_IO','Phi','Phim','eigenfrequencies','proportionalDampingVec');
    % eigenfrequencies in Hz
    fprintf('Loadinf SS model from \n%s\n', FileName);
end

%% Pick model IOs according to inputTable and outputTables
%%
% INPUTS
desiredInputLabels = FEM_IO.inputs_name(1:end);
% Remove OSS_DTA_Wind_6F input to keep compatibility
idx = strcmp(desiredInputLabels,'OSS_DTA_Wind_6F');
if (sum(idx)), desiredInputLabels(idx) = []; end
isDesired = zeros(size(Phim,1),1);
modelMuxDims = zeros(1,numel(desiredInputLabels));

cum_in_idxs = cumsum(FEM_IO.inputs_size);
FEM_In_ind_list = [cum_in_idxs-FEM_IO.inputs_size+1, cum_in_idxs];

for i1 = 1:numel(desiredInputLabels)
    aux = (FEM_In_ind_list(i1,1):FEM_In_ind_list(i1,2))';
    isDesired(aux) = 1;
    modelMuxDims(i1) = length(aux);
end
indDesInputs = find(isDesired ~= 0);
modelMuxDims(modelMuxDims == 0) = [];

% OUTPUTS
desiredOutputLabels = FEM_IO.outputs_name(1:end-9);
% Remove PMT_3D output to keep compatibility
% idx = strcmp(desiredOutputLabels,'PMT_3D');
% if (sum(idx)), desiredOutputLabels(idx) = []; end

isDesired = zeros(size(Phi,1),1);
modelDemuxDims = zeros(1,numel(desiredOutputLabels));

cum_out_idxs = cumsum(FEM_IO.outputs_size);
FEM_Out_ind_list = [cum_out_idxs-FEM_IO.outputs_size+1, cum_out_idxs];

for i1 = 1:numel(desiredOutputLabels)
    aux = (FEM_Out_ind_list(i1,1):FEM_Out_ind_list(i1,2))';
    isDesired(aux) = 1;
    modelDemuxDims(i1) = length(aux);
end
indDesOutputs = find(isDesired ~= 0);
modelDemuxDims(modelDemuxDims == 0) = [];



%% Load parameters of controllers and other subsystems
%%
% File with controller and interface parameters
ctrl_filename = sprintf('controls_5pt1g1K_z%s_llTT_oad.mat',sZa);
load(ctrl_filename,'m1sys','m2pos','m2pzt','tt7','mount','fem');

fprintf('\nLoading the IMS controller and interface file:\n%s\n',ctrl_filename);

% Extra settings
% Overwrite mount drive static friction flags
mount.az.Fr.iIsSig01Active = 0;
mount.el.Fr.iIsSig01Active = 0;
mount.gir.Fr.iIsSig01Active = 1;
girFr_Ts = fem.Ts;



%% Structural model discretization
%%
% Check whether the loaded controllers were discretized according to the
% assumed simulation frame rate
ERRMSG = sprintf('Actuator models assume a %gHz sampling frequency.\n',1/fem.Ts);
assert(fem.Ts == FEM_Ts,ERRMSG)

% Handle modal parameters
eigenfrequencies = reshape(eigenfrequencies,length(eigenfrequencies),1);
om2 = (2*pi*eigenfrequencies).^2;
twice_zom = 2*proportionalDampingVec(:).*(2*pi*eigenfrequencies(:));
adjDamp = sysDamp/(twice_zom(4)/2/sqrt(om2(4)));

% Adjust system damping
twice_zom = adjDamp*twice_zom;
om2(1:3) = 0;
twice_zom(1:3) = 0;
fprintf('Plant model damping ratio set to:%.2g\n',(twice_zom(4)/2/sqrt(om2(4))))
% Perform discretization and provide the 2nd order form DT simulation parameters
PG = zeros(length(om2),6);
for i = 1:length(om2)
    PhiGamma = expm([0 1 0; -om2(i) -twice_zom(i) 1; 0 0 0]*FEM_Ts);
    PG(i,:) = [PhiGamma(1,1:3) PhiGamma(2,1:3)];
end

phiB = Phim(indDesInputs,:)';
phiC = Phi(indDesOutputs,:);

fprintf('Telescope structural model sampling rate set to %gHz\n',1/FEM_Ts);

%% Open Simulink model
%%
open_system(simulink_fname);



%% Static gain mismatch compensation
%%

StaticModelFolder = ModelFolder;
staticFileName = fullfile(StaticModelFolder,"static_reduction_model.mat");

memory_label = simulink_fname+"/Telescope model/Structural Dynamics GMT/Psi_ss_memory";
matmult_label = simulink_fname+"/Telescope model/Structural Dynamics GMT/Psi_ss";
zoh_label = simulink_fname+"/Telescope model/Structural Dynamics GMT/Psi_ssZOH";
rt_label = simulink_fname+"/Telescope model/Structural Dynamics GMT/Psi_ssRT";

if dc_mismatch_comp
    try
        fprintf('Trying to load gain matrix from the FEM static solution from \n%s\n',...
            StaticModelFolder);
        try
            load(staticFileName,'gainMatrixMountControlled');
            gainMatrix = gainMatrixMountControlled;
        catch
            load(staticFileName,'gainMatrix');
        end

%         K_ss = phiC(:,1:end)* diag(1./((2*pi*eigenfrequencies(1:end)).^2))* phiB(1:end,:);
        K_ss = phiC(:,4:end)* diag(1./((2*pi*eigenfrequencies(4:end)).^2))* phiB(4:end,:);
        
        Psi_ss = gainMatrix(indDesOutputs,indDesInputs) - K_ss;
        for i1=1:3
            for i2=1:3
                Psi_ss(:,FEM_In_ind_list(i1,1):FEM_In_ind_list(i1,2)) = 0;
%                 Psi_ss(FEM_Out_ind_list(i2,1):FEM_Out_ind_list(i2,2),...
%                     FEM_In_ind_list(i1,1):FEM_In_ind_list(i1,2)) = 0;
            end
        end
        
        Psi_ssTs = 1/1000;        
        set_param(memory_label,'Commented','off');
        set_param(matmult_label,'Commented','off');
        set_param(zoh_label,'Commented','off');
        set_param(rt_label,'Commented','off');
    catch
        warning('Unable to compute static compensation matrix.');
        warning('Disabling static compensation matrix.');
        set_param(memory_label,'Commented','on');
        set_param(matmult_label,'Commented','on');
        set_param(zoh_label,'Commented','on');
        set_param(rt_label,'Commented','on');
        StaticModelFolder = [];
    end
else
    set_param(memory_label,'Commented','on');
    set_param(matmult_label,'Commented','on');
    set_param(zoh_label,'Commented','on');
    set_param(rt_label,'Commented','on');
    StaticModelFolder = [];
end

%% Mount trajectory profile
%%
if en_servo
    deg2rad = pi/180;   %[rad/deg] Unit convertion constant
    % Load ODC trajectories
    traj_ID = '101';
    sDir = '/Users/rromano/Workspace/mnt-odc/2022-12-20_ODC E2E Files';
    fprintf('Loading tj%s from folder\n%s\n',traj_ID,sDir);
    [mTj,~,odc_az0,odc_el0,odc_gir0] = ...
        fun_loadTjData(sDir,traj_ID,FEM_Ts);
    switch 0
        case 0
            % ODC TF filter
            f_Hz = 0.1;
            Htjf = tf(1,[1/(2*pi*f_Hz) 1]);
        case 1
            %  4th-order Bessel filter with 0.5Hz corner frequency
            Htjf = tf(97.4091, [1, 9.8141, 43.3429, 99.2538, 97.4091]);
            warning('Alternative trajectory filter selected.\n');
    end
    H_tj = c2d(Htjf,FEM_Ts, 'foh');
    
    vn = 1:length(mTj(:,1)); %1:4001;
    mount_tj.time = mTj(vn,1);
    % Compute relative trajectory
    mount_tj.signals.values = en_servo * deg2rad *...
        [mTj(vn,2), mTj(vn,3), mTj(vn,4)];
    mount_tj.signals.dimensions = 3;
else
    mount_tj.time = 0;
    mount_tj.signals.values = en_servo*[0e-6, 0e-6, 0e-6];
%     mount_tj.signals.values = [1*pi/180/3600, 0e-6, 0e-6];
    mount_tj.signals.dimensions = 3;
end

% Kinematic mount motion compensation
[Hkin_m12, Hkin_pmtnodes, Hkin_hp, Hkin_fsm, Hkin_m2p] = compute_Hkin(ModelFolder);

%% Load PMTs
%%
% Load TT PMT (PMT1)
fprintf("\nLoading performance matrix transformations...\n")
pmt1_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951 Rev B',...
    'GMT-DTA-190951_RevB_pmt1.csv');
% pmt1_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951',...
%     'GMT-DTA-190951_RevB_pmt1.csv');
pmt1 = dlmread(pmt1_fname,',',[14,3,27,302]);
if 1, fprintf("Size of PMT1:%ix%i\n",size(pmt1)); end

% Load differential piston PMT (PMT2)
pmt2_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951 Rev B',...
    'GMT-DTA-190951_RevB_pmt2.csv');
% pmt2_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951',...
%     'GMT-DTA-190951_RevB_pmt2.csv');
pmt2 = dlmread(pmt2_fname,',',[14,3,10,302]);
if 1, fprintf("Size of PMT2:%ix%i\n",size(pmt2)); end


%% Open Simulink model
%%
% Set decimation rate for Scopes & WS Logs - multiplicative factor denotes
% the sampling frequency
m1rbm_decim = 1/(fem.Ts*1000);
m2rbm_decim = 1/(fem.Ts*1000);
mount_decim = 1/(fem.Ts*1000);
tt_decim = 1/(fem.Ts*1000);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%%
return









%% Post-processing analysis
%% ...

if 0 %#ok<*UNRCH>
    figure(1);
    subplot(2,1,1);
    plot(M1rbm(:,1:3));
    ylabel('M1-S1 Txyz');
    xlabel('Samples');
    xlim([0,size(M1rbm,1)])
    grid on;
    subplot(2,1,2);
    plot(M2rbm(:,1:3));
    ylabel('M2-S1 Txyz');
    xlabel('Samples');
    xlim([0,size(M1rbm,1)])
    grid on;
end
    

%% Auxiliary functions
%%
%%%

%% ODC Functions - Taken from fun_runOne.m
%%
%% ********************************************
function [mTj,mTjInit,az0,el0,gir0] = fun_loadTjData(sDir,sidTj,t1,t2)
  load([sDir '/sim/input/tj_tracking/tg_tj' sidTj '.mat'],'W','Wtg');
  if(nargin<4), t2 = Wtg(end,1); end
  vn=fun_indexByLimits_v1(Wtg(:,1),t1,t2);
  mTj = Wtg(vn,1:7); %without acceleration <-- !!!
  mTjInit = W(vn,1:7);
  [mTj,az0,el0,gir0] = fun_remove_t0p0(mTj);
  [mTjInit,~,~,~] = fun_remove_t0p0(mTjInit);
end

%% ********************************************
function [mTj,az0,el0,gir0] = fun_remove_t0p0(mTj)
  mTj(:,1)=mTj(:,1)-mTj(1,1);
  az0=mTj(1,2); el0=mTj(1,3); gir0=mTj(1,4);
  for k=2:4
    mTj(:,k)=mTj(:,k)-mTj(1,k);
  end  
end

%% File: /base/util/fun_indexByLimits_v1.m
function vn_=fun_indexByLimits_v1(t_,t1,t2)
%Finds indexes in array t1_ between t1 and t2
%t1 and t2 can be -1, in which case they are ignored
n=length(t_);
if(t1<0)
  i1=1; 
else
  i1=find(t_>t1,1,'first');  
  if(isempty(i1)), i1=1; end
end
if(t2<0)
  tend = t_(end); tend = tend-abs(t2);
  i2=find(t_>=tend,1,'first');  
  if(isempty(i2)), i2=n; end
else
  i2=find(t_>t2,1,'first');  
  if(isempty(i2)), i2=n; end
end

vn_=i1:i2;
end