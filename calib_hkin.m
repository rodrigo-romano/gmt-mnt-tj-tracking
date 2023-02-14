%
% calib_hkin.m
%
% Script to compute Hkin matrices from simulation data
%

load('AZ_calib_za30.mat','mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
AZ_ = mountY.signals.values(end,2);
HkinAZ_30_pmt = [PMT1_DT(end,:), PMT2_DT(end,:)]' ./ AZ_;
HkinAZ_30_m12 = [M1rbm(end,:), M2rbm(end,:)]' ./ AZ_;

load('EL_calib_za30.mat','mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
EL_ = mountY.signals.values(end,4);
HkinEL_30_pmt = [PMT1_DT(end,:), PMT2_DT(end,:)]' ./ EL_;
HkinEL_30_m12 = [M1rbm(end,:), M2rbm(end,:)]' ./ EL_;

Hkin_m12_z30 = [HkinAZ_30_m12, HkinEL_30_m12, zeros(84,1)];
Hkin_pmt_z30 = [HkinAZ_30_pmt, HkinEL_30_pmt, zeros(21,1)];

% save("Hkin_z30","Hkin_m12_z30","Hkin_pmt_z30");


%% Hkin test data
%%

load('AZ_calib_za30.mat','mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
% load('EL_calib_za30.mat','mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');

%%
% Load Hkin
ModelFolder = fullfile(im.lfFolder,...
    "20220610_1023_MT_mount_zen_30_m1HFN_FSM_");
load(fullfile(ModelFolder,"pointingToRBM_FEM.mat"),...
    'az2M1','el2M1','az2M2','el2M2');

%% Test data mount motion plots
figure(874)
ax_labels = ["AZ","EL","GIR"];
for ik = 1:3
    subplot(3,1,ik)
    plot(mountY.time, mountY.signals.values(:,ik*2),'LineWidth',1.5);
    hold on;
    plot(mountY.time, mountY.signals.values(:,1+(ik-1)*2),'--','LineWidth',1.5);
    legend("Enc average","set-point",'Orientation','horizontal');
    ylabel(sprintf("%s (rad)",ax_labels(ik)));
    hold off;
end


%% Load LOM transformations
%%
% Segment TT optical sensitivity matrix
if(~exist('D_seg_tt','var') || false)
    ttOSMfile = fullfile(im.lfFolder,'LOM-data','lom_tt_dt.mat');
    load(ttOSMfile,'D_seg_tt');
    fprintf('\nTT sensitivity matrix loaded from \n%s\n',ttOSMfile);
end
% Piston optical sensitivity matrix
if(~exist('D_seg_piston','var') || false)
    pistonOSMfile = fullfile(im.lfFolder,'LOM-data','D_seg_piston_dt.mat');
    load(pistonOSMfile,'D_seg_piston');
    fprintf('\nPiston sensitivity matrix loaded from \n%s\n',pistonOSMfile);
end


%% Assess kinematic compensation

% Apply kinematic compensation
mnt_sp = mountY.signals.values(:,1:2:end);
pmt12 = [PMT1_DT, PMT2_DT] - mnt_sp*Hkin_pmt_z30';
pmt12_ = [PMT1_DT, PMT2_DT] - mnt_sp*Hkin_pmt_';
% M1 and M2 RBM after kinematic compensation
m1m2 = [M1rbm,M2rbm] - mnt_sp*Hkin_m12_z30';
m1m2_ = [M1rbm,M2rbm] - mnt_sp*Hkin_m12_';
% m1m2_ = [M1rbm,M2rbm] - mnt_sp*[[az2M1;az2M2],[el2M1;el2M2],zeros(84,1)]';

vn = 10001:20001;
figure(875)
subplot(3,1,1)
plot(mountY.time(vn), 1e3*10.2*pmt12_(vn,2:2:14))
ylabel('Seg TT (R_x) (nm)'); grid on;
subplot(3,1,2)
plot(mountY.time(vn), 1e3*10.2*pmt12_(vn,1:2:14))
ylabel('Seg TT (R_y) (nm)'); grid on;
subplot(3,1,3)
plot(mountY.time(vn), 1e9*pmt12_(vn,15:end))
ylabel('Seg piston (nm)'); grid on;
xlabel("Time (s)");

% PTT computed from M1 and M2 rigid-body motions
segtt = m1m2_ * D_seg_tt';
segp = m1m2_ * D_seg_piston';
segdp = segp - mean(segp,2);

figure(876)
subplot(3,1,1)
plot(mountY.time(vn), 10.2*segtt(vn,1:7)*180/pi*3600*1e3);
ylabel('Seg TT (R_x) (nm)'); grid on;
subplot(3,1,2)
plot(mountY.time(vn), 10.2*segtt(vn,8:14)*180/pi*3600*1e3);
ylabel('Seg TT (R_y) (nm)'); grid on;
subplot(3,1,3)
plot(mountY.time(vn), 1e9*segdp(vn,:));
ylabel('Seg piston (nm)'); grid on;
xlabel("Time (s)");

%%
return

% %% ODC's Hkin
% load('Hkin_v20p9.mat', 'Hkin_za30_HcTp19')
% Hkin_pmt = cell2mat(tfdata(Hkin_za30_HcTp19(4:end-1,:)));

