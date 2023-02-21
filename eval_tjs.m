%
% eval_tjs.m
% Script to evaluate the effect of mount trajectories on the acquired
% image.


%% Preamble
%%
% radians to mas conversion factor
rad2mas = 1e3*(180/pi)*3600;
% Root of the mean squared value function
rms = @(x,dir)squeeze(sqrt(mean(x.^2,dir)));
% Flag to recompute Hkin
recomputeHkin = true;%false;
% Figure number index
figidx = 0;

% ODC folder of utilities
odc_file_folder = '/Users/rromano/Workspace/mnt-odc';
odc_base_folder = fullfile(odc_file_folder,...
    '/2022-06-09_ODC Dynamic Wind Simulations_Simulink & Matlab files',...
    '/base/util');
addpath(odc_base_folder);


%% Load simulation data
%%
clear mountY
% sim_label = 'za00';
sim_label = 'za30';
% sim_label = 'm1ofl_m2fsm';  % M1OFL:ON, M2FSM/POS:ON
switch sim_label
    case 'm1ofl_m2fsm'
        FEM_LABEL = "20220611_1945_MT_mount_zen_00_m1HFN_FSM_";
        load("tj101_za00_m1ofl_m2fsm.mat",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj102_za00_m1ofl_m2fsm.mat",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
    case 'za00'
        FEM_LABEL = "20220611_1945_MT_mount_zen_00_m1HFN_FSM_";
%         load("tj101_za00_mount_withff",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
        load("tj102_za00_mount_noff",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
    case 'za30'
        FEM_LABEL = "20220610_1023_MT_mount_zen_30_m1HFN_FSM_";
        load("tj103_za30_mount_noff.mat",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj101_za30_mount_noff.mat",'mount','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj101_za30_mount_withff.mat",'mount','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj101_za30_limnt_noff.mat",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj101_za30_4thbessel_noff",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
end

try mount = mountY;
catch, fprintf('Variable mountY is not available!\n');
end


%% Load PMTs

if(~exist('pmt1','var'))
    % Load TT PMT (PMT1)
    fprintf("\nLoading performance matrix transformations:\n")
    pmt1_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951 Rev B',...
        'GMT-DTA-190951_RevB_pmt1.csv');
    % pmt1_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951',...
    %     'GMT-DTA-190951_RevB_pmt1.csv');
    pmt1 = dlmread(pmt1_fname,',',[14,3,27,302]); %#ok<*DLMRD> 
    if 1, fprintf("Size of PMT1:%ix%i\n",size(pmt1)); end
end
if(~exist('pmt2','var'))
    % Load differential piston PMT (PMT2)
    pmt2_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951 Rev B',...
        'GMT-DTA-190951_RevB_pmt2.csv');
    % pmt2_fname = fullfile(im.lfFolder,'PMTs','GMT-DTA-190951',...
    %     'GMT-DTA-190951_RevB_pmt2.csv');
    pmt2 = dlmread(pmt2_fname,',',[14,3,10,302]);
    if 1, fprintf("Size of PMT2:%ix%i\n",size(pmt2)); end
end

%% Compute Hkin
if( ~exist('Hkin','var')|| ~exist('Hkin_pmt','var')|| recomputeHkin || 0)
    ModelFolder = fullfile(im.lfFolder, FEM_LABEL);
    [Hkin, Hkin_pmtnodes] = compute_Hkin(ModelFolder);
    Hkin_pmt = [pmt1;pmt2]*Hkin_pmtnodes;
end

t = mount.time;     % Time vector [s]
Ts = t(2)-t(1);



%% Compare Hkin with the one from ODC
if(false)
    load('Hkin_v20p9.mat', 'Hkin_za00_HcTp19');
    odcHkin = cell2mat(tfdata(Hkin_za00_HcTp19(4:end-1,:)));
    figure(figidx+800)
    subplot(221), plotyy(1:14,-Hkin_pmt(1:14,1),1:14,odcHkin(1:14,1)); %#ok<*PLOTYY> 
    subplot(222), plotyy(1:7,Hkin_pmt(15:21,1),1:7,odcHkin(1:7,1));
    subplot(223), plotyy(1:14,-Hkin_pmt(1:14,2),1:14,odcHkin(1:14,2));
    subplot(224), plotyy(1:7,Hkin_pmt(15:21,2),1:7,odcHkin(1:7,2));
end


%% Kinematic compensation
%%

% M1 and M2 RBM after kinematic compensation
m1m2RBM = [M1rbm,M2rbm] - mount.signals.values(:,1:2:6) * Hkin';
% PMT1/2 outputs after kinematic compensation
pmt1pmt2 = [PMT1_DT,PMT2_DT] - mount.signals.values(:,1:2:6) * Hkin_pmt';


%% Apply LOM transformations
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

segtt = m1m2RBM * D_seg_tt';
segp = m1m2RBM * D_seg_piston';
segdp = segp - mean(segp,2);


%% Apply Rejection Transfer Functions (RTFs)
%%
plot_rtf = false;
% rtf = oad_AOTT_rtf();
hLtao = fun_rtf('LTAO', false, false);
% rtfP = getRTF_DP(0);
p_rtf_label = 'OIWFS125'; %'AGWS'; %
switch p_rtf_label
    case 'AGWS', rtfP = getRTF_DP(0);
    case 'OIWFS10', rtfP = getRTF_DP(1,1/10,50e-3);
    case 'OIWFS50', rtfP = getRTF_DP(1,1/50,10e-3);
    case 'OIWFS125', rtfP = getRTF_DP(1,1/125,10e-3);
end
fprintf('LTAO piston RTF variant: %s\n',p_rtf_label);

if true
    hLtaoP = fun_rtf('LTAOPiston', false, false);
    hbode = bodeoptions;
    hbode.FreqUnits = 'Hz'; hbode.grid = 'on';
    figure(figidx+110); w = 2*pi*logspace(-3,3,2001);
%     subplot(1,2,1); 
    bodemag(rtfP,getRTF_DP(1,1/50,20e-3),w,hbode);
%     subplot(1,2,2); impulse(pade(rtfP,10),pade(getRTF_DP(1,1/50,20e-3),10));
    legend('T_p = 1/50s, \tau_p = 10e-3','T_p = 1/50s, \tau_p = 20e-3'); grid on;
    set(gca,'FontSize',12);
end

segttf = fun_applyRtfAndCut(t, segtt ,hLtao, 1:length(t));
% segdpf = fun_applyRtfAndCut(t, segdp ,hLtaoP, 1:length(t));
segdpf = fun_applyRtfAndCut(t, segdp ,rtfP, 1:length(t));

%% PMT

% pmt1f = zeros(size(segtt,1),14);
pmt1f = fun_applyRtfAndCut(t, pmt1pmt2(:,1:14) ,hLtao, 1:length(t));
% pmt2f = fun_applyRtfAndCut(t, pmt1pmt2(:,15:21) ,hLtaoP, 1:length(t));
pmt2f = fun_applyRtfAndCut(t, pmt1pmt2(:,15:21) ,rtfP, 1:length(t));


%% Time-domain plots
%%

% PLot initial time instant
t1 = 16;

% GMTO LOM output plot
if 0
    figure(figidx+15) %#ok<*UNRCH>
    set(gcf,'position',[123   230   640   400])
%     segttp = [rad2mas*segttf,1e9*segdpf];
    segttp = [rad2mas*segtt,1e9*segdp];
    tt_labels= ["Seg X-motion (mas)","Seg Y-motion (mas)","Seg Piston (nm)"];
    for ik = 1:3
        subplot(3,1,ik)
        plot(t,segttp(:,(1:7)+(ik-1)*7));
        ylabel(tt_labels{ik});
        grid on; xlim([t1, t(end)]);
        if(ik == 1), title('GMTO-LOM TT and Piston'); end
    end
    
    xlabel("Time (s)");
    legend("S1","S2","S3","S4","S5","S6","S7",...
        'Location','southwest','Orientation','horizontal','NumColumns',7);
    legend boxoff
end

% PMT output plot
if 0
    figure(figidx+16)
    set(gcf,'position',[123   230   640   400])
%     segttp = [1e3*pmt1f,1e9*pmt2f];
    segttp = [1e3*pmt1pmt2(:,1:14),1e9*pmt1pmt2(:,15:end)];
    tt_labels= ["Seg X-motion (mas)","Seg Y-motion (mas)","Seg Piston (nm)"];
    for ik = 1:3
        subplot(3,1,ik)
        plot(t,segttp(:,(1:7)+(ik-1)*7));
        ylabel(tt_labels{ik});
        grid on; xlim([t1, t(end)]);
        if(ik == 1), title('PMT1 (TT) and PMT2 (Piston)'); end
    end
    
    xlabel("Time (s)");
    legend("S1","S2","S3","S4","S5","S6","S7",...
        'Location','southwest','Orientation','horizontal','NumColumns',7);
    legend boxoff
end

%% PSD plots
%%
if 0
    % Length of data segment
    M = 128;
    nFFT = 4096;
    detrend = 'none'; %'mean'; %

    PpsdLOM = zeros(nFFT/2+1,7); PpsdLOMf = zeros(nFFT/2+1,7);
    PpsdPMT = zeros(nFFT/2+1,7); PpsdPMTf = zeros(nFFT/2+1,7);
    for ik = 1:7
        [PpsdLOM(:,ik),freqP] = utils.pwelch(1e9*segdp(:,ik)...
            ,2*M,[],nFFT,1/Ts,'onesided',detrend);
        [PpsdPMT(:,ik),~] = utils.pwelch(1e9*pmt1pmt2(:,ik+14)...
            ,2*M,[],nFFT,1/Ts,'onesided',detrend);
        [PpsdLOMf(:,ik),~] = utils.pwelch(1e9*segdpf(:,ik)...
            ,2*M,[],nFFT,1/Ts,'onesided',detrend);
        [PpsdPMTf(:,ik),~] = utils.pwelch(1e9*pmt2f(:,ik)...
            ,2*M,[],nFFT,1/Ts,'onesided',detrend);
    end
    segdp_rms = rms(segdp,2);
    [PpsdLOM_rms,~] = utils.pwelch(1e9*segdp_rms...
        ,2*M,[],nFFT,1/Ts,'onesided',detrend);
%     segdpf_rms = rms(segdpf,2);

    figure(figidx+17)
    set(gcf,'position',[123   230   1.3*640   400])
    subplot(121)
    loglog(freqP,PpsdLOM,'-');
    title("Segment Piston PSD (GMT-LOM)")
    grid on; hold on; set(gca,'FontSize',12);
%     loglog(freqP,PpsdLOMf,'--');
    loglog(freqP,PpsdLOM_rms,'k--');
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('Seg Piston (nm^2/Hz)')
    xlim([min(freqP) max(freqP)]); axis tight;
%     legend("S1","S2","S3","S4","S5","S6","S7",...
%         "S1","S2","S3","S4","S5","S6","S7",...
%         'Location','southwest','NumColumns',2);
    legend boxoff

    subplot(122)
    loglog(freqP,PpsdPMT,'-');
    title("Segment Piston PSD (PMT2)")
    grid on; hold on; set(gca,'FontSize',12);
    loglog(freqP,PpsdPMTf,'--');
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('Seg Piston (nm^2/Hz)')
    xlim([min(freqP) max(freqP)]); axis tight;
    legend("S1","S2","S3","S4","S5","S6","S7",...
        "S1","S2","S3","S4","S5","S6","S7",...
        'Location','southwest','NumColumns',2);
    legend boxoff
end

if(0)
    size(rms(1e9*segdp,2))
    % Length of data segment
    M = 128;
    nFFT = 4096;
    detrend = 'none'; %'mean'; %

    PpsdLOM = zeros(nFFT/2+1,7); PpsdLOMf = zeros(nFFT/2+1,7);
    PpsdPMT = zeros(nFFT/2+1,7); PpsdPMTf = zeros(nFFT/2+1,7);
    for ik = 1:7
        [PpsdLOM(:,ik),freqP] = utils.pwelch(1e9*segdp(:,ik)...
            ,2*M,[],nFFT,1/Ts,'onesided',detrend);
        [PpsdPMT(:,ik),~] = utils.pwelch(1e9*pmt1pmt2(:,ik+14)...
            ,2*M,[],nFFT,1/Ts,'onesided',detrend);
    end

    figure(figidx+17)
    set(gcf,'position',[123   230   1.3*640   400])
    subplot(121)
    loglog(freqP,PpsdLOM,'-');
    title("Segment Piston PSD (GMT-LOM)")
    grid on; hold on; set(gca,'FontSize',12);
    loglog(freqP,PpsdLOMf,'--');
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('Seg Piston (nm^2/Hz)')
    xlim([min(freqP) max(freqP)]); axis tight;
%     legend("S1","S2","S3","S4","S5","S6","S7",...
%         "S1","S2","S3","S4","S5","S6","S7",...
%         'Location','southwest','NumColumns',2);
    legend boxoff

    subplot(122)
    loglog(freqP,PpsdPMT,'-');
    title("Segment Piston PSD (PMT2)")
    grid on; hold on; set(gca,'FontSize',12);
    loglog(freqP,PpsdPMTf,'--');
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('Seg Piston (nm^2/Hz)')
    xlim([min(freqP) max(freqP)]); axis tight;
    legend("S1","S2","S3","S4","S5","S6","S7",...
        "S1","S2","S3","S4","S5","S6","S7",...
        'Location','southwest','NumColumns',2);
    legend boxoff
end


%% RTF comparison plots
%%

if true
    AGWSrtf = getRTF_DP(0);
    OIWFS10rtf = getRTF_DP(1,1/10,50e-3);
    OIWFS50rtf = getRTF_DP(1,1/50,10e-3);
    OIWFS125rtf = getRTF_DP(1,1/125,10e-3);
    hbode = bodeoptions;
    hbode.FreqUnits = 'Hz'; hbode.grid = 'on';
    hbode.Title.String = ''; hbode.XLabel.FontSize = 12;
    hbode.YLabel.FontSize = 12; hbode.TickLabel.FontSize = 12;    
    figure(111);
    bodemag(AGWSrtf,OIWFS10rtf,OIWFS50rtf,OIWFS125rtf,...
        2*pi*logspace(-3,3,501),hbode);
    legend('AGWS RTF','OIWFS-10 RTF','OIWFS-50 RTF','OIWFS-125 RTF',...
        'Location','southeast'); legend boxoff;
end

%% WFE
kt1 = 25*(1/Ts);
[wfe,TTwfe,Pwfe] = calc_WFE(rad2mas*segttf(kt1:end,:), 1e9*segdpf(kt1:end,:));
fprintf('GMTO-LOM: %.2gnm\t%.3gnm\t%.3gnm\t\n',TTwfe,Pwfe,wfe);
[wfe_,TTwfe_,Pwfe_] = calc_WFE(1e3*pmt1f(kt1:end,:), 1e9*pmt2f(kt1:end,:));
fprintf('PMT-WFE: %.2gnm\t%.3gnm\t%.3gnm\t\n',TTwfe_,Pwfe_,wfe_);

% [wfe_,TTwfe_,Pwfe_] = calc_WFE(1e3*fun_detrend_v1(pmt1f(kt1:end,:),1),...
%     1e9*fun_detrend_v1(pmt2f(kt1:end,:),1));
% fprintf('PMT-WFE (without linear term): %.3gnm\t%.3gnm\t%.3gnm\t\n',TTwfe_,Pwfe_,wfe_);


%%
rmpath(odc_base_folder);
return










%% Tests

load(fullfile('/Users/rromano/Workspace/mnt-odc/2022-12-20_ODC E2E Files',...
    'res_eval/tracking/input/v20.9/tj',...
    'HbTp19/za30_v20.9_HbTp19_wlc0_tj101_nl1.mat'),'r');

%% Mount TO/PO ODC-GMTO comparison

load("tj101_za30_mount_noff",'mount','mountTo');
mount_ = mount; mountTo_ = mountTo;
load("tj101_za30_4thbessel_noff",'mountY','mountTo');
mount = mountY; mountTo = mountTo(:,(1:3)+3*1); % 1:3 Includes friction and parasitic torques
                                                % 4:6 Controller output +
                                                % driver dynamics
switch 2
    case 0, figure(44); vn = 1:2001;
    case 1, figure(45); vn = 20001:25001;
    case 2, figure(46); vn = 25001:150001;
end
        
comp_legend = ["ODC",'4th Bessel','GMTO noFF'];
ax_labels = ["AZ","EL","GIR"];
for i_ = 1:3
    subplot(3,2,2*i_-1)
    plot(r.t(vn),r.MotTo(vn,i_)); hold on;
    plot(mount.time(vn), mountTo(vn,i_),...
        mount_.time(vn), mountTo_(vn,i_),'--');
    ylabel(sprintf('%s DRV TO',ax_labels(i_)));
    if(i_ == 1), legend(comp_legend); legend boxoff; end
    if(i_ == 3), xlabel('Time (s)'); end
    hold off; axis tight
    subplot(3,2,2*i_)
    plot(r.t(vn),r.LoadPo(vn,i_)-r.SetPo(vn,i_)); hold on;
    e = mount.signals.values(vn,2*i_-1) - mount.signals.values(vn,2*i_);
    e_ = mount_.signals.values(vn,2*i_-1) - mount_.signals.values(vn,2*i_);
    plot(mount.time(vn),e*180/pi, mount_.time(vn),e_*180/pi,'-.');
    ylabel(sprintf('%s Pos (deg)',ax_labels(i_)));
    if(i_ == 1), legend(comp_legend,'Location','southeast'); legend boxoff;
    end
    if(i_ == 3), xlabel('Time (s)'); end
    hold off; axis tight
end

if true
    figure(47); i_ = 1;
    plot(r.t(vn),r.MotTo(vn,i_), mount.time(vn), mountTo(vn,i_)); 
    ylabel(sprintf('%s DRV TO (Nm)',ax_labels(i_)));
    legend("ODC",'GMTO','Fontsize',11); legend boxoff;
    xlabel('Time (s)');
    grid on; axis tight
end

%% Compare piston and TT

vn = 8001:179001;
figure(111)
subplot(211)
plot(r.t(vn),r.Piston(vn,:)); axis tight; grid on;
title("Seg piston (PMT2*1e9) - without RTF filtering")
ylabel("ODC Piston (nm)")
subplot(212)
% plot(t(vn),1e9*fun_detrend_v1(pmt1pmt2(vn,15:21),1)); axis tight; grid on;
plot(t(vn),1e9*pmt1pmt2(vn,15:21)); axis tight; grid on;
ylabel("GMTO Piston (nm)")
xlabel("Time (s)")
% title("Seg piston (PMT2) - wo RTF & linear contribution discarded")

%% WFE

fprintf('\nODC data and code from GMTO:\n');
t0_ = r.t;
vn = 25001:length(t0_);

if(~exist('hLtaoP','var'))
    hLtaoP = fun_rtf('LTAOPiston', false, false);
end
if(~exist('nTrend','var')), nTrend = 0; end

% p_detrend = fun_detrend_v1(r.Piston(1:end,:), nTrend);
% rtfPiston = fun_applyRtfAndCut(t0_, p_detrend ,hLtaoP, vn);
rtfPiston = fun_applyRtfAndCut(t0_, r.Piston(1:end,:) ,hLtaoP, vn);


if(~exist('kt1','var')), kt1 = 25001; end

[wfe_,TTwfe_,Pwfe_] = calc_WFE(0*1e3*r.TipTilt(vn,:), rtfPiston);
fprintf('PMT-WFE: %.3gnm\t%.3gnm\t%.3gnm\t\n\n',TTwfe_,Pwfe_,wfe_);












%% Function to post-process piston and TT as GMT-DOC-02011
%%
function [wfe, wfeTT, wfeP] = calc_WFE(segTT,segPiston)

if (size(segTT,2) ~= 14) || (size(segPiston,2) ~= 7)
    error("Tip&tilt and piston data at instant k shall be stored in different columns!")
elseif (size(segTT,1) ~= size(segPiston,1))
    error("Tip&tilt and piston data time sampling does not match!")
else
    N = size(segTT,1);
%     fprintf("Processing piston, tip&tilt data with %d samples.\n",N);
end

wfeTT = 10.2* sqrt((1/7/N)* sum((segTT - mean(segTT,2)).^2,'all'));

% wfeP = sqrt((1/7/N)* sum((segPiston - mean(segPiston,2)).^2,'all'));
piston_vRMS = fun_calcSimpleRms(segPiston);
wfeP = sqrt(sum(piston_vRMS.^2)/7);
info = 'Rtf Piston [nm] = [%.2f %.2f %.2f %.2f %.2f %.2f %.2f] => RMS=%.2f';
fprintf("\n"+info+"\n", piston_vRMS, wfeP);
wfe = abs(complex(wfeTT,wfeP));

end

%% Function to compute Hkin based on a particular model 
%%
function [Hkin_m12, Hkin_pmtnodes] = compute_Hkin(ModelFolder)

fprintf('Wait... Computing Hkin from model\n%s\n',ModelFolder);
% Load model
load(fullfile(ModelFolder,"modal_state_space_model_2ndOrder.mat"),...
    'outputTable','modalDisp2Outputs');
% M1/M2 RBM
m1_out = outputTable{'OSS_M1_lcl','indices'}{1};
m2_out = outputTable{'MC_M2_lcl_6D','indices'}{1};
pmt_out = outputTable{'PMT_3D','indices'}{1};
% MNT enc output indices
mntAZ_out = outputTable{'OSS_AzEncoder_Angle',"indices"}{1};
mntEL_out = outputTable{'OSS_ElEncoder_Angle',"indices"}{1};
mntGIR_out = outputTable{'OSS_RotEncoder_Angle',"indices"}{1};

Phi = modalDisp2Outputs;

Phi_mnt_enc = [mean(Phi(mntAZ_out,1:3),1);...
    mean(Phi(mntEL_out,1:3),1);...
    mean(Phi(mntGIR_out,1:3),1)];

Hkin_m12 = [Phi(m1_out,1:3); Phi(m2_out,1:3)] * pinv(Phi_mnt_enc);
Hkin_pmtnodes = Phi(pmt_out,1:3) * pinv(Phi_mnt_enc);
end


%% Auxiliar functions providing the LTAO RTFs
%%
% Function to compute the GLAO RTF (according to REQ-L3-OAD-35398)
function [rtf,L] = oad_AOTT_rtf() %#ok<DEFNU> 

% LTAO RTF (according to REQ-L3-OAD-35398)
fz = 800;       % ASM closed-loop bandwidth (Hz)
delta = 0.75;   % ASM control loop damping
T = 2e-3;       % TT loop sampling period
tau = 1e-3;     % TT sensor delay
gi = 0.4;       % TT "integral gain"

omz = fz*2*pi;
s=tf('s');
L = exp(-tau*s) * (1-exp(-T*s))/(T*s) * (gi/T/s) *...
    (omz^2)/((omz^2)+2*omz*delta*s+(s^2));

rtf = 1/(1 + L);

end

%% Function to compute the differential piston RTF
%%
function RTF_DP = getRTF_DP(rtf_variant,Tp,taup)
% rtf_variant == 0 -> GMT-REQ-00506 (Rev.L) - Eq.(3.2)
% rtf_variant == 1 => OIWFS variants
s=tf('s');

if(nargin < 1), rtf_variant = 0; end
% Optical Piston sensor Integration time
if(nargin < 2 && ~rtf_variant), Tp = 30; end
if(nargin < 2 && rtf_variant), Tp = 100e-3; end
% Optical piston sensor latency
if(nargin < 3 && ~rtf_variant), taup= 6e-3; end
if(nargin < 3 && rtf_variant), taup= 50e-3; end

% Edge sensor latency
switch rtf_variant
    case 0, taue = .1/1000;      % AGWS RTF
    otherwise, taue = .2e-3;     % OIWFS RTF
end
        
omz = 800*2*pi; % ASM closed-loop bandwidth (Hz)
dz = .75;       % ASM control loop damping
gpi = .5;       % Optical piston feedback integrator gain
Te = 2e-3;      % Edge sensor integration time
geff = .8;      % Edge sensor feedforward gain

num = 1 + omz^2*geff*exp(-taue*s)*(exp(-Te*s)-1)/(Te*s*    (omz^2+2*omz*dz*s+s^2) );
den = 1 - omz^2*gpi *exp(-taup*s)*(exp(-Tp*s)-1)/(Tp^2*s^2*(omz^2+2*omz*dz*s+s^2) );
RTF_DP = num/den;
    
end

%% ODC's RTF function
%%
function h=fun_rtf(sWhat, bPade, bWithPlot)
% fun_rtf  Builds one of 3 RTFs as a h(s) transfer function  
% sWhat: one of { 'NS','LTAO','LTAOPiston'}
% bPade=[true|false]: if true the Pade approximation is used for the time delay
% bWithPlot=[true|false]: if true tranfer function will be ploted

  vNames ={ 'NS','LTAO','LTAOPiston'};
  nWhat = find(strcmp(vNames, sWhat));
  if(isempty(nWhat)), error([sWhat ' is not defined!']), end
  if(1==nWhat)
    fz=25; dz=0.60;
    tau=6e-3; T=5e-3; g=0.2; bIsBase=true;
    h1=fun_tfBase(fz, dz, tau, T, g, bIsBase);  
    h = 1/(1-h1);  

  elseif(2==nWhat)
    fz=800; dz=0.75;
    tau=5e-4; T=2e-3; g=0.4; bIsBase=true;
    h1=fun_tfBase(fz, dz, tau, T, g, bIsBase);  
    h = 1/(1-h1);  

  elseif(3==nWhat)
    fz=800;dz=0.75;
    
    tau=1e-4; T=2e-3; g=0.8; bIsBase=false;
    h1=fun_tfBase(fz, dz, tau, T, g, bIsBase);  

    tau=6e-3; T=30;  g=0.5; bIsBase=true;
    h2=fun_tfBase(fz, dz, tau, T, g, bIsBase);  
    
    h=(1+h1)/(1-h2);
  else
    error('nWhat is wrong!')
  end
  
  if(bPade)
    h= pade(h,10);
  end
  
  if(bWithPlot)
    %disp(['Use Pade: ' num2str(bPade)])
    fun_bode_v1(h,100,'off',1,sWhat,'dB',logspace(-3,3,1000));
  end
end

%% **************************************************************
function h=fun_tfBase(f,d,tau,T,g,bIsBase)
  s=tf('s');
  x1 = 2*pi*f;
  x2 = x1^2 + 2*x1*d*s + s^2;
  x3 = x1^2 + x1*d*s + s^2;
  if(~bIsBase)
    h = (x1^2*g*exp(-tau*s)*(exp(-T*s) - 1))/(T*s*x3);
  else
    h = (x1^2*g*exp(-tau*s)*(exp(-T*s) - 1))/(T^2*s^2*x2);
  end  
end

%% *********************************************
function X=fun_applyRtfAndCut(t0_,X,hRtf,vn)
%fun_applyRtfAndCut Apply RTF to all columns and truncate according to index vector vn
    m=size(X,2); Hrtf=ss(zeros(m,m)); 
    for k=1:m, Hrtf(k,k) = hRtf; end
    X = lsim(Hrtf,X,t0_,[],'foh');
    X=X(vn,:);
end

%% *********************************************
function r_=fun_calcSimpleRms(X)
  m=size(X,2);
  n=size(X,1);
  for k=1:m
    r_(k)=sqrt(sum((X(:,k)).^2)/n); %#ok<AGROW> 
  end
end


