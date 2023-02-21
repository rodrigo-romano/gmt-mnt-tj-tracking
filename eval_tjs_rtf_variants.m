%
% eval_tjs_rtf_variants.m
% 


%% Preamble
%%
clearvars
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

% Load PMTs
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

% Load Segment TT optical sensitivity matrix
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


%% Load simulation data
%%
clear mountY
sim_label = 'za00';       % M1OFL:OFF, M2FSM/POS:OFF
% sim_label = 'za30';       % M1OFL:OFF, M2FSM/POS:OFF
% sim_label = 'm1ofl_m2fsm';  % M1OFL:ON, M2FSM/POS:ON
switch sim_label
    case 'm1ofl_m2fsm'
        FEM_LABEL = "20220611_1945_MT_mount_zen_00_m1HFN_FSM_";
        load("tj101_za00_m1ofl_m2fsm.mat",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj102_za00_m1ofl_m2fsm.mat",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
    case 'za00'
        FEM_LABEL = "20220611_1945_MT_mount_zen_00_m1HFN_FSM_";
        load("tj101_za00_mount_withff",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj102_za00_mount_noff",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
    case 'za30'
        FEM_LABEL = "20220610_1023_MT_mount_zen_30_m1HFN_FSM_";
%         load("tj101_za30_mount_noff.mat",'mount','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
        load("tj103_za30_mount_noff.mat",'mountY','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
%         load("tj101_za30_mount_withff.mat",'mount','PMT1_DT','PMT2_DT','M1rbm','M2rbm');
end

try mount = mountY;
catch, fprintf('Variable mountY is not available!\n');
end

% Compute Hkin
if( ~exist('Hkin','var')|| ~exist('Hkin_pmt','var')|| recomputeHkin || 0)
    ModelFolder = fullfile(im.lfFolder, FEM_LABEL);
    [Hkin, Hkin_pmtnodes] = compute_Hkin(ModelFolder);
    Hkin_pmt = [pmt1;pmt2]*Hkin_pmtnodes;
end

t = mount.time;     % Time vector [s]
Ts = t(2)-t(1);


%% Kinematic compensation
%%

% M1 and M2 RBM after kinematic compensation
m1m2RBM = [M1rbm,M2rbm] - mount.signals.values(:,1:2:6) * Hkin';
% PMT1/2 outputs after kinematic compensation
pmt1pmt2 = [PMT1_DT,PMT2_DT] - mount.signals.values(:,1:2:6) * Hkin_pmt';


%% Apply LOM transformations
%%
segtt = m1m2RBM * D_seg_tt';
segp = m1m2RBM * D_seg_piston';
segdp = segp - mean(segp,2);


%% Apply Rejection Transfer Functions (RTFs)
%%

% RTF comparison plots

AGWSrtf = getRTF_DP(0); %#ok<*UNRCH>
OIWFS10rtf = getRTF_DP(1,1/10,50e-3);
OIWFS50rtf = getRTF_DP(1,1/50,10e-3);
OIWFS125rtf = getRTF_DP(1,1/125,10e-3);
if true
    hbode = bodeoptions; hbode.FreqUnits = 'Hz'; hbode.grid = 'on';
    hbode.XLabel.FontSize = 12; hbode.YLabel.FontSize = 12;
    hbode.TickLabel.FontSize = 12; hbode.Title.String = '';
    figure(figidx + 11);
    bodemag(AGWSrtf,OIWFS10rtf,OIWFS50rtf,OIWFS125rtf,...
        2*pi*logspace(-3,3,501),hbode);
    legend('AGWS RTF','OIWFS-10 RTF','OIWFS-50 RTF','OIWFS-125 RTF',...
        'Location','southeast'); legend boxoff;
end


segdpf = fun_applyRtfAndCut(t, segdp ,AGWSrtf, 1:length(t));
meandpfAGWS = mean(segdpf,2);
segdpf_rms(:,1) = rms(segdpf ,2);%- mean(segdpf,2)

segdpf = fun_applyRtfAndCut(t, segdp ,OIWFS10rtf, 1:length(t));
meandpfOIWFS10 = mean(segdpf,2);
segdpf_rms(:,2) = rms(segdpf ,2);%- mean(segdpf,2)

segdpf = fun_applyRtfAndCut(t, segdp ,OIWFS50rtf, 1:length(t));
meandpfOIWFS50 = mean(segdpf,2);
segdpf_rms(:,3) = rms(segdpf ,2);%- mean(segdpf,2)

segdpf = fun_applyRtfAndCut(t, segdp ,OIWFS125rtf, 1:length(t));
meandpfOIWFS125 = mean(segdpf,2);
segdpf_rms(:,4) = rms(segdpf ,2);%- mean(segdpf,2)

segdp_rms = rms(segdp,2); %- mean(segdp,2)

% pmt2f = fun_applyRtfAndCut(t, pmt1pmt2(:,15:21) ,rtfP, 1:length(t));

%%
if true
    figure(figidx+21)
    vn = 25*(1/Ts):length(t);
    plot(t(vn),[meandpfAGWS(vn),meandpfOIWFS10(vn),meandpfOIWFS50(vn),meandpfOIWFS125(vn)]);
    axis tight; grid on;
    ylabel('Avg seg Piston ({\bf after RTF})')
    legend("AGWS","OIWFS10","OIWFS50","OIWFS125","wo RTF",...
        'Location','southwest','NumColumns',4); legend boxoff;
end


if true
    figure(figidx+22)
    vn = 25*(1/Ts):48000;%length(t);
    plot(t(vn),1e9*segdpf_rms(vn,:));
    axis tight; grid on;
    ylabel('Seg piston (nm)')
    legend("AGWS","OIWFS10","OIWFS50","OIWFS125",...
        'Location','southwest','NumColumns',4); legend boxoff;%,"wo RTF"
end

%% Compute PSDs
%%
% Length of data segment
nFFT = 2^19;
M = 2^17;
detrend = 'none'; %'mean'; %

PpsdLOMf_rms = zeros(max(nFFT,M)/2+1,4);
[PpsdLOM_rms,freqP] = utils.pwelch(1e9*segdp_rms(:,1),M,[],nFFT,1/Ts,'onesided',detrend);
[PpsdLOMf_rms(:,1),~] = utils.pwelch(1e9*segdpf_rms(:,1),M,[],nFFT,1/Ts,'onesided',detrend);
[PpsdLOMf_rms(:,2),~] = utils.pwelch(1e9*segdpf_rms(:,2),M,[],nFFT,1/Ts,'onesided',detrend);
[PpsdLOMf_rms(:,3),~] = utils.pwelch(1e9*segdpf_rms(:,3),M,[],nFFT,1/Ts,'onesided',detrend);
[PpsdLOMf_rms(:,4),~] = utils.pwelch(1e9*segdpf_rms(:,4),M,[],nFFT,1/Ts,'onesided',detrend);


%% Plot results
%%
if 1
    figure(figidx+19)
    set(gcf,'position',[123   230   1.3*640   400])
    subplot(121)
    loglog(freqP,PpsdLOMf_rms,'-');
    title("Segment Piston RMS PSD (GMT-LOM)")
    grid on; hold on; set(gca,'FontSize',12);
    loglog(freqP,PpsdLOM_rms,'k--');
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('RMS Seg Piston (nm^2/Hz)')
    xlim([min(freqP) max(freqP)]); axis tight;
    legend("AGWS","OIWFS10","OIWFS50","OIWFS125","wo RTF",...
        'Location','southwest','NumColumns',1);
    legend boxoff

    subplot(122)    
%     loglog(freqP,PpsdPMT,'-');
    deltaFs = diff(freqP(1:2));
    semilogx(freqP, sqrt(deltaFs*cumsum(PpsdLOMf_rms)));
    title("Cumulative RMS Piston PSD")
    grid on; hold on; set(gca,'FontSize',12);
    
    hold off;
    xlabel('Frequency (Hz)')
    ylabel('Seg Piston (nm^2/Hz)')
    xlim([min(freqP) max(freqP)]); axis tight;
    legend("AGWS","OIWFS10","OIWFS50","OIWFS125",...
        'Location','southwest','NumColumns',1);
    legend boxoff
end


%%
rmpath(odc_base_folder);
return









%% Function to post-process piston and TT as GMT-DOC-02011
%%
function [wfe, wfeTT, wfeP] = calc_WFE(segTT,segPiston) %#ok<*DEFNU> 

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

%% Auxiliar functions providing the LTAO RTFs
%%
% Function to compute the GLAO RTF (according to REQ-L3-OAD-35398)
function [rtf,L] = oad_AOTT_rtf() 

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
    case 0, taue = .1e-3;      % AGWS RTF
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


