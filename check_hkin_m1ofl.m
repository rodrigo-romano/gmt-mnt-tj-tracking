%
% eval_tjs_rtf_variants.m
% 


%% Preamble
%%
% clearvars
% radians to mas conversion factor
rad2mas = 1e3*(180/pi)*3600;
% Root of the mean squared value function
rms = @(x,dir)squeeze(sqrt(mean(x.^2,dir)));
% Flag to recompute Hkin
recomputeHkin = true;%false;
% Figure number index
figidx = 0;

% Piston optical sensitivity matrix
if(~exist('D_seg_piston','var') || false)
    pistonOSMfile = fullfile(im.lfFolder,'LOM-data','D_seg_piston_dt.mat');
    load(pistonOSMfile,'D_seg_piston');
    fprintf('\nPiston sensitivity matrix loaded from \n%s\n',pistonOSMfile);
end

%% Load simulation data
%%

% Kinematic compensation matrix
FEM_LABEL = "20220611_1945_MT_mount_zen_00_m1HFN_FSM_";
% Compute Hkin
if( ~exist('Hkin_hp','var')|| recomputeHkin || 0)
    ModelFolder = fullfile(im.lfFolder, FEM_LABEL);
    [Hkin, ~, Hkin_hp, ~, ~] = compute_Hkin(ModelFolder);
end

% Exp1 - data
load("tj101_za00_m1ofl_m2fsm_debug","mountY",'M1rbm','M2rbm','m1HP_D');
t = mountY.time;     % Time vector [s]
Ts = t(2)-t(1);
relHP_D = (m1HP_D - mountY.signals.values(:,1:2:6)*Hkin_hp')...
    *kron(eye(7),[eye(6),-eye(6)])';
m1m2RBM = [M1rbm,M2rbm] - mountY.signals.values(:,1:2:6) * Hkin';
segp = m1m2RBM * D_seg_piston';
segdp = segp - mean(segp,2);

% Exp2 - data
% load("tj101_za00_mntFBFF_debug","mountY",'M1rbm','M2rbm','m1HP_D');
load("tj101_za00_mntFBFF_DCcomp_debug","mountY",'M1rbm','M2rbm','m1HP_D');
relHP_D_ = (m1HP_D - mountY.signals.values(:,1:2:6)*Hkin_hp')...
    *kron(eye(7),[eye(6),-eye(6)])';
m1m2RBM_ = [M1rbm,M2rbm] - mountY.signals.values(:,1:2:6) * Hkin';
segp_ = m1m2RBM_ * D_seg_piston';
segdp_ = segp_ - mean(segp_,2);


%% Plot results
%%

iseg = 2;
% vn = 25*(1/Ts):length(t);
vn = 1001:50*(1/Ts);
% M1 HP displacement
figure(figidx + iseg*10);
plot(t(vn),relHP_D(vn,(1:6)+(iseg-1)*6)); hold on;
set(gca,'ColorOrderIndex',1)
plot(t(vn),relHP_D_(vn,(1:6)+(iseg-1)*6),'--');
xlabel("Times (s)");
ylabel(sprintf("M1S%d - HP displacements",iseg));
hold off; % axis tight;
legend("HP1","HP2","HP3","HP4","HP5","HP6",...
        "HP1","HP2","HP3","HP4","HP5","HP6",...
        'Location','southeast','NumColumns',2);

% create a new pair of axes inside current figure
axes('position',[.4 .52 .45 .35])
box on % put box around new pair of axes
vn = 36001:38001;
plot(t(vn),relHP_D(vn,(1:6)+(iseg-1)*6)); hold on;
set(gca,'ColorOrder','factory')
plot(t(vn),relHP_D_(vn,(1:6)+(iseg-1)*6),':');
set(gca,'Fontsize',8);
grid on; axis tight; hold off;
axis tight;


%% M1 Tz motion
%%
vn = 58000:60001;
% M1-Tz displacement
figure(figidx + 1+ iseg*10);
plot(t(vn),1e6*m1m2RBM(vn,3:6:42)); hold on;
plot(t(vn),1e6*m1m2RBM_(vn,3:6:42),'--');
xlabel("Times (s)");
ylabel("M1 Tz (um)");
hold off; % axis tight;
legend("S1","S2","S3","S4","S5","S6","S7",...
        "S1","S2","S3","S4","S5","S6","S7",...
        'Location','southwest','NumColumns',2);
    

%% Segment differential piston plot
%%
vn = 36001:37001;
% M1 HP displacement
figure(figidx + 2+ iseg*10);
plot(t(vn),segdp(vn,:)); hold on;
plot(t(vn),segdp_(vn,:),'--');
xlabel("Times (s)");
ylabel("Seg piston");
hold off; % axis tight;
legend("S1","S2","S3","S4","S5","S6","S7",...
        "S1","S2","S3","S4","S5","S6","S7",...
        'Location','southwest','NumColumns',2);