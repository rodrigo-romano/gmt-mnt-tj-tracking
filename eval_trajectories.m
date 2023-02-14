%
% Plot mount set-point trajectories 
% extrapolate to 1000Hz and plot
%
% Flag to remove the initial offset
soft_start = false;
% Trajectory label
traj_ID = '101';
% Trajectory file folder
ffolder = '/Users/rromano/Workspace/gmt-data/mnt_tj_dt/GMT-DTA-193917';
fname = fullfile(ffolder,...
    sprintf('GMT-DTA-193917_RevA_Trajectory%s.csv',traj_ID));
tj_DT = csvread(fname,10,0);

%%

t = tj_DT(:,1);

% dt = t(2)-t(1); % extrapolated sampling time.
dt = 1e-3; % extrapolated sampling time.
ts = dt; % start time
te = t(end);

t1 = (ts:dt:te);
N = length(t1);

az = zeros(N,1);
el = az;
gir = az;

for k = 1:N
   [~,id] = max(t(t<(t1(k)-1e-3)));
   if(isempty(id)), id = 1; end
   
   az0 = tj_DT(1,2);
   el0 = tj_DT(1,5);
   gir0 = tj_DT(1,8);
   
   delt = t1(k)-t(id);
   az(k) = tj_DT(id,2:4)*[1;delt;delt^2] - soft_start*az0;
   el(k) = tj_DT(id,5:7)*[1; delt; delt^2] - soft_start*el0;
   gir(k) = tj_DT(id,8:10)*[1; delt; delt^2] - soft_start*gir0;
   
end

%% ODC trajectories

try
    sDir = '/Users/rromano/Workspace/mnt-odc/2022-12-20_ODC E2E Files';
    [mTj,~,odc_az0,odc_el0,odc_gir0] = ...
        fun_loadTjData(sDir,traj_ID,ts,te);
    H_tj = tf(1,[2*pi*0.1, 1]);
    odc_tj = true;
catch
    odc_tj = false;
end


%% Position command trajectory

figure(1000 + str2double(traj_ID));

subplot(311); hold off;
stairs(t1,az,'k'); hold on;
stairs(t1,interp1(t,tj_DT(:,2),t1) - soft_start*az0,'.-');
if(odc_tj)
    stairs(mTj(:,1),mTj(:,2)+(~soft_start)*odc_az0,'--');
    tjf = lsim(H_tj,mTj(:,2),mTj(:,1));
    stairs(mTj(:,1),tjf+(~soft_start)*odc_az0,'k-.');
end
xlim([min(t1) 10]);
xlabel('time (seconds)'); ylabel('AZ deg');
grid on;
title('Azimuth axis');
if(odc_tj), legend('GMTO','GMTO cte term','ODC');
else, legend('GMTO','GMTO cte term');
end
legend('box','off','orientation','Horizontal','Location','southeast');

subplot(312); hold off;
stairs(t1,el,'k'); hold on;
stairs(t1,interp1(t,tj_DT(:,5),t1)- soft_start*el0,'.-');
if(odc_tj)
    stairs(mTj(:,1),mTj(:,3)+(~soft_start)*odc_el0,'--');
    tjf = lsim(H_tj,mTj(:,3),mTj(:,1));
    stairs(mTj(:,1),tjf+(~soft_start)*odc_el0,'k-.');
end
% xlim([min(t1) max(t1)]);
xlabel('time (seconds)'); ylabel('EL deg');
grid on;
title('Elevation axis')

subplot(313); hold off;
stairs(t1,gir,'k'); hold on;
stairs(t1,interp1(t,tj_DT(:,8),t1)- soft_start*gir0,'.-');
if(odc_tj)
    stairs(mTj(:,1),mTj(:,4)+(~soft_start)*odc_gir0,'--');
    tjf = lsim(H_tj,mTj(:,4),mTj(:,1));
    stairs(mTj(:,1),tjf+(~soft_start)*odc_gir0,'k-.');
end
% xlim([min(t1) max(t1)]);
xlabel('time (seconds)'); ylabel('GIR deg');
grid on;
title('GIR')


%% Approximation error plot

err = [az, el, gir] - interp1(t,tj_DT(:,[2,5,8]),t1);
if(1)
    figure(3000 + str2double(traj_ID));
    
    ax_label = {'Azimuth axis','Elevation axis','GIR'};
    for iax = 1:3
        subplot(3,1,iax); hold off;
        stairs(t1,err(:,iax),'k'); hold on;
        xlim([min(t1) 3]);
        xlabel('time (seconds)'); ylabel('Approx error (deg)');
        grid on;
        title(ax_label{iax});
    end
end
% if(odc_tj), legend('GMTO','GMTO cte term','ODC');
% else, legend('GMTO','GMTO cte term');
% end
% legend('box','off','orientation','Horizontal','Location','southeast');



%% Derivative of the trajectories

iax = 1;
% Original ODC's trajectory data
load(fullfile('/Users/rromano/Workspace/mnt-odc/2022-12-20_ODC E2E Files',...
    'res_eval/tracking/input/v20.9/tj',...
    'HcTp19/za30_v20.9_HcTp19_wlc0_tj101_nl1.mat'),'r');
% Numerical derivative
w=2*pi*200;  s=tf('s'); h=s/(s/w+1);
rvel = lsim(h, r.SetPo(:,iax), r.t); 
racc = lsim(h, rvel, r.t); 
tj = [r.SetPo(:,iax), rvel, racc];

% Alternative trajectory
switch 1
    case 0
        % 1st order finite-differences
        rvel_ = [0; diff(r.SetPo(:,iax))];
        racc_ = [0; 0; diff(r.SetPo(:,iax), 2)];        
        Ts_ = diff(r.t(1:2));
        tj_ = [r.SetPo(:,iax), (1/Ts_)*rvel_, (1/Ts_)^2*racc_];
    case 1
        % 4th-order Bessel filter
%         Hb4 = tf(2.4937, [1, 3.9257, 6.9349, 6.3522, 2.4937]); %0.2Hz
        Hb4 = tf(97.4091, [1, 9.8141, 43.3429, 99.2538, 97.4091]); %0.5Hz
%         Hb4 = tf(1.5585e3, [1, 19.6283, 173.3715, 794.03, 1.5585e3]); %1Hz
        rpos_ = lsim(Hb4, mTj(:,iax+1), mTj(:,1));
        rvel_ = lsim(h, rpos_, mTj(:,1));
        racc_ = lsim(h, rvel_, mTj(:,1));
        tj_ = [rpos_, rvel_, racc_];
        
end


figure(23);
vn = 25001:140001; %18001;%
plot_ylabel_units = ["Tj POS (deg)","Tj VEL (deg/s)","Tj ACC (deg/s^2)"];
ax_labels = ["AZ","EL","GIR"];
for ik = 1:3
    subplot(3,1,ik);
    plot(r.t(vn),tj(vn,ik),'Linewidth',2); hold on;
    plot(r.t(vn),tj_(vn,ik),'-.');
    if(ik == 1), plot(mTj(vn,1),mTj(vn,iax+1),'k--'); end
    ylabel(plot_ylabel_units(ik))
    grid on; axis tight; hold off;
    if(ik == 1)
        title(sprintf("tj101 - %s trajectory and its derivatives",ax_labels(iax)));
        legend('ODC data','4th-order Bessel @0.5Hz','raw tj data');
        legend boxoff
    end
    if(ik == 3)
        legend('ODC data','4th-order Bessel @0.5Hz'); legend boxoff
        xlabel("Time (s)");
    end
end

% create a new pair of axes inside current figure
axes('position',[.16 .135 .3 .12])
box on % put box around new pair of axes
vn = 86001:93001; %18001;%
plot(r.t(vn),tj(vn,ik),'Linewidth',2); hold on;
plot(r.t(vn),tj_(vn,ik),'-.');
set(gca,'Fontsize',8);
grid on; axis tight; hold off;
axis tight;






%% ODC Functions - Taken from fun_runOne.m
%%
%% ********************************************
function [mTj,mTjInit,az0,el0,gir0] = fun_loadTjData(sDir,sidTj,t1,t2)
  load([sDir '/sim/input/tj_tracking/tg_tj' sidTj '.mat'],'W','Wtg');
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