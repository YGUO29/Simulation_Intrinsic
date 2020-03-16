% depth sensitivity, Sassaroli 2011 method
% uses the pMC function
% ==========================================
f2 = cell(1,length(ind));
det2 = f2;
traj2 = f2;
NewDetWeight_depth = f2;
NewDetWeight_hor = f2;
DetWeight = zeros(1,6);
profile_depth = f2;
profile_hor = cell(1,6);
SimVol = 200; % "replay" simulation volume, does not have to be the entire volume (to save time)
SimCenter = cfg.size(1)/2;

%% plot along z axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
z = [1:SimVol].*cfg.unitinmm;
profile_depth_peak = zeros(1:6);
for i = 1:6
    [maxy,ii] = max(profile_depth{i}); profile_depth_peak(i) = z(ii);
    % h1(i) = plot(z, profile_depth{i}, 'color',COLOR(i,:)); hold on
    h1(i) = plot(z, profile_depth{i}./maxy, 'color',COLOR(i,:)); hold on
    set(gca,'XLim',[0 4], 'YLim',[0 1])
    xlabel('distance from surface (mm)')
    ylabel('\DeltaSignal (%) when \Delta\mu_a=10%')
end

for k = 1:length(bounds1)
plot([bounds1(k) bounds1(k)].*cfg.unitinmm,get(gca,'YLim'),'LineStyle',':','color',[0.5 0.5 0.5]);
end

s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,6); s2(:) = {'nm, peak at '};
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeastoutside')
%% normalize with delta mua (yaxis = percent change in signal when mua increase by 1)
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
z = [1:SimVol].*cfg.unitinmm;
profile_depth_norm = profile_depth;
for i = 1:6
    if strcmp(cfg.config,'without skull')
        delta_mua = [0.1*skull_mua(ind(i)).*ones(1,bounds1(1)), 0.1*brain_mua(ind(i)).*ones(1,SimVol-bounds1(1))];
    else
        delta_mua = [0.1*skull_mua(ind(i)).*ones(1,bounds1(1)), 0.1*brain_mua(ind(i)).*ones(1,SimVol-bounds1(1))];
    end
    profile_depth_norm{i} = profile_depth_norm{i}./delta_mua;
    [maxy,ii] = max(profile_depth_norm{i}); profile_depth_peak(i) = z(ii);
%     h1(i) = plot(z, profile_depth_norm{i}./maxy, 'color',COLOR(i,:)); hold on
    h1(i) = plot(z, profile_depth_norm{i}, 'color',COLOR(i,:)); hold on
    set(gca,'yscale','log')
set(gca,'XLim',[0 4])
xlabel('distance from surface (mm)')
ylabel('\DeltaSignal (%) / \Delta\mu_a')

end

for k = 1:length(bounds1)
plot([bounds1(k) bounds1(k)].*cfg.unitinmm,get(gca,'YLim'),'LineStyle',':','color',[0.5 0.5 0.5]);
end

s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,6); s2(:) = {'nm, peak at '};
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeastoutside')
%% plot along x axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
x = [-99:100].*cfg.unitinmm;
profile_hor_dia = zeros(1,6);
for i = 2
    profile = profile_hor{i};
    h2(1) = plot(x, profile./max(profile), 'color',COLOR(i,:)); hold on
%     [~, ii] = min(abs((profile./max(profile) - 1/exp(1))));
    [~, ii] = min(abs((profile./max(profile) - 1/2)));
    profile_hor_dia(i) = 2*abs(x(ii));
    
    profile = profile_hor_straight{i};
    h2(2) = plot(x, profile./max(profile), 'LineStyle','--', 'color',COLOR(i,:)); hold on
%     [~, ii] = min(abs((profile./max(profile) - 1/exp(1))));
    [~, ii] = min(abs((profile./max(profile) - 1/2)));
    profile_hor_dia2 = 2*abs(x(ii));

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
end
s2 = cell(1,2); 
s2(1) = {'Angled illumination (NA = 0.8), \newline 1/e width = '};
s2(2) = {'Straight illumination, \newline 1/e width = '};

s3 = arrayfun(@num2str,[profile_hor_dia(i), profile_hor_dia2],'UniformOutput',false);
legend(h2, strcat(s2,s3), 'location', 'northeastoutside')
if isnan(cfg.perturb)
    title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change 10%'],'fontsize',14);
else
    title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change ', num2str(cfg.perturb)],'fontsize',14);
end

% s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
% s2 = cell(1,6); s2(:) = {'nm, diameter = '};
% s3 = arrayfun(@num2str,profile_hor_dia,'UniformOutput',false);
% legend(h2, strcat(s1,s2,s3), 'location', 'northeastoutside')

%% plot example trajectories
traj = traj2{2};
AxLim = [1 300 0 60];
XINTRINSIC_SIM_PlotTraj(traj, bounds1, AxLim)
%% normalize with delta mua (yaxis = percent change in signal when mua increase by 1)
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
x = [-SimVol/2+1:SimVol/2].*cfg.unitinmm;
profile_hor_dia = zeros(1,6);
for i = 1:6
    profile = profile_hor{i};
    h2(i) = plot(x, profile./max(profile), 'color',COLOR(i,:)); hold on
    [~, ii] = min(abs((profile./max(profile) - 1/exp(2))));
    profile_hor_dia(i) = 2*abs(x(ii));
% set(gca,'yscale','log')
% set(gca,'XLim',[0 4])
    xlabel('distance from the center (mm)')
end
s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,6); s2(:) = {'nm, diameter = '};
s3 = arrayfun(@num2str,profile_hor_dia,'UniformOutput',false);
legend(h2, strcat(s1,s2,s3), 'location', 'northeastoutside')
%% import figure and plot data together
fig2 = gcf;
%%
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
x = [-99:100].*cfg.unitinmm;
% profile_hor_dia = zeros(1,6);
for i = 1:6
    axObjs = fig1.Children; 
    dataObjs = axObjs(2).Children;
    profile1 = dataObjs(i).YData;
    axObjs = fig2.Children; 
    dataObjs = axObjs(2).Children;
    profile2 = dataObjs(i).YData;
%     profile2 = profile_hor{i}./max(profile_hor{i});
subplot(2,3,i)
    h2(i) = plot(x, profile1, 'color',COLOR(i,:)); hold on
    h2(i) = plot(x, profile2, 'LineStyle','--', 'color',COLOR(i,:)); hold on

%     [~, ii] = min(abs((profile./max(profile) - 1/exp(1))));
%     profile_hor_dia(i) = 2*abs(x(ii));
% set(gca,'yscale','log')
% set(gca,'XLim',[0 4])
    xlabel('distance from the center (mm)')
end
% s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
% s2 = cell(1,6); s2(:) = {'nm, diameter = '};
% s3 = arrayfun(@num2str,profile_hor_dia,'UniformOutput',false);
% legend(h2, strcat(s1,s2,s3), 'location', 'northeastoutside')
