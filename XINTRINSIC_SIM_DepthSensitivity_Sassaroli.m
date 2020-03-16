% depth sensitivity, Sassaroli 2011 method
% ==========================================
f2 = cell(1,length(ind));
det2 = f2;
traj2 = f2;
NewDetWeight_depth = f2;
NewDetWeight_hor = f2;
DetWeight = zeros(1,6);
profile_depth = f2;
profile_hor = cell(1,6);
SimVol = 10; % "replay" simulation volume, does not have to be the entire volume (to save time)
SimCenter = cfg.size(1)/2;

%% get initial total photon weights
for i = 2


det_temp = det1{i};
det_temp.ppath = det_temp.ppath(det_temp.DetInd,:); % select photons by its exit angle
DetWeight(i) = sum(mcxdetweight(det_temp,cfg.prop));
% change properties
% switch cfg.config
%     case {'without skull', 'without skull off focus'}
%         if isnan(cfg.perturb)
%         cfg.prop = [0 0 1 1            % medium 0: the environment
%         brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
%         brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37   % medium 2: white matter
%         brain_mua(ind(i))*1.1 gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 3: gray matter, perturbed
%         brain_mua(ind(i))*1.1 white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 4: white matter, pertubed
%         else
%         cfg.prop = [0 0 1 1            % medium 0: the environment
%         brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
%         brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37   % medium 2: white matter
%         brain_mua(ind(i))+cfg.perturb gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 3: gray matter, perturbed
%         brain_mua(ind(i))+cfg.perturb white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 4: white matter, perturbed
%         end
%     case {'with skull', 'with skull off focus'}
%         if isnan(cfg.perturb)
%         cfg.prop=[0 0 1 1            % medium 0: the environment
%         skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
%         brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
%         brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37   % medium 3: white matter
%         skull_mua(ind(i))*1.1 skull_mus(ind(i)) 0.9337, 1.56 % medium 4: skull, perturbed
%         brain_mua(ind(i))*1.1 gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 5: gray matter, perturbed
%         brain_mua(ind(i))*1.1 white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 6: white matter, perturbed
%         else
%         cfg.prop=[0 0 1 1            % medium 0: the environment
%         skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
%         brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
%         brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37   % medium 3: white matter
%         skull_mua(ind(i))+cfg.perturb skull_mus(ind(i)) 0.9337, 1.56 % medium 4: skull, perturbed
%         brain_mua(ind(i))+cfg.perturb gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 5: gray matter, perturbed
%         brain_mua(ind(i))+cfg.perturb white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 6: white matter, perturbed
%         end
%     otherwise 
%         disp('Check property setting (cfg.prop)')
% end

           
% calculate total photon weights after perturbation   
newcfg              = cfg;
newcfg.seed         = seeds{i}.data(:,det1{i}.DetInd);
newcfg.detphotons   = det1{i}.data(:,det1{i}.DetInd);
% newcfg.outputtype   = 'jacobian';

%  ========== change properties in Z direction  ==========
NewDetWeight_depth{i} = zeros(1,SimVol);
for zz = 1:SimVol
    % ========== define tissue structure ==========
%     newcfg.vol = cfg.vol;
%     switch cfg.config
%         case 'without skull'
%             if zz <= bounds
%                 newcfg.vol(:,:,zz) = 3; % gray matter, perturbed
%             else
%                 newcfg.vol(:,:,zz) = 4; % white matter, perturbed
%             end
%             newcfg.vol=uint8(newcfg.vol);
%         case 'with skull'
%             if zz <= bounds(1)
%                 newcfg.vol(:,:,zz) = 4; % skull, perturbed
%             elseif zz > bounds(1) && zz <= bounds(2)
%                 newcfg.vol(:,:,zz) = 5; % gray matter, perturbed           
%             else
%                 newcfg.vol(:,:,zz) = 6; % white matter, perturbed
%             end
%             newcfg.vol=uint8(newcfg.vol);
%         otherwise
%     end
    tic, [f2{i}, det2{i}, vol2, seeds2, traj2{i}] = mcxlab(newcfg); toc
    det_temp = det2{i};
    NewDetWeight_depth{i}(zz) = sum(mcxdetweight(det_temp,newcfg.prop));
    
end
profile_depth{i} = -(NewDetWeight_depth{i} - DetWeight(i)).*100./DetWeight(i);

% ========== change properties in X direction  ==========
% NewDetWeight_hor{i} = zeros(1,SimVol);
% fwait = waitbar(0,'Perturbation simulation started ...');
% 
% for xx = SimCenter-SimVol/2+1: SimCenter+SimVol/2
%     waitbar((xx-(SimCenter-SimVol/2))/SimVol,fwait,['Perturbation simulation processing ',num2str((xx-(SimCenter-SimVol/2))),'/',num2str(SimVol)]);
%     
%     % ========== define tissue structure ==========
%     newcfg.vol = cfg.vol;
%     switch cfg.config
%         case 'without skull'
%             newcfg.vol(xx,SimCenter,1:bounds) = 3; % gray matter, perturbed
%             newcfg.vol(xx,SimCenter,bounds+1:end) = 4; % white matter, perturbed
%             newcfg.vol=uint8(newcfg.vol);
%         case 'without skull off focus'
%             newcfg.vol(xx,SimCenter,bounds(1)+1:bounds(2)) = 3; % gray matter, perturbed
%             newcfg.vol(xx,SimCenter,bounds(2)+1:end) = 4; % white matter, perturbed
%             newcfg.vol(xx,SimCenter,bounds+1:end) = 3; % gray matter, perturbed
%             newcfg.vol=uint8(newcfg.vol);
%         case 'with skull'
%             newcfg.vol(xx,:,1:bounds(1)) = 4; % skull, perturbed
%             newcfg.vol(xx,:,bounds(1)+1:bounds(2)) = 5; % gray matter, perturbed
%             newcfg.vol(xx,:,bounds(2)+1:end) = 6; % gray matter, perturbed
%             newcfg.vol(xx,SimCenter,bounds(1)+1:bounds(2)) = 5; % gray matter, perturbed
%             newcfg.vol(xx,SimCenter,bounds(2)+1:end) = 6; % gray matter, perturbed
%             newcfg.vol=uint8(newcfg.vol);
%         case 'with skull off focus'
%             newcfg.vol(xx,SimCenter,bounds(1)+1:bounds(2)) = 5; % gray matter, perturbed
%             newcfg.vol(xx,SimCenter,bounds(2)+1:end) = 5; % gray matter, perturbed
%             newcfg.vol=uint8(newcfg.vol);
%         otherwise
%     end
%     tic, [f2{i}, det2{i}, vol2, seeds2, traj2{i}] = mcxlab(newcfg); toc
%     det_temp = det2{i};
%     NewDetWeight_hor{i}(xx-(SimCenter-SimVol/2)) = sum(mcxdetweight(det_temp,newcfg.prop));
% end
% profile_hor{i} = -(NewDetWeight_hor{i} - DetWeight(i)).*100./DetWeight(i);


% ========== change properties in R direction  ==========
% NewDetWeight_r{i} = zeros(1,SimVol/2);
% fwait3 = waitbar(0,'Horizontal perturbation simulation started ...');
% theta = 0:pi/50:2*pi;
% profile_r{i} = zeros(1,SimVol/2+1);
% % ========== assign r values to x and y ==========
% R = zeros(SimVol/2, SimVol/2);
% for xx = 1:SimVol/2
%     for yy = 1:SimVol/2
%         R(xx,yy) = round(vecnorm([xx,yy]));
%     end
% end
% temp = 1:SimVol/2;
% RR = [rot90(R,2), flipud(temp'), flipud(R); fliplr(temp), 0, temp; fliplr(R), temp', R ];
% r_real = zeros(1, SimVol+1);
% test = zeros(400, 400);
% for rr = 0 : SimVol/2
%     waitbar(rr/(SimVol/2),fwait3,['Horizontal perturbation simulation processing ',num2str(rr),'/',num2str(SimVol/2)]);
%     % ========== define tissue structure ==========
%     newcfg.vol = cfg.vol;
%     [xx, yy, ~] = find(RR == rr); 
%     xx = xx - SimVol/2 - 1 + SimCenter; 
%     yy = yy - SimVol/2 - 1 + SimCenter; 
%     switch cfg.config
%         case 'without skull'
%             xx = rr.*cos(theta);
%             yy = rr.*sin(theta);
%             newcfg.vol(150+round(xx),150+round(yy),1:bounds) = 3; % gray matter, perturbed
%             newcfg.vol(150+round(xx),150+round(yy),bounds+1:end) = 4; % white matter, perturbed
%             newcfg.vol=uint8(newcfg.vol);
%         case 'without skull off focus'
%             for k = 1:length(xx)
% %                 newcfg.vol(xx(k),yy(k),1:bounds) = 3; % gray matter, perturbed
%                 newcfg.vol(xx(k),yy(k),bounds+1:end) = 3; % gray matter, perturbed
%                 test(xx(k),yy(k)) = rr;
%             end         
%         case {'with skull', 'with skull off focus'}
%             for k = 1:length(xx)
% %                 newcfg.vol(xx(k),yy(k),bounds(1)+1:bounds(2)) = 4; % skull matter, perturbed
%                 newcfg.vol(xx(k),yy(k),bounds(2)+1:end) = 5; % gray matter, perturbed
%                 test(xx(k),yy(k)) = rr;
%             end
%         otherwise
%     end
%     newcfg.vol=uint8(newcfg.vol);
%     xx = xx - SimCenter; yy = yy - SimCenter;
%     r_real(rr+1) = mean(vecnorm([xx,yy]'));

%     VolDiff = newcfg.vol(:,:,bounds(1)+1) - cfg.vol(:,:,bounds(1)+1);
%     nChange(rr+1) = length(find(VolDiff(:)));
%     nChange(rr+1) = length(xx); % how many columes were perturbed
%     tic, [f2{i}, det2{i}, vol2, seeds2, traj2{i}] = mcxlab(newcfg); toc
%     det_temp = det2{i};
%     NewDetWeight_r{i}(rr+1) = sum(mcxdetweight(det_temp,newcfg.prop));

end
% r_real = [-fliplr(r_real(2:101)), r_real(1:101)];
% profile_r{i} = -(NewDetWeight_r{i} - DetWeight(i))./nChange.*100./DetWeight(i);
% profile_r_skull{i} = [fliplr(profile_r{i}(2:end)) profile_r{i}];

%% plot along z axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
z = [1:SimVol].*cfg.unitinmm;

for i = 1:6
[maxy,ii] = max(profile_depth{i}); profile_depth_peak(i) = z(ii);
% h1(i) = plot(z, profile_depth{i}, 'color',COLOR(i,:)); hold on
h1(i) = plot(z, profile_depth{i}./maxy, 'color',COLOR(i,:)); hold on
% set(gca,'yscale','log')
set(gca,'XLim',[0 4], 'YLim',[0 1])
if i == 6
    for k = 1:length(bounds)
    plot([bounds(k) bounds(k)].*cfg.unitinmm,get(gca,'YLim'),'LineStyle',':','color',[0.5 0.5 0.5]);
    end
end
xlabel('distance from surface (mm)')
ylabel('\DeltaSignal (%) when \Delta\mu_a=10%')

end
s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,6); s2(:) = {'nm, peak at '};
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeastoutside')
%% normalize with delta mua (yaxis = percent change in signal when mua increase by 1)
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
z = [1:SimVol].*cfg.unitinmm;
profile_depth_norm = profile_depth1;
for i = 2
    if strcmp(cfg.config,'without skull')
        delta_mua = [0.1*skull_mua(ind(i)).*ones(1,cfg.bounds(1)), 0.1*brain_mua(ind(i)).*ones(1,SimVol-cfg.bounds(1))];
    else
        delta_mua = [0.1*skull_mua(ind(i)).*ones(1,cfg.bounds(1)), 0.1*brain_mua(ind(i)).*ones(1,SimVol-cfg.bounds(1))];
    end
    profile_depth_norm{i} = profile_depth_norm{i}./delta_mua;
    [maxy,ii] = max(profile_depth_norm{i}); profile_depth_peak(i) = z(ii);
%     h1(i) = plot(z, profile_depth_norm{i}./maxy, 'color',COLOR(i,:)); hold on
    h1(i) = plot(z, profile_depth_norm{i}, 'color',COLOR(i,:)); hold on
    set(gca,'yscale','log')
set(gca,'XLim',[0 4])
if i == 6
    for k = 1:length(cfg.bounds)
    plot([cfg.bounds(k) cfg.bounds(k)].*cfg.unitinmm,get(gca,'YLim'),'LineStyle',':','color',[0.5 0.5 0.5]);
    end
end
xlabel('distance from surface (mm)')
ylabel('\DeltaSignal (%) / \Delta\mu_a')

end
s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,6); s2(:) = {'nm, peak at '};
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeastoutside')
%% plot along x axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
x = [-SimVol/2:SimVol/2].*cfg.unitinmm;
profile_hor_dia = zeros(1,6);
for i = 2
    profile = profile_hor{i};
    h2(1) = plot(x, profile./sum(profile), 'color',COLOR(i,:)); hold on
%     [~, ii] = min(abs((profile./max(profile) - 1/exp(1))));
    [~, ii] = min(abs((profile./sum(profile) - max(profile./sum(profile))./2)));
    profile_hor_dia(i) = 2*abs(x(ii));
    
%     profile = profile_hor_straight{i};
%     h2(2) = plot(x, profile./max(profile), 'LineStyle','--', 'color',COLOR(i,:)); hold on
% %     [~, ii] = min(abs((profile./max(profile) - 1/exp(1))));
%     [~, ii] = min(abs((profile./max(profile) - 1/2)));
%     profile_hor_dia2 = 2*abs(x(ii));

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
    xlim([-0.4 0.4])
end
% s2 = cell(1,2); 
% s2(1) = {'Angled illumination (NA = 0.8), \newline 1/e width = '};
% s2(2) = {'Straight illumination, \newline 1/e width = '};
% 
% s3 = arrayfun(@num2str,[profile_hor_dia(i), profile_hor_dia2],'UniformOutput',false);
% legend(h2, strcat(s2,s3), 'location', 'northeastoutside')

i = 2;
s1 = arrayfun(@num2str,WLs(i),'UniformOutput',false);
s2 = cell(1,length(i)); s2(:) = {'nm, diameter = '};
s3 = arrayfun(@num2str,profile_hor_dia(i),'UniformOutput',false);
legend(h2, strcat(s1,s2,s3), 'location', 'northeastoutside')

if isnan(cfg.perturb)
    title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change 10%'],'fontsize',14);
else
    title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change ', num2str(cfg.perturb)],'fontsize',14);
end

%% plot along R axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
profile_r_dia = zeros(1,3);
for i = 2
    profile = profile_r5{i};
%     h2(i) = plot(x, profile./brain_mua(ind(i))/0.1, 'color',COLOR(i,:)); hold on
    h2 = plot(cfg.unitinmm.*r_real, profile, 'color',COLOR(i,:)); hold on
    [~, ii] = min(abs((profile - max(profile)./2)));
    profile_r_dia(i) = 2*abs(r_real(ii)).*cfg.unitinmm;

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
    xlim([-5 5])
end

i = 2;
s1 = arrayfun(@num2str,WLs(i),'UniformOutput',false);
s2 = cell(1,length(i)); s2(:) = {'nm, diameter = '};
s3 = arrayfun(@num2str,profile_r_dia(i),'UniformOutput',false);
legend(h2, strcat(s1,s2,s3), 'location', 'northeastoutside')

% s2 = cell(1,3); s2(:) = {'Straight illumination, \newline FWHM = '};
% s3 = arrayfun(@num2str,profile_r_dia ,'UniformOutput',false);
% legend(h2, strcat(s2,s3), 'location', 'northeastoutside')
if isnan(cfg.perturb)
    title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change 10%'],'fontsize',14);
else
    title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change ', num2str(cfg.perturb)],'fontsize',14);
end

%% plot example trajectories
traj = traj2{2};
AxLim = [1 300 1 100];
XINTRINSIC_SIM_PlotTraj(traj, cfg.bounds, AxLim)
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
fig{3} = gcf;
%%
profile = cell(1,3);

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
% profile_hor_dia = zeros(1,6);
for i = 1:3
    axObjs = fig{i}.Children; 
    dataObjs = axObjs(2).Children;
    profile{i} = dataObjs.YData;
    x = dataObjs.XData;
% subplot(2,3,i)
    h2(i) = plot(x, profile{i}, 'color',COLOR(i,:)); hold on
    
%     [~, ii] = min(abs((profile./max(profile) - 1/exp(1))));
%     profile_hor_dia(i) = 2*abs(x(ii));
% set(gca,'yscale','log')
set(gca,'XLim',[-0.4 0.4])
    xlabel('distance from the center (mm)')
end
% s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
% s2 = cell(1,6); s2(:) = {'nm, diameter = '};
% s3 = arrayfun(@num2str,profile_hor_dia,'UniformOutput',false);
% legend(h2, strcat(s1,s2,s3), 'location', 'northeastoutside')
