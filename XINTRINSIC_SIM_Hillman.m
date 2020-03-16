%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cfg;
% set tissue parameters
load('D:\=code=\Simulation_Intrinsic\spectra.mat')
lambda = spectra(:,1);

WLs = [488 530 630];
brain_mua = [0.33 0.55 0.024];
gray_matter_mus = [21 21 24];
gray_matter_g = [0.8 0.82 0.87];
ind = [1 2 3];
% COLOR = hex2rgb(['1B97F7'; '23C249'; 'FAB738'; 'E92716'; '874040'; '000000']);
COLOR = [   0       0       1;
            0       0.5     0;
            1       0       0];
%% preparing the input data

    % set seed to make the simulation repeatible
    cfg.seed=hex2dec('623F9A9E'); 

    cfg.nphoton=1e7;
    cfg.unitinmm = 0.01; % define units, 10 um
    % time-domain simulation parameters
    cfg.tstart=0;
    cfg.tend=5e-9;
    cfg.tstep=5e-9;
    cfg.issrcfrom0 = 1;
%     cfg.outputtype = 'fluence';


% ========== define source and detector ==========
%     cfg.srctype = 'pencil';
    NA = 0.12;
    cfg.srctype = 'cone';
    cfg.srcparam1=[asin(NA) 0 0 0];
    % cfg.srcparam1 = 100; %radius for the disk
    cfg.size = 300.*ones(1,3); % 15mm cube
%     cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcpos=[149.5, cfg.size(2)/2, 1];

    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/sqrt(2), cfg.size(2)/sqrt(2))]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
% ========== define a 3 layer structure ==========
cfg.config = 'without skull';
switch cfg.config
    case 'without skull'
        cfg.vol = ones(cfg.size); % 15 mm volume
        cfg.vol = uint8(cfg.vol);
    case 'with skull'
    otherwise
end

    cfg.isreflect = 1; % disable reflection at exterior boundary
    % GPU thread configuration
    cfg.autopilot=1;
    cfg.gpuid=1;
% ========== output control ==========
    cfg.savedetflag = 'dpxv';
    cfg.debuglevel = 'M';
    cfg.maxjumpdebug = cfg.nphoton; % number of trajectory positions saved

%% forward simulation
f1 = cell(1,length(ind));
seeds = f1;
det1 = f1;
traj = f1;
Phi_s = f1;
%%
cfg.perturb = nan; % a number or nan (perturb by 10%)
cfg.exitangle = [0 0.03];
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')
% for i = 1:length(ind)
for i = 2:3
    switch cfg.config
        case 'without skull'
            if isnan(cfg.perturb)
            cfg.prop = [0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i))*1.1 gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37 ];  % medium 3: gray matter, perturbed
            else
            cfg.prop = [0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i))+cfg.perturb gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37 ];  % medium 3: gray matter, perturbed
            end           
       case 'with skull'
        otherwise 
            disp('Check property setting (cfg.prop)')
    end
   
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    % f1=mcxlab(cfg);
    [f1{i},det1{i},vol,seeds{i},traj{i}]=mcxlab(cfg);
    toc;
% ========== detected photons ==========
    ExitAngle           = vecnorm(det1{i}.v(:,1:2)'); % sin_theta (i.e. NA), v is a normal vector
    det1{i}.DetInd      = find(ExitAngle < cfg.exitangle(2) & ExitAngle >= cfg.exitangle(1)); 
    det1{i}.DetNumber   = length(det1{i}.DetInd);
    det1{i}.DetExitPos  = det1{i}.p(det1{i}.DetInd,1:2);
% ========== plot figure ==========
    subplot(1,3,i);
    I1 = log10(squeeze(sum(f1{i}.data(:,size(cfg.vol,1)./2+1,:,:),4))');
    I1(I1<0) = 0;
%         if i == 1
        cmax1 = max(I1(:)); cmin1 = min(I1(:));
%         end
    contourf(I1, cmin1:0.5:cmax1);
    axis square
    title(['flux with reflection, no skull, wavelength ',num2str(WLs(i)),'nm'],'fontsize',14);
    set(gca,'clim',[cmin1 cmax1]);
    colorbar
    drawnow
    
    Phi_s{i} = squeeze(sum(f1{i}.data,4));
end

%% reversed simulation
f2 = cell(1,length(ind));
seeds2 = f2;
det2 = f2;
traj2 = f2;
Phi_d = f2; % flux map from reversed simulation
% ========== switch source and detector, straight disk illumination ==========
    cfg.srctype = 'disk';
    cfg.srcparam1=[150*sqrt(2) 0 0 0];
    % cfg.srcparam1 = 100; %radius for the disk
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcdir=[0 0 1];
    cfg.srcdir=cfg.srcdir./vecnorm(cfg.srcdir);
%     cfg.srcdir=[0.8 0 0.6];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, 1]; % [x y z radius]

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
for i = 1:3
    switch cfg.config
        case 'without skull'
            if isnan(cfg.perturb)
            cfg.prop = [0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i))*1.1 gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37 ];  % medium 3: gray matter, perturbed
            else
            cfg.prop = [0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i))+cfg.perturb gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37 ];  % medium 3: gray matter, perturbed
            end           
       case 'with skull'
        otherwise 
            disp('Check property setting (cfg.prop)')
    end
   
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    % f1=mcxlab(cfg);
    [f2{i},det2{i},vol,seeds2{i},traj2{i}]=mcxlab(cfg);
    toc;
% ========== detected photons ==========
%     ExitAngle   = sin(vecnorm(det1{i}.v(:,1:2)')./vecnorm(det1{i}.v')); % sin_theta (NA)
%     det1{i}.DetInd      = find(ExitAngle < NA);
%     det1{i}.DetNumber   = length(det1{i}.DetInd);
%     det1{i}.DetExitPos  = det1{i}.p(det1{i}.DetInd,1:2);
% ========== plot figure ==========
    subplot(2,3,i);
%     I1 = log10(squeeze(sum(f1{i}.data(:,size(cfg.vol,1)./2+1,:,:),4))');
    I1 = squeeze(sum(f2{i}.data,4));
    I1 = log10(squeeze( sum(I1,2))');
    I1(I1<0) = 0;
        if i == 1
        cmax1 = max(I1(:)); cmin1 = min(I1(:));
        end
    contourf(I1, cmin1:0.1:cmax1);
%     imagesc(flipud(I1),[cmin, cmax])
    hold on
    for k = 1:length(bounds1)
        plot([0 size(cfg.vol,1)],[bounds1(k) bounds1(k)],'--');
    end
    title(['fluence without reflection, with skull, wavelength ',num2str(WLs(i)),'nm'],'fontsize',14);
    set(gca,'clim',[cmin1 cmax1]);
    colorbar
    drawnow
    
    Phi_d{i} = squeeze(sum(f2{i}.data,4));
end
%% Jacobian
fig3 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
fig4 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
fig5 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);


Jacobian = cell(1,3);
profile_depth_peak = zeros(1,3);
profile_hor_dia = zeros(1,3);
for i = 1:3
    
det_temp = det1{i};
det_temp.ppath = det_temp.ppath(det_temp.DetInd,:); % select photons by its exit angle
DetWeight(i) = sum(mcxdetweight(det_temp,cfg.prop));

Jacobian{i} = cfg.unitinmm.^3.*Phi_s{i}.*Phi_d{i} ./ (DetWeight(i)/cfg.nphoton);
Jacobian_2d = squeeze( sum(Jacobian{i},2) ); % sum across y direction

figure(fig3), subplot(2,3,i), 
Jacobian_temp = Jacobian_2d(131:170, 1:16);
cmax1 = max(Jacobian_temp(:)); cmin1 = min(Jacobian_temp(:));
% imagesc(log(flipud(Jacobian_temp')), [-17 -9])
contourf(log(flipud(Jacobian_temp')))
axis image
title(['Jacobian', num2str(WLs(i)),' nm']), colorbar

figure(fig4), hold on, 
xx = [1:300].*cfg.unitinmm;
% profile_depth = squeeze(sum(Jacobian_2d,1));
profile_depth = squeeze(sum(Jacobian_2d,1))./max(squeeze(sum(Jacobian_2d,1))); % sum across x direction
[~,ii] = max(profile_depth); profile_depth_peak(i) = xx(ii);
h1(i) = plot(xx, profile_depth, 'color',COLOR(i,:))
set(gca,'XLim',[2*cfg.unitinmm 1])
axis square
xlabel('distance from surface (mm)')


Jacobian_2d = squeeze(Jacobian{i}(:,150,:)); % take the middle slice
figure(fig5), hold on, 
% h(i) = plot([1:300].*cfg.unitinmm, squeeze(sum(Jacobian_2d,1)), 'color',COLOR(i,:))
xx = [-149:150].*cfg.unitinmm;
profile_hor = squeeze(sum(Jacobian_2d(:,bounds1(1):end),2))./max( squeeze(sum(Jacobian_2d(:,bounds1(1):end,:),2)) );
profile_hor_dia(i) = 2*interp1(profile_hor(150:end),xx(150:end),1/exp(1));
h2(i) = plot(xx(51:250), profile_hor(51:250), 'color',COLOR(i,:))
set(gca,'XLim',[-1 1])
xlabel('distance from the center (mm)')

end


s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,3); s2(:) = {'nm, peak at '};
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeastoutside')

s2 = cell(1,3); s2(:) = {'nm, diameter = '};
s3 = arrayfun(@num2str,profile_hor_dia,'UniformOutput',false);
legend(h2, strcat(s1,s2,s3),'location','northeastoutside')
%% perturbation simulation
f2 = cell(1,length(ind));
det2 = f2;
traj2 = f2;
NewDetWeight_depth = f2;
NewDetWeight_hor = f2;
NewDetWeight_r = f2;
DetWeight = zeros(1,3);
profile_depth = f2;
profile_hor = f2;
profile_r = f2;
SimVol = 200; % "replay" simulation volume, does not have to be the entire volume (to save time)
%% get initial total photon weights
for i = 1:3

% change properties
switch cfg.config
case 'without skull'
    if isnan(cfg.perturb)
    cfg.prop = [0 0 1 1            % medium 0: the environment
    brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
    brain_mua(ind(i))*1.1 gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37 ];  % medium 3: gray matter, perturbed
    else
    cfg.prop = [0 0 1 1            % medium 0: the environment
    brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
    brain_mua(ind(i))+cfg.perturb gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37 ];  % medium 3: gray matter, perturbed
    end           
case 'with skull'
otherwise 
    disp('Check property setting (cfg.prop)')
end

det_temp = det1{i};
det_temp.ppath = det_temp.ppath(det_temp.DetInd,:); % select photons by its exit angle
DetWeight(i) = sum(mcxdetweight(det_temp,cfg.prop));
           
% calculate total photon weights after perturbation   
newcfg              = cfg;
newcfg.seed         = seeds{i}.data(:,det1{i}.DetInd);
newcfg.detphotons   = det1{i}.data(:,det1{i}.DetInd);
newcfg.outputtype   = 'jacobian';

%  ========== change properties in Z direction  ==========
% NewDetWeight_depth{i} = zeros(1,SimVol);
% fwait1 = waitbar(0,'Vertical perturbation simulation started ...');
% for zz = 1:SimVol
%     waitbar(zz/SimVol,fwait1,['Vertical perturbation simulation processing ',num2str(zz),'/',num2str(SimVol)]);
% 
%     % ========== define tissue structure ==========
%     newcfg.vol = cfg.vol;
%     switch cfg.config
%         case 'without skull'
%             newcfg.vol(:,:,zz) = 2; % gray matter, perturbed
%             newcfg.vol=uint8(newcfg.vol);
%         case 'with skull'
%         otherwise
%     end
%     tic, [f2{i}, det2{i}, vol2, seeds2, traj2{i}] = mcxlab(newcfg); toc
%     det_temp = det2{i};
%     NewDetWeight_depth{i}(zz) = sum(mcxdetweight(det_temp,newcfg.prop));
%     
% end
% profile_depth{i} = -(NewDetWeight_depth{i} - DetWeight(i)).*100./DetWeight(i);

% ========== change properties in X direction  ==========
% NewDetWeight_hor{i} = zeros(1,SimVol);
% fwait2 = waitbar(0,'Horizontal perturbation simulation started ...');
% 
% for xx = 150-SimVol/2+1: 150+SimVol/2
%     waitbar((xx-(150-SimVol/2))/SimVol,fwait2,['Horizontal perturbation simulation processing ',num2str((xx-(150-SimVol/2))),'/',num2str(SimVol)]);
%     % ========== define tissue structure ==========
%     newcfg.vol = cfg.vol;
%     switch cfg.config
%         case 'without skull'
%             newcfg.vol(xx,150,:) = 2; % gray matter, perturbed
%             newcfg.vol=uint8(newcfg.vol);
%         case 'with skull'
%         otherwise
%     end
%     tic, [f2{i}, det2{i}, vol2, seeds2, traj2{i}] = mcxlab(newcfg); toc
%     det_temp = det2{i};
%     NewDetWeight_hor{i}(xx-(150-SimVol/2)) = sum(mcxdetweight(det_temp,newcfg.prop));
% end
% profile_hor{i} = -(NewDetWeight_hor{i} - DetWeight(i)).*100./DetWeight(i);


% ========== change properties in R direction  ==========
NewDetWeight_r{i} = zeros(1,SimVol/2);
fwait3 = waitbar(0,'Horizontal perturbation simulation started ...');
theta = 0:pi/50:2*pi;
profile_r{i} = zeros(1,SimVol/2+1);
for rr = 0 : SimVol/2
    waitbar(rr/(SimVol/2),fwait3,['Horizontal perturbation simulation processing ',num2str(rr),'/',num2str(SimVol/2)]);
    % ========== define tissue structure ==========
    newcfg.vol = cfg.vol;
    switch cfg.config
        case 'without skull'
            xx = rr.*cos(theta);
            yy = rr.*sin(theta);
            newcfg.vol(150+round(xx),150+round(yy),:) = 2; % gray matter, perturbed
            newcfg.vol=uint8(newcfg.vol);
        case 'with skull'
        otherwise
    end
%     plot(150+round(xx),150+round(yy)), hold on, pause
    VolDiff = newcfg.vol(:,:,1) - cfg.vol(:,:,1);
    nChange(rr+1) = length(find(VolDiff(:)));
    tic, [f2{i}, det2{i}, vol2, seeds2, traj2{i}] = mcxlab(newcfg); toc
    det_temp = det2{i};
    NewDetWeight_r{i}(rr+1) = sum(mcxdetweight(det_temp,newcfg.prop));
end
profile_r{i} = -(NewDetWeight_r{i} - DetWeight(i))./nChange.*100./DetWeight(i);
profile_r{i} = [fliplr(profile_r{i}(2:end)) profile_r{i}];
end

%% plot along z axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
z = [1:SimVol].*cfg.unitinmm;
profile_depth_peak = zeros(1,3);

for i = 1
    [maxy,ii] = max(profile_depth{i}); profile_depth_peak(i) = z(ii);
    % h1(i) = plot(z, profile_depth{i}, 'color',COLOR(i,:)); hold on
    h1(i) = plot(z, profile_depth{i}, 'color',COLOR(i,:)); hold on
    set(gca,'XLim',[0 4], 'YLim',[0 1])
    xlabel('distance from surface (mm)')
    ylabel('\DeltaSignal (%) when \Delta\mu_a=10%')
%     xlim([0 1])
    axis square

end

s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,3); s2(:) = {'nm, peak at '};
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeastoutside')

%% plot along x axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
x = [-SimVol/2+1:SimVol/2].*cfg.unitinmm;
profile_hor_dia = zeros(1,3);
for i = 1:3
    profile = profile_hor{i};
%     h2(i) = plot(x, profile./brain_mua(ind(i))/0.1, 'color',COLOR(i,:)); hold on
    h2(i) = plot(x, profile./max(profile), 'color',COLOR(i,:)); hold on

%     [~, ii] = min(abs((profile./max(profile) - 1/exp(1))));
    [~, ii] = min(abs((profile./max(profile) - 1/2)));
    profile_hor_dia(i) = 2*abs(x(ii));

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
    xlim([-1 1])
    axis square
end
s2 = cell(1,3); s2(:) = {'Straight illumination, \newline FWHM = '};
s3 = arrayfun(@num2str,profile_hor_dia ,'UniformOutput',false);
legend(h2, strcat(s2,s3), 'location', 'northeastoutside')
if isnan(cfg.perturb)
%     title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change 10%'],'fontsize',14);
    title([cfg.config, '\newline absorption coeff. change 1(/mm)'],'fontsize',14);
else
    title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change ', num2str(cfg.perturb)],'fontsize',14);
end

%% plot along R axis
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
x = [-SimVol/2:SimVol/2].*cfg.unitinmm;
profile_hor_dia = zeros(1,3);
for i = 1:3
    profile = profile_r{i};
%     h2(i) = plot(x, profile./brain_mua(ind(i))/0.1, 'color',COLOR(i,:)); hold on
    h2(i) = plot(x, profile./max(profile), 'color',COLOR(i,:)); hold on
    [~, ii] = min(abs((profile./max(profile) - 1/2)));
    profile_r_dia(i) = 2*abs(x(ii));

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
    xlim([-1 1])
    axis square
end
s2 = cell(1,3); s2(:) = {'Straight illumination, \newline FWHM = '};
s3 = arrayfun(@num2str,profile_r_dia ,'UniformOutput',false);
legend(h2, strcat(s2,s3), 'location', 'northeastoutside')

% title([cfg.config, ', wavelength ',num2str(WLs(i)),'nm \newline absorption coeff. change 10%'],'fontsize',14);
% title([cfg.config, '\newline absorption coeff. change 1(/mm)'],'fontsize',14);
title([cfg.config, ', normalized'],'fontsize',14);
