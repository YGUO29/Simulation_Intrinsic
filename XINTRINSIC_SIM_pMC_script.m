%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temp: convert data format from XS to YG
wavelength = spectra(:,1);
skull_g = 0.9;
skull_mus = interp1(D.WaveLength2, D.MuSctPSkull, wavelength)./(1-skull_g);
skull_mua = interp1(D.MuAbsSkull(:,1),D.MuAbsSkull(:,2),wavelength,'linear','extrap');
brain_mua = D.MuAbsHbT;
gray_matter_mus = interp1(D.WaveLength2, D.MuSctBrainGM, wavelength);
white_matter_mus = interp1(D.WaveLength2, D.MuSctBrainWM, wavelength);

optical_para_xs = optical_para;
optical_para_xs(:,1) = skull_mus;
optical_para_xs(:,2) = skull_mua;
optical_para_xs(:,7) = brain_mua;

%% load updated parameters
addpath(genpath(cd));
clear cfg;
% set tissue parameters
load('optical_para_xs.mat')
load('spectra.mat')
lambda = spectra(:,1);
% figure,
for i = 1:length(optical_para_label)
    eval(['para.',optical_para_label{i}, '= optical_para_xs(:,i);']);
%     subplot(1,7,i), plot(eval(['para.',optical_para_label{i}])), hold on
end
para.skull_g = 0.92;
para.skull_thick = 0.5; 
para.gm_thick = 1.3;
para.WLs = [470 530 590 625 730 850];
para.ind = floor(interp1(lambda,1:length(lambda), para.WLs));
% COLOR = hex2rgb(['1B97F7'; '23C249'; 'FAB738'; 'E92716'; '874040'; '000000']);
para.COLOR = [   0       0       1;
            0       0.5     0;
            1       0.5     0;
            1       0       0;
            0.5     0       0;
            0.25    0       0];
        
%% 
% [cfg, det1, seeds] = MC_test(config, [0 1], para);
% newcfg              = cfg;
% newcfg.seed         = seeds{i}.data(:,det1{i}.DetInd);
% newcfg.detphotons   = det1{i}.data(:,det1{i}.DetInd);
% [S.f2{i}, S.det2{i}, S.vol2, S.seeds2, S.traj2{i}] = mcxlab(newcfg);
profile_x1 = cell(1,6);
tEnd = zeros(1, 50);
for i = 1
    tStart = tic;
    para.iRep = i;
    exitangle = [0 0.03];
    config = 'with skull';
    para.profile = profile_x1;
%     para.profile = [];
    para.SimVol = 122;
%     [cfg, det1, seeds] = MC(config, exitangle, para);
    [cfg, det1, seeds] = MC_test0806(config, exitangle, para);
    [profile_x1, x1, S1]  = pMC(cfg, det1, seeds, 'x', para);
    close all
    tEnd(i) = toc(tStart);
end

%%
% [0.8548 0.8844] angled with 0.87NA, deviation = asin(0.03) (�1.7191
% degrees), 0.87 determined by objective radius = 16.32, working distance =
% 9.25
% [0.5757 0.6237] angled with 0.6NA, deviation = asin(0.03) (�1.7191 degrees)
exitangle = [0 0.03];
config = 'without skull';
[cfg, det1, seeds] = MC1(config, exitangle, para);
[profile_depth1, z1, S1]  = pMC(cfg, det1, seeds, 'z', para);

%
exitangle = [0.8548 0.8844];
config = 'without skull';
[cfg, det1, seeds] = MC1(config, exitangle, para);
[profile_depth2, z2, S2]  = pMC(cfg, det1, seeds, 'z', para);

%
exitangle = [0 0.03];
config = 'without skull off focus';
[cfg, det1, seeds] = MC1(config, exitangle, para);
[profile_depth3, z3, S3]  = pMC(cfg, det1, seeds, 'z', para);

%
exitangle = [0.8548 0.8844];
config = 'without skull off focus';
[cfg, det1, seeds] = MC1(config, exitangle, para);
[profile_depth4, z4, S4]  = pMC(cfg, det1, seeds, 'z', para);

%
exitangle = [0.8548 0.8844];
config = 'with skull';
[cfg, det1, seeds] = MC(config, exitangle, para);
[profile_depth6, z6, S6]  = pMC(cfg, det1, seeds, 'z', para);

%%
exitangle = [0 0.03];
config = 'with skull off focus';
[cfg, det1, seeds] = MC(config, exitangle, para);
[profile_x7, x7, S7]  = pMC(cfg, det1, seeds, 'x', para);
%
exitangle = [0.8548 0.8844];
config = 'with skull off focus';
[cfg, det1, seeds] = MC(config, exitangle, para);
[profile_x8, x8, S8]  = pMC(cfg, det1, seeds, 'x', para);
%% plot along z axis

nProfile = 2;
profile = cell(1,nProfile); x = profile;
profile{1} = profile_depth3{2};
x{1} = z3;
profile{2} = profile_depth4{2};
x{2} = z4;

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);

for i = 1:nProfile
    
if contains(cfg.config,'without skull')
    delta_mua = 0.1*para.brain_mua(para.ind(2)).*ones(1,length(profile{i}));
%     .*(cfg.unitinmm^3*cfg.size(1)*cfg.size(2))
    profile{i} = profile{i}./delta_mua;
else
    delta_mua = [0.1*para.skull_mua(para.ind(2)).*ones(1,cfg.bounds(1)), 0.1*para.brain_mua(para.ind(i)).*ones(1,length(profile{i})-cfg.bounds(1))];
    profile{i} = profile{i}./delta_mua;
end

[maxy,ii] = max(profile{i}); profile_depth_peak(i) = x{i}(ii);
% h1(i) = plot(z, profile_depth{i}, 'color',COLOR(i,:)); hold on
% h1(i) = plot(x, profile{i}, 'color',COLOR(2,:)); hold on
h1(i) = plot(x{i}, profile{i}); hold on

% set(gca,'yscale','log')
set(gca,'XLim',[0 max(x{i})], 'YLim',[min(profile{i}) max(profile{i})])

for k = 1:length(cfg.bounds)
plot([cfg.bounds(k)+1 cfg.bounds(k)+1].*cfg.unitinmm,get(gca,'YLim'),'LineStyle',':','color',[0.5 0.5 0.5]);
end

xlabel('distance from surface (mm)')
ylabel('\DeltaSignal (%) / \Delta\mu_a')

end
s1 = arrayfun(@num2str,para.WLs(2),'UniformOutput',false); s1 = repmat(s1,1,nProfile);
s2 = [{'nm,straight illu., peak at '}, {'nm, angled illu., peak at '}];
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeast')
title(['Depth sensitivity, ', cfg.config])

%% plot along R axis
cfg.config = 'without skull off focus';

nProfile = 2;
profile = cell(1,nProfile); x = profile;
profile{1} = profile_r3{2};
x{1} = x3.*cfg.unitinmm;
profile{2} = profile_r4{2};
x{2} = x4.*cfg.unitinmm;

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
profile_r_dia = zeros(1,nProfile);
for i = 1:nProfile
    delta_mua = 0.1*para.brain_mua(para.ind(2)).*ones(1,length(profile{i}));
    profile{i} = profile{i}./delta_mua;
    h1(i) = plot(x{i}, profile{i}); hold on

    [~, ii] = min(abs((profile{i} - max(profile{i})./2)));
    profile_r_dia(i) = 2*abs(x{i}(ii));

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
    xlim([-2.5 2.5])
end
s1 = arrayfun(@num2str,para.WLs(2),'UniformOutput',false); s1 = repmat(s1,1,nProfile);
s2 = [{'nm,straight illu., diameter = '}, {'nm, angled illu., diameter = '}];
s3 = arrayfun(@num2str,profile_r_dia,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeast')
title(['Spatial Resolution (along R direction) ', cfg.config])
%% plot along X axis

nProfile = 2;
profile = cell(1,nProfile); x = profile;
profile{1} = profile_x1{2};
x{1} = x1;
profile{2} = profile_x4{2};
x{2} = x4;

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
profile_r_dia = zeros(1,nProfile);
for i = 1:nProfile
    delta_mua = 0.1*para.brain_mua(para.ind(2)).*ones(1,length(profile{i}));
    profile{i} = profile{i}./delta_mua;
    h1(i) = plot(x{i}, profile{i}); hold on

    [~, ii] = min(abs((profile{i} - max(profile{i})./2)));
    profile_r_dia(i) = 2*abs(x{i}(ii));

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
    xlim([-2.5 2.5])
end
s1 = arrayfun(@num2str,para.WLs(2),'UniformOutput',false); s1 = repmat(s1,1,nProfile);
s2 = [{'nm,straight illu., diameter = '}, {'nm, angled illu., diameter = '}];
s3 = arrayfun(@num2str,profile_r_dia,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeast')
title(['Spatial Resolution (along X direction) ', cfg.config])

%% plot along X axis, multiple wavelengths
nProfile = 6;
% profile = profile_x1;
% x = x1;


figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440         918        1013         420]);
profile_r_dia = zeros(1,nProfile);
for i = 1:nProfile
    delta_mua   = 0.1*para.brain_mua(para.ind(2)).*ones(1,length(profile{i}));
    profile{i}  = mean(profile{i},1)./delta_mua;
%     profile{i}  = profile{i}(3,:)./delta_mua;
%     h1(i)       = plot(x, profile{i}, 'color', para.COLOR(i,:)); hold on
    h1(i)       = plot(x, profile{i}./max(profile{i}), 'color', para.COLOR(i,:)); hold on

    [~, ii] = min(abs((profile{i} - max(profile{i})./2)));
    profile_r_dia(i) = 2*abs(x(ii));

    xlabel('distance from the center (mm)')
    ylabel('reflectance change (%)')
    xlim([-2.5 2.5])
end
s1 = arrayfun(@num2str,para.WLs,'UniformOutput',false);
s2 = {'nm, diameter = '}; s2 = repmat(s2,1,nProfile);
s3 = arrayfun(@num2str,profile_r_dia,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeast')
title(['Spatial Resolution (along X direction) ', cfg.config])
%% plot example trajectories
traj = traj2{2};
AxLim = [1 400 1 500];
XINTRINSIC_SIM_PlotTraj(traj, cfg.bounds, AxLim)