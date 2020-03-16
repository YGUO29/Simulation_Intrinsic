%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
%
% simulate a 4-layer brain model using MCXLAB.
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cfg;
%% set tissue parameters
load('D:\=code=\Simulation_Intrinsic\optical_para.mat')
load('D:\=code=\Simulation_Intrinsic\spectra.mat')
lambda = spectra(:,1);
for i = 1:length(optical_para_label)
    eval([optical_para_label{i}, '= optical_para(:,i);']);
end

WLs = [470 530 590 625 730 850];
ind = floor(interp1(lambda,1:length(lambda), WLs));
% COLOR = hex2rgb(['1B97F7'; '23C249'; 'FAB738'; 'E92716'; '874040'; '000000']);
COLOR = [   0       0       1;
            0       0.5     0;
            1       0.5     0;
            1       0       0;
            0.5     0       0;
            0.25    0       0];
%% preparing the input data
    % set seed to make the simulation repeatible
    cfg.seed=hex2dec('623F9A9E'); 

    cfg.nphoton=1e7;
    cfg.unitinmm = 0.05; % define units, 10 um
    % time-domain simulation parameters
    cfg.tstart=0;
    cfg.tend=5e-9;
    cfg.tstep=5e-9;
    cfg.issrcfrom0 = 1;

% ========== define source and detector ==========
%     cfg.srctype = 'pencil';
    NA = 0.03;
    cfg.srctype = 'cone';
    cfg.srcparam1=[asin(NA) 0 0 0];
    % cfg.srcparam1 = 100; %radius for the disk
    cfg.size = [300 300 300]; % 15mm cube
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/2, cfg.size(2)/2)]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
    cfg.maxexitangle = NA;
    cfg.outputtype = 'fluence';

% ========== define a 3 layer structure ==========
cfg.config = 'without skull';
switch cfg.config
    case 'without skull'
        cfg.vol = ones(cfg.size); % 15 mm volume
        bounds1 = floor(1.327/cfg.unitinmm); % gray matter,
        cfg.vol(:,:,bounds1+1:end)=2; % the rest: white matter
        cfg.vol=uint8(cfg.vol);
    case 'with skull'
        cfg.vol = ones(cfg.size); % 15 mm volume
        bounds1 = [floor(0.7/cfg.unitinmm) floor(0.7/cfg.unitinmm)+floor(1.327/cfg.unitinmm)]; % skull, gray matter, white matter
        cfg.vol(:,:,bounds1(1)+1:bounds1(2))=2; % 0.7mm: skull
        cfg.vol(:,:,bounds1(2)+1:end)=3; % the rest: white matter
        cfg.vol=uint8(cfg.vol);
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
seeds1 = f1;
det1 = f1;
traj1 = f1;
Phi_s = f1;
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
% for i = 1:length(ind)
for i = 2
    switch cfg.config
        case 'without skull'
            cfg.prop=[0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 2: white matter
        case 'with skull'
            cfg.prop=[0 0 1 1            % medium 0: the environment
            skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 3: white matter
        otherwise 
            disp('Check property setting (cfg.prop)')
    end
   
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    % f1=mcxlab(cfg);
    [f1{i},det1{i},vol,seeds1{i},traj1{i}]=mcxlab(cfg);
    toc;
% ========== detected photons ==========
    ExitAngle   = vecnorm(det1{i}.v(:,1:2)'); % sin_theta (i.e. NA), v is a normal vector
    det1{i}.DetInd      = find(ExitAngle < NA);
%     det1{i}.DetInd      = find(ExitAngle < 0.6237 & ExitAngle > 0.5757); % take sine(angle) = 0.6, within a ±1.7191 degrees range
    det1{i}.DetNumber   = length(det1{i}.DetInd);
    det1{i}.DetExitPos  = det1{i}.p(det1{i}.DetInd,1:2);
% ========== plot figure ==========
    subplot(2,3,i);
%     I1 = log10(squeeze(sum(f1{i}.data(:,size(cfg.vol,1)./2+1,:,:),4))');
    I1 = squeeze(sum(f1{i}.data,4));
    I1 = log10(squeeze( sum(I1,2))');
    I1(I1<0) = 0;
        if i == 1
        cmax1 = max(I1(:)); cmin1 = min(I1(:));
        end
    contourf(I1, cmin1:0.25:cmax1);
%     imagesc(flipud(I1),[cmin, cmax])
    hold on
    for k = 1:length(bounds1)
        plot([0 size(cfg.vol,1)],[bounds1(k) bounds1(k)],'--');
    end
    title(['fluence with reflection, wavelength ',num2str(WLs(i)),'nm'],'fontsize',14);
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
% ========== switch source and detector, ring illumination ==========
% r1      = 200;
% r2      = 190;
% psize   = 400;
% pat     = zeros(psize,psize);
% [xi,yi] = meshgrid(0:psize-1,0:psize-1);
% dx      = (xi-(psize-1)/2);
% dy      = (yi-(psize-1)/2);
% pat(dx.*dx+dy.*dy<=r1*r1 & dx.*dx+dy.*dy>=r2*r2) = 1;
% cfg.issrcfrom0 = 1;
% cfg.srcpos = [0 0 -100];
% cfg.srctype = 'pattern';
% cfg.srcpattern = pat;
% cfg.srcparam1 = [cfg.size(1) 0 0 psize];
% cfg.srcparam2=[0 cfg.size(2) 0 psize];
% cfg.srcdir=[0 0 1 cfg.size(1)/2*tan(asin(0.8))];



figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
% for i = 1:length(ind)
for i = 2
    switch cfg.config
        case 'without skull'
            cfg.prop=[0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 2: white matter
        case 'with skull'
            cfg.prop=[0 0 1 1            % medium 0: the environment
            skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 3: white matter
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
%         if i == 1
        cmax1 = max(I1(:)); cmin1 = min(I1(:));
%         end
    contourf(I1, cmin1:0.1:cmax1);
%     imagesc(flipud(I1),[cmin, cmax])
    hold on
    for k = 1:length(bounds1)
        plot([0 size(cfg.vol,1)],[bounds1(k) bounds1(k)],'--');
    end
    title(['fluence with reflection, wavelength ',num2str(WLs(i)),'nm'],'fontsize',14);
    set(gca,'clim',[cmin1 cmax1]);
    colorbar
    drawnow
    
    Phi_d{i} = squeeze(sum(f2{i}.data,4));
end
%% Jacobian
fig3 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
fig4 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);
fig5 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w','Position',[1440 918 798 420]);



Jacobian = cell(1,6);
profile_depth_peak = zeros(1,6);
profile_hor_dia = zeros(1,6);
for i = 2
    
det_temp = det1{i};
det_temp.ppath = det_temp.ppath(det_temp.DetInd,:); % select photons by its exit angle
DetWeight(i) = sum(mcxdetweight(det_temp,cfg.prop));

Jacobian{i} = cfg.unitinmm.^3.*Phi_s{i}.*Phi_d{i}./ (DetWeight(i)/cfg.nphoton);
Jacobian_2d = squeeze( sum(Jacobian{i},2) );

figure(fig3), subplot(2,3,i), 
Jacobian_temp = Jacobian_2d(:, 1:80);
cmax1 = max(Jacobian_temp(:)); cmin1 = min(Jacobian_temp(:));
% imagesc(log(flipud(Jacobian_temp')), [-40 -8])
contourf(log(Jacobian_temp'), -40:1:0)
hold on
    for k = 1:length(bounds1)
        plot([0 size(cfg.vol,1)],[bounds1(k) bounds1(k)],'--');
    end
axis image
title(['Jacobian', num2str(WLs(i)),' nm']), colorbar


figure(fig4), hold on, 
xx = [1:300].*cfg.unitinmm;
% profile_depth = squeeze(sum(Jacobian_2d,1));
profile_depth = squeeze(sum(Jacobian_2d,1))./max(squeeze(sum(Jacobian_2d,1)));
[~,ii] = max(profile_depth); profile_depth_peak(i) = xx(ii);
h1(i) = plot(xx, profile_depth, 'color',COLOR(i,:))
set(gca,'XLim',[2*cfg.unitinmm 4])
if i == 6
    for k = 1:length(bounds1)
    plot([bounds1(k) bounds1(k)].*cfg.unitinmm,get(gca,'YLim'),'LineStyle',':','color',[0.5 0.5 0.5]);
    end
end
xlabel('distance from surface (mm)')

Jacobian_2d = squeeze(Jacobian{i}(:,150,:));
figure(fig5), hold on, 
% h(i) = plot([1:300].*cfg.unitinmm, squeeze(sum(Jacobian_2d,1)), 'color',COLOR(i,:))
xx = [-149:150].*cfg.unitinmm;
profile_hor = squeeze(sum(Jacobian_2d(:,bounds1(1):end),2))./max( squeeze(sum(Jacobian_2d(:,bounds1(1):end,:),2)) );
profile_hor_dia(i) = 2*interp1(profile_hor(150:end),xx(150:end),1/exp(1));
h2(i) = plot(xx(51:250), profile_hor(51:250), 'color',COLOR(i,:))
xlabel('distance from the center (mm)')

end


s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
s2 = cell(1,6); s2(:) = {'nm, peak at '};
s3 = arrayfun(@num2str,profile_depth_peak,'UniformOutput',false);
legend(h1, strcat(s1,s2,s3),'location','northeastoutside')

s2 = cell(1,6); s2(:) = {'nm, diameter = '};
s3 = arrayfun(@num2str,profile_hor_dia,'UniformOutput',false);
legend(h2, strcat(s1,s2,s3),'location','northeastoutside')
