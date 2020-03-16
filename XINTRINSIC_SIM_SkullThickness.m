%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
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
    cfg.tstep=5e-10;
    cfg.issrcfrom0 = 1;

% ========== define source and detector ==========
%     cfg.srctype = 'pencil';
    NA = 0.05;
    cfg.srctype = 'cone';
    cfg.srcparam1=[asin(NA) 0 0 0];
    % cfg.srcparam1 = 100; %radius for the disk
    cfg.size = [300 300 300]; % 15mm cube
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 2];
    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/2, cfg.size(2)/2)]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
    cfg.maxexitangle = NA;
% ========== define a 3 layer structure ==========
cfg.config = 'with skull';
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
%     cfg.debuglevel = 'M';
    cfg.maxjumpdebug = 5e7; % number of trajectory positions saved
%%    
SkullThick = 0.7:0.3:2.2;
i = 2; % ith wavelength


f1 = cell(1,length(SkullThick));
seeds = f1;
det1 = f1;
traj = f1;
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')
for j = 1:length(SkullThick)
    switch cfg.config
        case 'without skull'
            cfg.prop=[0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 2: white matter
        
            cfg.vol = ones(cfg.size); % 15 mm volume
            bounds1 = floor(1.327/cfg.unitinmm); % gray matter,
            cfg.vol(:,:,bounds1+1:end)=2; % the rest: white matter
            cfg.vol=uint8(cfg.vol);
        case 'with skull'
            cfg.prop=[0 0 1 1            % medium 0: the environment
            skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 3: white matter
        
            cfg.vol = ones(cfg.size); % 15 mm volume
            bounds1 = [floor(SkullThick(j)/cfg.unitinmm) floor(SkullThick(j)/cfg.unitinmm)+floor(1.327/cfg.unitinmm)]; % skull, gray matter, white matter
            cfg.vol(:,:,bounds1(1)+1:bounds1(2))=2; % 0.7mm: skull
            cfg.vol(:,:,bounds1(2)+1:end)=3; % the rest: white matter
            cfg.vol=uint8(cfg.vol);
        otherwise 
            disp('Check property setting (cfg.prop)')
    end
   
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    % f1=mcxlab(cfg);
    [f1{j},det1{j},vol,seeds{j}]=mcxlab(cfg);
    toc;
% ========== detected photons ==========
    ExitAngle   = sin(vecnorm(det1{j}.v(:,1:2)')./vecnorm(det1{j}.v')); % sin_theta (NA)
    det1{j}.DetInd      = find(ExitAngle < NA);
    det1{j}.DetNumber   = length(det1{j}.DetInd);
    det1{j}.DetExitPos  = det1{j}.p(det1{j}.DetInd,1:2);
% ========== plot figure ==========
    subplot(2,3,j);
    I1 = log10(squeeze(sum(f1{j}.data(:,size(cfg.vol,1)./2+1,:,:),4))');
    I1(I1<0) = 0;
        if j == 1
        cmax1 = max(I1(:)); cmin1 = min(I1(:));
        end
    contourf(I1, cmin1:0.5:cmax1);
    axis square
    hold on
    for k = 1:length(bounds1)
        plot([0 size(cfg.vol,1)],[bounds1(k) bounds1(k)],'--');
    end
    title(['flux with reflection, no skull, wavelength ',num2str(WLs(i)),'nm'],'fontsize',14);
    set(gca,'clim',[cmin1 cmax1]);
    colorbar
    drawnow
    
%     subplot(3,6,6+i);
%     h1 = histogram2(DetExitPosition(:,1),DetExitPosition(:,2),100,...
%         'XBinLimits',[cfg.detpos(1)-10, cfg.detpos(1)+10],'YBinLimits',[cfg.detpos(2)-10, cfg.detpos(2)+10],...
%         'DisplayStyle','tile','ShowEmptyBins','on');
%     DetExitProfile = sum(h1.Values,2);
%         [mm, ii] = max(DetExitProfile); cutoff = mm./2;
%         [~,xx1] = min(abs(DetExitProfile(1:ii)-cutoff));
%         [~,xx2] = min(abs(DetExitProfile(ii:end)-cutoff));    
%     DetExitSize = (xx2 + ii - xx1).*cfg.unitinmm;
%     title(['Exit position, FWHM: ',num2str(DetExitSize,3),'mm'])
%     
%     subplot(3,6,12+i);
%     DetWeight = mcxdetweight(det1,cfg.prop);
%     DetPathlength = det1.ppath.*repmat(DetWeight(:),1,size(det1.ppath,2));
%     DetMeanPathlength = cfg.unitinmm*sum(DetPathlength(:,2)) / sum(DetWeight(:));
%     h2 = histogram(DetPathlength);
%     title(['Pathlength distribution, mean: ',num2str(DetMeanPathlength,3),'mm'])
end
    
%% calculate signal from skull & from cortex
DetWeight_skull = zeros(length(det1),length(ind)); % rows: skull thickness; colunms: wavelengths
DetWeight_brain = DetWeight_skull;
for i = 1:length(ind) % wavelength
    cfg.prop=[0 0 1 1            % medium 0: the environment
    skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
    brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
    brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 3: white matter

    for j = 1:length(det1) % skull thickness
        det_temp = det1{j}; det_part = det_temp;
        det_temp.ppath = det_temp.ppath(det_temp.DetInd,:);
%         DetWeight = sum(mcxdetweight(det_temp,cfg.prop));
        pskull = find(det_temp.ppath(:,2) == 0); % photons that only enters the skull
        det_part.ppath = det_temp.ppath(pskull,:);
        DetWeight_skull(j,i) = sum(mcxdetweight(det_part,cfg.prop)); % total weight of detected photons that only entered the skull
        
        pbrain = setdiff(1:size(det_temp.ppath,1),pskull); 
        det_part.ppath = det_temp.ppath(pbrain,:);
        DetWeight_brain(j,i) = sum(mcxdetweight(det_part,cfg.prop));  % total weight of detected photons that entered the brain
        
%         DetPhoton_skull(i,j) = length(pskull)/size(det_temp.ppath,1);
%         DetPhoton_brain(i,j) = length(pbrain)/size(det_temp.ppath,1);

    end
end
%% changing signal amplitude
fig2 = figure('DefaultAxesFontSize',12, 'DefaultLineLineWidth', 2,'color','w');
set(fig2,'position',get(fig2,'position')+[0 -50 50 50]);
clear h1
for k = 1:length(ind)
    subplot(2,3,k), hold on
    
% ===== ratio between skull & brain signal amplitude =====
%     plot(SkullThick, DetWeight_skull(:,k)./(DetWeight_skull(:,k) + DetWeight_brain(:,k)).*100, 'marker','o'), 
%     plot(SkullThick, DetWeight_brain(:,k)./(DetWeight_skull(:,k) + DetWeight_brain(:,k)).*100, 'marker','o')
%     ylabel('Relative amp. (%)')
% ===== ratio of skull & brain signal amplitude / incident amplitude =====
%     plot(SkullThick, DetWeight_skull(:,k)./cfg.nphoton.*100, 'marker','o'), 
%     plot(SkullThick, DetWeight_brain(:,k)./cfg.nphoton.*100, 'marker','o')
%     ylabel('signal amplitude (%)')
% ===== noise and signal ratio =====
%     plot(SkullThick, sqrt(DetWeight_skull(:,k))./cfg.nphoton.*100, 'marker','o'), 
%     plot(SkullThick, DetWeight_brain(:,k)./cfg.nphoton.*100, 'marker','o')
%     ylabel('signal and noise (%)')
% ===== noise and signal ratio =====
    h1(1) = plot(SkullThick, sqrt(DetWeight_skull(:,k))./cfg.nphoton.*100, 'marker','o');
    h1(2) = plot(SkullThick, 0.02.*DetWeight_brain(:,k)./cfg.nphoton.*100, 'marker','o');
    h1(3) = plot(SkullThick, 0.06.*DetWeight_brain(:,k)./cfg.nphoton.*100, 'marker','o');
    h1(4) = plot(SkullThick, 0.1.*DetWeight_brain(:,k)./cfg.nphoton.*100, 'marker','o');
    h1(5) = plot(SkullThick, 0.14.*DetWeight_brain(:,k)./cfg.nphoton.*100, 'marker','o');
    set(gca,'yscale','log')
    ylabel('signal and noise (%)')
    t = title(['wavelength = ',num2str(WLs(k)), ' nm','\newline',' ']);
%     set(t,'position',get(t,'position')+[0 0 15])
    xlabel('skull thickness (mm)')
    if k == 6
        s1 = arrayfun(@num2str,[2 6 10 14],'UniformOutput',false); s1 = [' ',s1];
        s2 = cell(1,4); s2(:) = {'% signal change '}; s2 = ['noise from skull', s2];
        legend(h1, strcat(s1,s2),'location','northeast')
    end
end
%% changing total #photons

fig2 = figure('DefaultAxesFontSize',12, 'DefaultLineLineWidth', 2,'color','w');
set(fig2,'position',get(fig2,'position')+[0 -50 50 50]);
clear h1
for k = 1:length(ind)
    subplot(2,3,k), hold on
% ===== noise and signal ratio, change #photons =====
    h1(1) = plot(SkullThick, sqrt(DetWeight_skull(:,k))./cfg.nphoton.*100.*sqrt(2), 'marker','o'); % 0.5x N_Photons
    h1(2) = plot(SkullThick, sqrt(DetWeight_skull(:,k))./cfg.nphoton.*100, 'marker','o'); % 1x N_Photons
    h1(3) = plot(SkullThick, sqrt(DetWeight_skull(:,k))./cfg.nphoton.*100./sqrt(2), 'marker','o'); % 2x N_Photons
    h1(4) = plot(SkullThick, sqrt(DetWeight_skull(:,k))./cfg.nphoton.*100./sqrt(3), 'marker','o'); % 3x N_Photons
    h1(5) = plot(SkullThick, 0.1.*DetWeight_brain(:,k)./cfg.nphoton.*100, 'marker','o');
    set(gca,'yscale','log')
    ylabel('signal and noise (%)')
    t = title(['wavelength = ',num2str(WLs(k)), ' nm','\newline',' ']);
%     set(t,'position',get(t,'position')+[0 0 15])
    xlabel('skull thickness (mm)')
    if k == 6
        legend(h1, {'#photons = 0.5N', '#photons = N (1e7)', '#photons = 2N', '#photons = 3N', '10% signal change'},'location','northeast')
    end
end

%% plot pathlength
clear DetMeanPathlength
DetMeanPathlength = zeros(6,length(lambda));
for k = 1:6 % rows: different skull thickness
    det2 = det1{k};
    det2.ppath = det1{k}.ppath(det1{k}.DetInd,:);

    for i = 1:length(lambda) % columes: wavelengths
        cfg.prop=[0 0 1 1            % medium 0: the environment
        brain_mua(i) gray_matter_mus(i) gray_matter_g(i), 1.37   % medium 1: gray matter
        brain_mua(i) white_matter_mus(i) white_matter_g(i), 1.37];   % medium 2: white matter

        avgpath = mcxmeanpath(det2,cfg.prop);
        DetMeanPathlength(k,i) = avgpath(2); % pathlength in skull
%         DetMeanPathlength(k,i) = sum(avgpath(2:3)); % pathlength in brain

    end
end

%% plot skull thickness vs L1 (pathlengths in skull)
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w'),
for i = 1:length(ind)
    h1(i) = plot(SkullThick, DetMeanPathlength(:,ind(i)),'marker','o','color',COLOR(i,:)); hold on
end
% set(gca,'yscale','log')
xlabel('skull thickness (mm)')
ylabel('mean pathlength (mm)')
s1 = arrayfun(@num2str,WLs,'UniformOutput',false);
legend(h1, s1,'location','northeastoutside')


