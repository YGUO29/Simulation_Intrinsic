%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
%
% simulate a 4-layer brain model using MCXLAB.
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cfg;
%% set tissue parameters
load('C:\Users\guoyu\Documents\MATLAB\Simulation_Intrinsic\optical_para.mat')
load('C:\Users\guoyu\Documents\MATLAB\Simulation_Intrinsic\spectra.mat')
lambda = spectra(:,1);
for i = 1:length(optical_para_label)
    eval([optical_para_label{i}, '= optical_para(:,i);']);
end

WLs = [470 530 590 625 730 850];
ind = floor(interp1(lambda,1:length(lambda), WLs));
%% preparing the input data

    % set seed to make the simulation repeatible
    cfg.seed=hex2dec('623F9A9E'); 

    cfg.nphoton=5e7;
    cfg.unitinmm = 0.05; % define units, 10 um
    % time-domain simulation parameters
    cfg.tstart=0;
    cfg.tend=5e-9;
    cfg.tstep=5e-10;
    
% ========== define source and detector ==========
%     cfg.srctype = 'pencil';
    NA = 0.05;
    cfg.srctype = 'cone';
    cfg.srcparam1=[asin(NA) 0 0 0];
    % cfg.srcparam1 = 100; %radius for the disk
    cfg.size = [300 300 300]; % 15mm cube
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/2, cfg.size(2)/2)]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
    cfg.maxjumpdebug = cfg.nphoton;
    cfg.maxexitangle = 0.05;
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

    cfg.isreflect=0; % disable reflection at exterior boundary
    % GPU thread configuration
    cfg.autopilot=1;
    cfg.gpuid=1;
% ========== output control ==========
    cfg.issaveexit = 1;
%%
f1 = cell(1,length(ind));
seeds = f1;
det1 = f1;
traj = f1;
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')
for i = 1:length(ind)
% for i = 1
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
    [f1{i},det1{i},vol,seeds{i},traj{i}]=mcxlab(cfg);
    toc;
% ========== detected photons ==========
    ExitAngle   = sin(vecnorm(det1{i}.v(:,1:2)')./vecnorm(det1{i}.v')); % sin_theta (NA)
    det1{i}.DetInd      = find(ExitAngle < NA);
    det1{i}.DetNumber   = length(det1{i}.DetInd);
    DetExitPosition = det1{i}.p(det1{i}.DetInd,1:2);
% ========== plot figure ==========
    subplot(2,3,i);
    I1 = log10(squeeze(sum(f1{i}.data(:,size(cfg.vol,1)./2+1,:,:),4))');
    I1(I1<0) = 0;
        if i == 1
        cmax1 = max(I1(:)); cmin1 = min(I1(:));
        end
    contourf(I1, cmin1:0.5:cmax1);
    hold on
    for k = 1:length(bounds1)
        plot([0 size(cfg.vol,1)],[bounds1(k) bounds1(k)],'--');
    end
    title(['flux with no reflection, no skull, wavelength ',num2str(WLs(i)),'nm'],'fontsize',14);
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

%%
f2 = cell(1,length(ind));
det2 = cell(1,length(ind));
jac = cell(1,length(ind));
cc = jac;
fig1 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
figure(fig1), hold on
% for k = 1:length(bounds1)
%     plot([bounds1(k) bounds1(k)].*cfg.unitinmm,get(gca,'YLim'),'--');
% end

newcfg=cfg;
newcfg.outputtype='jacobian';

% for i = 1
for i = 1:length(ind)
newcfg.seed=seeds{i}.data(:,det1{i}.DetInd);
newcfg.detphotons=det1{i}.data(:,det1{i}.DetInd);

cfg.prop=[0 0 1 1            % medium 0: the environment
skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 3: white matter

[f2{i}, det2{i}, vol2, seeds2]=mcxlab(newcfg);
jac{i}=sum(f2{i}.data,4);

cc{i} = log10(abs(squeeze(sum(jac{i}(:,size(cfg.vol,1)./2+1,:,:),4))));
cc{i}(isinf(cc{i})) = nan;
cc{i} = cc{i}';
figure(fig1),
if i == 1
    cmin = min(cc{i}(:));
    cmax = max(cc{i}(:));
end
subplot(2,3,i), contourf(cc{i},cmin:0.5:cmax)
hold on
for k = 1:length(bounds1)
    plot([0 size(cfg.vol,1)],[bounds1(k) bounds1(k)],'--');
end
title(['',num2str(WLs(i)),'nm'],'fontsize',14);
set(gca,'clim',[cmin cmax]);
colorbar
drawnow
end
%%
fig2 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
figure(fig2), hold on
for k = 1:length(bounds1)
plot([bounds1(k) bounds1(k)].*cfg.unitinmm,[-5 -3],'--k');
% plot([bounds1(k) bounds1(k)].*cfg.unitinmm,[0 1],'--k');
end

for i = 1:6
hold on,
% data_temp = mean(cc{i},2,'omitnan');
data_temp = cc{i}(:,151);

% data_temp = data_temp./max(data_temp);
h(i) = plot((1:cfg.size).*cfg.unitinmm, data_temp);
set(gca,'XLim',[0 2.5],'YLim',[-5 -3])
% if i == 1
%     for k = 1:length(bounds1)
%     plot([bounds1(k) bounds1(k)].*cfg.unitinmm,get(gca,'YLim'),'--');
%     end
% end
end
legend(h,arrayfun(@num2str,WLs,'UniformOutput',false),'location','northeastoutside')
%%
% fig3 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
% figure(fig3), hold on
% for k = 1:length(bounds1)
% plot([bounds1(k) bounds1(k)].*cfg.unitinmm,[-5.6 -5],'--k');
% end
% for i = 1:6
% hold on,
% h(i) = plot((1:cfg.size).*cfg.unitinmm, mean(cc{i},2,'omitnan') );
% set(gca,'XLim',[0 2.5],'YLim',[-5.6 -5])
% % if i == 1
% %     for k = 1:length(bounds1)
% %     plot([bounds1(k) bounds1(k)].*cfg.unitinmm,get(gca,'YLim'),'--');
% %     end
% % end
% end
% legend(h,arrayfun(@num2str,WLs,'UniformOutput',false),'location','northeastoutside')
