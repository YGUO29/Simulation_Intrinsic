%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    info=mcxlab('gpuinfo')

clear cfg;
% set tissue parameters
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
    cfg.size = 300.*ones(1,3); % 15mm cube
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/sqrt(2), cfg.size(2)/sqrt(2))]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
% ========== define a 3 layer structure ==========
cfg.config = 'with skull';
switch cfg.config
    case 'without skull'
        cfg.vol = ones(cfg.size); % 15 mm volume
        bounds = floor(1.327/cfg.unitinmm); % gray matter,
        cfg.vol(:,:,bounds+1:end)=2; % the rest: white matter
        cfg.vol=uint8(cfg.vol);
    case 'with skull'
        cfg.vol = ones(cfg.size); % 15 mm volume
        bounds = [floor(0.7/cfg.unitinmm) floor(0.7/cfg.unitinmm)+floor(1.327/cfg.unitinmm)]; % skull, gray matter, white matter
        cfg.vol(:,:,bounds(1)+1:bounds(2))=2; % 0.7mm: skull
        cfg.vol(:,:,bounds(2)+1:end)=3; % the rest: white matter
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
seeds = f1;
det1 = f1;
traj = f1;
%%
cfg.perturb = nan; % a number or nan (perturb by 10%)
cfg.exitangle = [0 0.03];
% cfg.exitangle = [0.8548 0.8844]; 
% [0.8548 0.8844] angled with 0.8NA, deviation = asin(0.03) (±1.7191 degrees)
% [0.5757 0.6237] angled with 0.6NA, deviation = asin(0.03) (±1.7191 degrees)
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')


% for i = 1:length(ind)
for i = 2
    switch cfg.config
        case 'without skull'
            if isnan(cfg.perturb)
            cfg.prop = [0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37   % medium 2: white matter
            brain_mua(ind(i))*1.1 gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 3: gray matter, perturbed
            brain_mua(ind(i))*1.1 white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 4: white matter, pertubed
            else
            cfg.prop = [0 0 1 1            % medium 0: the environment
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 1: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37   % medium 2: white matter
            brain_mua(ind(i))+cfg.perturb gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 3: gray matter, perturbed
            brain_mua(ind(i))+cfg.perturb white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 4: white matter, perturbed
            end
        case 'with skull'
            if isnan(cfg.perturb)
            cfg.prop=[0 0 1 1            % medium 0: the environment
            skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 3: white matter
%             skull_mua(ind(i))*1.1 skull_mus(ind(i)) 0.9337, 1.56 % medium 4: skull, perturbed
%             brain_mua(ind(i))*1.1 gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 5: gray matter, perturbed
%             brain_mua(ind(i))*1.1 white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 6: white matter, perturbed
            else
            cfg.prop=[0 0 1 1            % medium 0: the environment
            skull_mua(ind(i)) skull_mus(ind(i)) 0.9337, 1.56 % medium 1: skull
            brain_mua(ind(i)) gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 2: gray matter
            brain_mua(ind(i)) white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37   % medium 3: white matter
            skull_mua(ind(i))+cfg.perturb skull_mus(ind(i)) 0.9337, 1.56 % medium 4: skull, perturbed
            brain_mua(ind(i))+cfg.perturb gray_matter_mus(ind(i)) gray_matter_g(ind(i)), 1.37   % medium 5: gray matter, perturbed
            brain_mua(ind(i))+cfg.perturb white_matter_mus(ind(i)) white_matter_g(ind(i)), 1.37];   % medium 6: white matter, perturbed
            end
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
    det1{i}.DetInd1     = find(det1{i}.p(:,3)<-6e-5);
    det1{i}.DetInd2     = find(ExitAngle < cfg.exitangle(2) & ExitAngle >= cfg.exitangle(1)); 
    det1{i}.DetInd      = intersect(det1{i}.DetInd1, det1{i}.DetInd2);
    det1{i}.DetNumber   = length(det1{i}.DetInd);
    det1{i}.DetExitPos  = det1{i}.p(det1{i}.DetInd,1:2);
% ========== plot figure ==========
    subplot(2,3,i);
    I1 = log10(squeeze(sum(f1{i}.data(:,size(cfg.vol,1)./2+1,:,:),4))');
    I1(I1<0) = 0;
%         if i == 1
        cmax1 = max(I1(:)); cmin1 = min(I1(:));
%         end
    contourf(I1, cmin1:0.5:cmax1);
    axis square
    hold on
    for k = 1:length(bounds)
        plot([0 size(cfg.vol,1)],[bounds(k) bounds(k)],'--');
    end
    title(['flux ', cfg.config, ', wavelength ',num2str(WLs(i)),'nm'],'fontsize',14);
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
%% plot example trajectories

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
for i = 2
    photon_id = unique(traj2{i}.id);
    hold on
    for b = 1:length(bounds)
        plot([0 size(cfg.vol,1)],[bounds(b) bounds(b)],'--');
    end
    title([' photon trajectories, no skull \newline wavelength = ', num2str(WLs(i)),', #DetPhoton = ',num2str(length(photon_id)) ])
    for k = 200:210
        ind_traj = find(traj2{i}.id == photon_id(k));
        % plot3(traj{i}.pos(ind_traj,1), traj{i}.pos(ind_traj,2), traj{i}.pos(ind_traj,3)), hold on
        % axis([0 300 0 300 0 300])
        plot(traj2{i}.pos(ind_traj,1), traj2{i}.pos(ind_traj,3))
    %     hold on
        scatter(traj2{i}.pos(ind_traj(1),1), traj2{i}.pos(ind_traj(1),3),50)
        axis([1 300 0 60])
        xticklabels({'5' '6' '7' '8' '9' '10'})
        yticks([0 20 40 60]), yticklabels({'0' '1' '2' '3'})
        xlabel('lateral distance, mm')
        ylabel('depth, mm')
        pause
    end
end


