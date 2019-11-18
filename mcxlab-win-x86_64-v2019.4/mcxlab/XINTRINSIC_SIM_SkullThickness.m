%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
% 
% Depth sensitivity based on 2011 Tian paper
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cfg;
%% set tissue parameters
load('C:\Users\guoyu\Documents\MATLAB\Simulation_Intrinsic\optical_para.mat')
for i = 1:length(optical_para_label)
    eval([optical_para_label{i}, '= optical_para(:,i);']);
end
%% preparing the input data

    % set seed to make the simulation repeatible
    cfg.seed=hex2dec('623F9A9E'); 

    cfg.nphoton=5e7;
    cfg.unitinmm = 0.05; % define units, 50 um
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
    
% ========== define properties ==========
    cfg.prop=[0 0 1 1            % medium 0: the environment
       0.019 mus_skull(2)   0.89 1.55     % medium 1: skin & skull (n changed from 1.37 to 1.55)
%        0.004 0.009 0.89 1.37     % medium 2: CSF
       mua_HbT(ind(2)) 21 0.82, 1.37    % medium 3: gray matter
       mua_HbT(ind(2)) 40.9 0.84 1.37];   % medium 4: white matter

%     bounds1 = [floor(0.7/cfg.unitinmm) floor(0.7/cfg.unitinmm)+floor(0.1/cfg.unitinmm)]
%     cfg.vol(:,:,bounds1(1)+1:bounds1(2))=2; % 0.7mm: skin & skull, 0.1mm CSF
%     cfg.vol(:,:,bounds1(2)+1:end)=3; % the rest: gray matter
%     cfg.vol=uint8(cfg.vol);
    cfg.isreflect=0; % disable reflection at exterior boundary
    % GPU thread configuration
    cfg.autopilot=1;
    cfg.gpuid=1;
% ========== output control ==========
    cfg.issaveexit = 1;
%%    
f1 = cell(1,length(ind));
det1 = cell(1,length(ind));
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')
% for i = 1
for i = 2:length(ind)
% ========== define a 3 layer structure ==========
    cfg.vol = ones(cfg.size); % 15 mm volume
    bounds = [i.*round(0.7/cfg.unitinmm), i.*round(0.7/cfg.unitinmm) + round(1.3/cfg.unitinmm)]; 
    cfg.vol(:,:,bounds(1)+1:bounds(2)) = 2;
    cfg.vol(:,:,bounds(2)+1:end) = 3;
   
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    % f1=mcxlab(cfg);
    [f1{i},det1{i},vol,seeds,traj]=mcxlab(cfg);
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
    plot([0 size(cfg.vol,1)],[bounds(1) bounds(1)],'--',...
        [0 size(cfg.vol,1)],[bounds(2) bounds(2)],'--');
    title(['flux with no reflection, no skull, wavelength ',num2str(WLs(i)),'nm']);
    set(gca,'clim',[cmin1 cmax1]);
    colorbar
    drawnow
% ========== plot depth sensitivity ==========
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
DetWeight_skull = zeros(length(det1),length(ind));
DetWeight_brain = DetWeight_skull;
for j = 1:length(ind)
    cfg.prop=[0 0 1 1            % medium 0: the environment
    0.019 mus_skull(j)   0.89 1.55     % medium 1: skin & skull (n changed from 1.37 to 1.55)
    %        0.004 0.009 0.89 1.37     % medium 2: CSF
    mua_HbT(ind(j)) 21 0.82, 1.37    % medium 3: gray matter
    mua_HbT(ind(j)) 40.9 0.84 1.37];   % medium 4: white matter
    for i = 1:length(det1)
        det_temp = det1{i}; det2 = det_temp;
        det_temp.ppath = det_temp.ppath(det_temp.DetInd,:);
        DetWeight = sum(mcxdetweight(det_temp,cfg.prop));
        pskull = find(det_temp.ppath(:,2) == 0); 
        det2.ppath = det_temp.ppath(pskull,:);
        DetWeight_skull(i,j) = sum(mcxdetweight(det2,cfg.prop))/DetWeight;
        
        pbrain = setdiff(1:size(det_temp.ppath,1),pskull); 
        det2.ppath = det_temp.ppath(pbrain,:);
        DetWeight_brain(i,j) = sum(mcxdetweight(det2,cfg.prop))/DetWeight;
        
%         DetPhoton_skull(i,j) = length(pskull)/size(det_temp.ppath,1);
%         DetPhoton_brain(i,j) = length(pbrain)/size(det_temp.ppath,1);

    end
end

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')
for k = 1:length(ind)
    subplot(2,3,k), hold on
    plot(DetWeight_skull(:,k)), plot(DetWeight_brain(:,k))
end

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w')
for k = 1:length(ind)
    subplot(2,3,k), hold on
    plot(DetPhoton_skull(:,k)), plot(DetPhoton_brain(:,k))
end

set(gca,'XTickLabels',{'0-400um', '400-800um', '800-1200um', '1200-1600um', '1600-15000um'})

