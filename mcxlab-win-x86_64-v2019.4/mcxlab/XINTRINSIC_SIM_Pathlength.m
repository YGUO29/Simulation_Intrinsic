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
for i = 1:length(optical_para_label)
    eval([optical_para_label{i}, '= optical_para(:,i);']);
end
%% preparing the input data
figure
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
    cfg.size = [300 300 300];
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/2, cfg.size(2)/2)]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
    cfg.maxjumpdebug = cfg.nphoton;
    cfg.maxexitangle = NA;
    
% ========== define medium optical properties ==========
    % http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/CollinsAtlasMesh
    % format: [mua(1/mm) mus(1/mm) g n]
    cfg.prop=[0 0 1 1            % medium 0: the environment
%        0.019 7.8   0.89 1.55     % medium 1: skin & skull (n changed from 1.37 to 1.55)
%        0.004 0.009 0.89 1.37     % medium 2: CSF
       mua_HbT(ind(1)) 21 0.82, 1.37   ];   % medium 3: gray matter

    
    % define a 3 layer structure
    bounds1 = [1 1];
    cfg.vol = ones(cfg.size); % 3x3x3 mm volume
%     cfg.vol(:,:,bounds1(1)+1:bounds1(2))=2; % 0mm: skin & skull, 0.1mm CSF
%     cfg.vol(:,:,bounds1(2)+1:end)=2; % 2.9mm gray matter
    cfg.vol=uint8(cfg.vol);
    cfg.isreflect=0; % disable reflection at exterior boundary
    % GPU thread configuration
    cfg.autopilot=1;
    cfg.gpuid=1;

% ========== output control ==========
    cfg.issaveexit = 1;
    
% ========== Run simulation ==========
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    [f1,det1,vol,seeds,traj]=mcxlab(cfg);
    det1.ppath = det1.ppath.*cfg.unitinmm; % pathlength in units of mm
    toc;
    
    ExitAngle   = sin(vecnorm(det1.v(:,1:2)')./vecnorm(det1.v')); % sin_theta (NA)
    DetInd      = find(ExitAngle < cfg.maxexitangle);
    DetNumber   = length(DetInd);
    DetExitPosition = det1.p(DetInd,1:2);
    DetPathlength = det1.ppath(DetInd,1);
%% ========== plot figure ==========
    figure,
    I1 = log10(squeeze(sum(f1.data(:,size(cfg.vol,1)./2+1,:,:),4))');
    I1(I1<0) = 0;
    cmax1 = max(I1(:)); cmin1 = min(I1(:));
    contourf(I1, cmin1:0.5:cmax1);
    hold on
    plot([0 size(cfg.vol,1)],[bounds1(1) bounds1(1)],'--',...
        [0 size(cfg.vol,1)],[bounds1(2) bounds1(2)],'--');
    title(['flux with no reflection, no skull, wavelength ',num2str(WLs(i)),'nm']);
    set(gca,'clim',[cmin1 cmax1]);
    colorbar
    drawnow
    
    figure,
    h1 = histogram2(DetExitPosition(:,1),DetExitPosition(:,2),100,...
        'XBinLimits',[cfg.detpos(1)-10, cfg.detpos(1)+10],'YBinLimits',[cfg.detpos(2)-10, cfg.detpos(2)+10],...
        'DisplayStyle','tile','ShowEmptyBins','on');
    DetExitProfile = sum(h1.Values,2);
        [mm, ii] = max(DetExitProfile); cutoff = mm./2;
        [~,xx1] = min(abs(DetExitProfile(1:ii)-cutoff));
        [~,xx2] = min(abs(DetExitProfile(ii:end)-cutoff));    
    DetExitSize = (xx2 + ii - xx1).*cfg.unitinmm;
    title(['Exit position, FWHM: ',num2str(DetExitSize,3),'mm'])
%% plot pathlength
clear DetMeanPathlength
det2 = det1{1};
det2.ppath = det1{1}.ppath(det1{1}.DetInd,:).*cfg.unitinmm;

for i = 1:length(lambda)
    cfg.prop=[0 0 1 1            % medium 0: the environment
    brain_mua(i) gray_matter_mus(i) gray_matter_g(i), 1.37   % medium 1: gray matter
    brain_mua(i) white_matter_mus(i) white_matter_g(i), 1.37];   % medium 2: white matter
   
%       DetWeight = exp( -mua_HbT(ind(i)) *DetPathlength);
%       DetMeanPathlength(i) = sum(DetPathlength.*DetWeight) / sum(DetWeight(:));

    avgpath = mcxmeanpath(det2,cfg.prop);
    DetMeanPathlength(i) = sum(avgpath);
%     DetWeight = mcxdetweight(det1,cfg.prop);
%     DetPathlength = det1.ppath.*repmat(DetWeight(:),1,size(det1.ppath,2));
%     DetMeanPathlength(i) = cfg.unitinmm*sum(DetPathlength(:,2)) / sum(DetWeight(:));
%     h2 = histogram(DetPathlength);
end
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w'),
plot(lambda,DetMeanPathlength)
hold on, scatter(lambda(ind),DetMeanPathlength(ind),'k*')
xlabel('wavelength (nm)')
ylabel('mean pathlength (mm)')


