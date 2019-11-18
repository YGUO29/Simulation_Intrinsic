%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
%
% simulate a 4-layer brain model using MCXLAB.
% 
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cfg;
%% set tissue parameters
load('C:\Users\guoyu\Documents\MATLAB\Simulation_Intrinsic\spectra.mat')
lambda = spectra(:,1); e_HbO = spectra(:,2); e_HbR = spectra(:,3); 
c_Hb   = 129; % hemoglobin (g/L) used in Hilman 2016

perc_blood = 0.03; % blood volume / tissue volume
sat_O2 = 0.75;
c_HbO = c_Hb*sat_O2; c_HbR = c_Hb*(1-sat_O2);
mua_HbR = perc_blood.*(e_HbR.*c_HbR.*2.303./64500)./10; % unitsL 1/mm
mua_HbO = perc_blood.*(e_HbO.*c_HbO.*2.303./64500)./10;
mua_HbT = mua_HbO+mua_HbR;

WLs = [470 530 590 625 730 850];
ind = floor(interp1(lambda,1:length(lambda), WLs));

figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2),

subplot(121), hold on, set(gca, 'YScale', 'log')
semilogy(lambda, mua_HbR), semilogy(spectra(:,1),mua_HbO),
semilogy(lambda,mua_HbT,'k')
scatter(lambda(ind),mua_HbT(ind),32)
for i = 1:length(ind)
text(lambda(ind(i))+5,mua_HbT(ind(i)),num2str(mua_HbT(ind(i)), 2))
end
title('Absorption coefficient')
legend('mua-HbR', 'mua-HbO2', 'mua-HbT')
subplot(122), hold on, set(gca, 'YScale', 'log')
semilogy(spectra(:,1),e_HbR), hold on, semilogy(spectra(:,1),e_HbO)
for i = 1:length(ind); xline(lambda(ind(i)),'--'); end
title('Moler extinction coefficient')
legend('\epsilon-HbR', '\epsilon-HbO2')


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
    cfg.size = [300 300 300];
    cfg.srcpos=[cfg.size(1)/2, cfg.size(2)/2, 1];
    cfg.srcdir=[0 0 1];
    cfg.detpos=[cfg.size(1)/2, cfg.size(2)/2, 0, min(cfg.size(1)/2, cfg.size(2)/2)]; % [x y z radius]
    cfg.maxdetphoton = cfg.nphoton;
    cfg.maxjumpdebug = cfg.nphoton;
    cfg.maxexitangle = 0.05;
% ========== define a 3 layer structure ==========
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
    
fig1 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
fig2 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');
fig3 = figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w');

for i = 1:length(ind)
% ========== define medium optical properties ==========
    % http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/CollinsAtlasMesh
    % format: [mua(1/mm) mus(1/mm) g n]
    cfg.prop=[0 0 1 1            % medium 0: the environment
%        0.019 7.8   0.89 1.55     % medium 1: skin & skull (n changed from 1.37 to 1.55)
%        0.004 0.009 0.89 1.37     % medium 2: CSF
       mua_HbT(ind(i)) 21 0.82, 1.37   ];   % medium 3: gray matter

% ========== Run 1st simulation ==========
    fprintf('running simulation .... this takes about 35 seconds on a GTX 470\n');
    tic;
    [f1,det1,vol,seeds,traj]=mcxlab(cfg);
    det1.ppath = det1.ppath.*cfg.unitinmm; % pathlength in units of mm
    toc;
    
    mua_change = cfg.prop(2,1).*0.1; % change absorption coefficient by 10% [1/mm]
    cfg.prop(2,1) = cfg.prop(2,1) .*0.9;
% ========== Run 2nd simulation ==========
    fprintf('running simulation ... this takes about 35 seconds on a GTX 470\n');
    tic;
    [f2,det2,vol,seeds,traj] = mcxlab(cfg);
    det2.ppath = det2.ppath.*cfg.unitinmm; % pathlength in units of mm
    toc;
    
%%
x = -cfg.size(1)*cfg.unitinmm/2 : cfg.unitinmm : cfg.size(1).*cfg.unitinmm/2 - cfg.unitinmm;
z = 0 : cfg.unitinmm : cfg.size(3)*cfg.unitinmm - cfg.unitinmm;
[X, Z] = meshgrid(x,z); 

figure(fig1),
subplot(2,3,i)
I1 = log10(squeeze(sum(f1.data(:,size(cfg.vol,1)./2+1,:,:),4))');
I1(I1<0) = 0;
if i == 1
cmax1 = max(I1(:)); cmin1 = min(I1(:));
end
contourf(X, Z, I1, cmin1:0.5:cmax1);
hold on
plot([0 size(cfg.vol,1)],[bounds1(1) bounds1(1)],'--',...
    [0 size(cfg.vol,1)],[bounds1(2) bounds1(2)],'--');
title(['flux with no reflection, no skull, wavelength ',num2str(WLs(i)),'nm']);
set(gca,'clim',[cmin1 cmax1]);
colorbar
drawnow
    
figure(fig2)
subplot(2,3,i)
I1 = squeeze(sum(f1.data(:,size(cfg.vol,1)./2+1,:,:),4))';
I2 = squeeze(sum(f2.data(:,size(cfg.vol,1)./2+1,:,:),4))';
II = log10(abs(I1 - I2)./mua_change); 
II(II<0) = 0;

% if i == 1
    cmax1 = max(II(:)); cmin1 = min(II(:));
% end
contourf(X, Z, II, cmin1:1:cmax1, 'ShowText','On');

xlabel('x direction (mm)'), ylabel('y direction (mm)')
title(['Depth sensitivity, wavelength ',num2str(WLs(i)),'nm']);
set(gca,'clim',[cmin1 cmax1]); colorbar, drawnow

figure(fig3), hold on
plot(mean(II,2))
title(['Depth sensitivity, wavelength ',num2str(WLs(i)),'nm']);
drawnow

end
ExitAngle   = sin(vecnorm(det1.v(:,1:2)')./vecnorm(det1.v')); % sin_theta (NA)
DetInd      = find(ExitAngle < cfg.maxexitangle);
DetNumber   = length(DetInd);
DetExitPosition = det1.p(DetInd,1:2);
DetPathlength = det1.ppath(DetInd,1);


%%
%% plot pathlength
clear DetMeanPathlength
det2 = det1;
det2.ppath = det1.ppath(DetInd,:);
cfg.vol_original = cfg.vol;
cfg.prop=[0 0 1 1            % medium 0: the environment
%        0.019 7.8   0.89 1.55     % medium 1: skin & skull (n changed from 1.37 to 1.55)
%        0.004 0.009 0.89 1.37     % medium 2: CSF
mua_HbT(i) 21 0.82, 1.37   
mua_HbT(i)*0.8 21 0.82, 1.37];   % medium 3: gray matter
for m = 1:300
    for n = 1:300
        cfg.vol = cfg.vol_original;
        cfg.vol(m, :, n) = 3;
%       DetWeight = exp( -mua_HbT(ind(i)) *DetPathlength);
%       DetMeanPathlength(i) = sum(DetPathlength.*DetWeight) / sum(DetWeight(:));
        avgpath = mcxmeanpath(det2,cfg.prop);
        DetMeanPathlength(i) = avgpath;
%     DetWeight = mcxdetweight(det1,cfg.prop);
%     DetPathlength = det1.ppath.*repmat(DetWeight(:),1,size(det1.ppath,2));
%     DetMeanPathlength(i) = cfg.unitinmm*sum(DetPathlength(:,2)) / sum(DetWeight(:));
%     h2 = histogram(DetPathlength);
    end
end
figure('DefaultAxesFontSize',18, 'DefaultLineLineWidth', 2,'color','w'),
plot(lambda,DetMeanPathlength)
hold on, scatter(lambda(ind),DetMeanPathlength(ind),'k*')
xlabel('wavelength (nm)')
ylabel('mean pathlength (mm)')