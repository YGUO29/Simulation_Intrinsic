%% set tissue parameters
load('D:\=code=\Simulation_Intrinsic\spectra.mat')
lambda = spectra(:,1); e_HbO = spectra(:,2); e_HbR = spectra(:,3); 
c_Hb   = 129.645; % hemoglobin (g/L) used in Hilman 2016
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