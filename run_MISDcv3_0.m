clear, clc, close all
name_bas= 'Po';

%%
load input_Po.mat

%% load basin configuration data
load PoBasin_config.mat

%%
sez_outlet = 15;
bas_check = 89;
%% run MISDc model
tic, [NS_sez,KGE_sez,KGE_out_IRR,Qsim_out,WW,WW2,VolIRR,SWE,PERC,EVAP]=MISD_v3_0('input_Po',BAS_PAR,EBRR_BASPAR,load('X_opt.txt'),sez_outlet,bas_check,ID_bas_app,1);toc

%% figure
close all
figure
set(gcf,'position',[246,250,1158,220])
hold on
plot(D,Q(:,bas_check),'g-','linew',3)
plot(D,Qsim_out(:,sez_outlet),'k.-')
legend('Obs','Sim')
axis([D(1) D(end) 0 8000])
datetick('x',12,'keeplimits')
box on, grid on

[~,RMSE_irr,~,RRQ_irr]=perf(Qsim_out(:,sez_outlet),Q(:,bas_check));
title(['Po Pontelagoscuro: - KGE=',num2str(klinggupta(Qsim_out(:,sez_outlet),Q(:,bas_check)),'%3.2f'),...
    ' | R= ',num2str(RRQ_irr,'%3.2f'),' | rRMSE= ',num2str(RMSE_irr,'%3.2f')])

ylabel('Discharge (m^3/sec)')

