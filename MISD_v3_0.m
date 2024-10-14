% -------------------------------------------------------------------------------
%   MISD_07(name,[Cf(Nbas),lambda(Nbas),eta(Nbas)],C,D,delta_T,FIG)
%   Nbas : number of subcatchments
%   Nsez : number of sections
%
%   MISD: SEMIDISTRIBUTED EVENT-BASED RAINFALL-RUNOFF MODEL
%   Author: Luca Brocca
%   Research Institute for Geo-Hydrological Protection, National Research Council
% -------------------------------------------------------------------------------

function [NS_sez,KGE_sez,KGE_out,Qsim_out,WW,WW2,VolIRR,SWE_pack,PERC2,E]=MISD_v3_0(input,BAS_PAR,EBRR_BASPAR,PAR,sez_outlet,bas_check,ID_bas_app,FIG,W2ini,Qchan_ini)

% Loading data

load([input,'.mat'])
TEMPER=T;
PIO_=P;


% Loading basin parameters
Nbas      = BAS_PAR(1); % number of subcatchments
Nsez      = BAS_PAR(2); % number of sections
Ninf      = BAS_PAR(3); % number of upstream inflows


delta_T=round(nanmean(diff(D))*24*10000)/10000;
dt      = 0.2;          % computation time step in hour


% Morphological data
DIST=EBRR_BASPAR(:,2:Nsez+1);               % catchment distance to outlet sections
Ab=EBRR_BASPAR(:,2+Nsez);Ab(1:Ninf)=[];     % catchments area
A_DD=EBRR_BASPAR(:,3+Nsez);A_DD(1:Ninf)=[]; % catchments type (1: concentrater, 2: distributed)

% Initialization
M = size(D,1);
QQsim=zeros(M,Nbas);            % catchments runoff (mm/delta_T)
QQQsim=zeros(M,Nbas+Ninf,Nsez); % sections runoff (mm/delta_T)

W_max     = 200;     % FIXED WATER CAPACITY 1st LAYER

% Main ROUTINE: SWB
W_p    = PAR(1,:);  % initial conditions, fraction of W_max (0-1)
W_max2 = PAR(2,:);  % total water capacity of 2nd layer
m2     = PAR(3,:);  % exponent of drainage for 1st layer
Ks     = PAR(4,:);  % hydraulic conductivity for 1st layer
Kc     = PAR(6,:);  % parameter of potential evapotranspiration
alpha  = PAR(7,:);  % exponent runoff
Cm     = PAR(8,:);  % Snow module parameter degree-day
m22    = PAR(9,:);  % exponent of drainage for 2nd layer
Ks2    = PAR(10,:); % hydraulic conductivity for 2nd layer

Ks      = Ks.*delta_T;  % mm/h --> mm/delta_T
Ks2     = Ks2.*delta_T; % mm/h --> mm/delta_T


%     D=basin_data{i}{:,1}; PIO_= basin_data{i}{:,2}; TEMPER=basin_data{i}{:,3}; Q =  basin_data{i}{:,4};
MESE=month(D);


% Snow Module
[PIO,SWE,SWE_pack,snowfall]=snow_model(PIO_, TEMPER, -0.5, 0.5, Cm);

% Potential Evapotranspiration parameter
L=[0.2100;0.2200;0.2300;0.2800;0.3000;0.3100;
    0.3000;0.2900;0.2700;0.2500;0.2200;0.2000];
Ka=1.26;
EPOT=(TEMPER>0).*(Kc.*(Ka*L(MESE).*(0.46.*TEMPER+8)-2))./(24/delta_T);
% EPOT=Kc.*EPOT;

% Initialization
BF=zeros(size(PIO));
QS=NaN(size(PIO));
WW=zeros(size(PIO));
WW2=zeros(size(PIO));
IE=zeros(size(PIO));
PERC=zeros(size(PIO));
PERC2=zeros(size(PIO));
SE =zeros(size(PIO));
SE2 = zeros(size(PIO));
VolIRR = zeros(size(PIO));

W=zeros(size(PIO));
W2=zeros(size(PIO));
W(1,:)=W_p*W_max;
if nargin>8
    W2= W2ini;
else
    W2(1,:)=W_p.*W_max2;
end

for t=1:M
    
    if t>1, W(t,:)=W(t-1,:);  W2(t,:)=W2(t-1,:); end
    IE(t,:)=PIO(t,:).*((W(t,:)/W_max).^alpha);
    E(t,:)= EPOT(t,:).*W(t,:)/W_max;
    
    if  sum(W2(t,:)<W_max2)>=1
        PERC(t,:)=Ks.*(W(t,:)/W_max).^(m2);
    else
        PERC(t,:)=0;
    end
    
    if sum(PERC(t,:)>0.75*W_max)>=1
        PERC(t,find(PERC(t,:)>0.75*W_max))=0.75*W_max;
    end
    
    
    PERC2(t,:)=Ks2.*(W2(t,:)./W_max2).^(m22);
    WUse = load('prelievi_Po.txt');

    
    W(t,:)=max(0,W(t,:)+(PIO(t,:)-IE(t,:)-PERC(t,:)-E(t,:))+SWE(t,:));
    W2(t,:)=max(0,W2(t,:)+PERC(t,:)-PERC2(t,:)-WUse(t,7));
    
    if sum(W(t,:)>=W_max)>=1
        SE(t,find(W(t,:)>=W_max))=W(t,find(W(t,:)>=W_max))-W_max;
        W(t,find(W(t,:)>=W_max))=repmat(W_max,1,sum(W(t,:)>=W_max));
    else
        SE(t,:)=0;
    end
    
    if sum(W2(t,:)>=W_max2)>=1
        SE2(t,find(W2(t,:)>=W_max2))=W2(t,find(W2(t,:)>=W_max2))-W_max2(find(W2(t,:)>=W_max2));
         W2(t,find(W2(t,:)>=W_max2))=W_max2(find(W2(t,:)>=W_max2));
    else
        SE2(t,:)=0;
    end
    
    
    WW(t,:)= W(t,:)./W_max;
    WW2(t,:)=W2(t,:)./W_max2;
    
    
    % irrigation contribution 
%     DISTR_Irr= [4;5;8;13;14;15;17;18;19;20;22;23;24;29;30;32;33;34;37;38;39;40;42;43;46;49;50;52;54;55;57;58;59;60;61;62;63;65;66;68;69;70;71;72;75;77;78;79;80;81;82;83;85;87;88;89];
    DISTR_Irr = find(ID_irr==1);

    Thres_irr = 0.5;

    Thres_ini = 0.23;
    
for ll=1:length(DISTR_Irr)
    
    if WW(t,DISTR_Irr(ll))<Thres_ini && MESE(t)==5|MESE(t)==6|MESE(t)==7|MESE(t)==8
%         disp('IRRIGATION')
        if WW2(t,DISTR_Irr(ll))*W_max2(DISTR_Irr(ll))>(Thres_irr-WW(t,DISTR_Irr(ll)))*W_max
            W(t,DISTR_Irr(ll))= Thres_irr*W_max;
            W2(t,DISTR_Irr(ll))= WW2(t,DISTR_Irr(ll)).*W_max2(DISTR_Irr(ll))-((Thres_irr-WW(t,DISTR_Irr(ll))).*W_max);
            
            WW2(t,DISTR_Irr(ll))= W2(t,DISTR_Irr(ll))./W_max2(DISTR_Irr(ll));
            WW(t,DISTR_Irr(ll))= Thres_irr;
        else
            W(t,DISTR_Irr(ll))= 0.1*(WW2(t,DISTR_Irr(ll))*W_max2(DISTR_Irr(ll)));
            W2(t,DISTR_Irr(ll))= 0.9*(WW2(t,DISTR_Irr(ll))*W_max2(DISTR_Irr(ll)));
            
%             WW2(t,DISTR_Irr(ll))= W2(t,DISTR_Irr(ll))./W_max2;
%             WW(t,DISTR_Irr(ll))= W(t,DISTR_Irr(ll))./W_max;
            WW2(t,DISTR_Irr(ll))= W2(t,DISTR_Irr(ll))./W_max2(DISTR_Irr(ll));
            WW(t,DISTR_Irr(ll))= W(t,DISTR_Irr(ll))./W_max;

        end
        VolIRR(t,DISTR_Irr(ll))= W(t,DISTR_Irr(ll));
    end
end



    
    if W(t,:)<0,  W(t,:)=0; end
    if W2(t,:)<0, W2(t,:)=0; end

    
    % Runoff contribution
    BF(t,:)=Ks2.*((W(t,:)+W2(t,:))./(W_max+W_max2)).^(m22);
    QS(t,:)=IE(t,:)+SE(t,:)+SE2(t,:);       % surface flow
    
end

if size(QS,2)>Nbas
    for i=1:Nbas
        
        peffmean(:,i) = nanmean(QS(:,find(ID_PTbas(1,:)==i)),2);
        BF_mean(:,i) =  nanmean(BF(:,find(ID_PTbas(1,:)==i)),2);
    end
        
else
    peffmean = QS;
    BF_mean  = BF;
    
end



for i=1:Nbas
    gamma  = PAR(5,i);  % coefficient lag-time relationship
    C      = PAR(11,i); % Celerity
    Diff   = PAR(12,i); % Diffusivity
    
    % Convolution (GIUH and NASH)
    if Ab(i)>0.0
        if A_DD(i)==1
            IUH=GIUH(gamma,Ab(i),dt,delta_T)*dt;
        else
            IUH=NASH(2*gamma,Ab(i),dt,delta_T,1)*dt;
        end
        % interpolation for dt time step
        peffint=interp1(1:M,peffmean(:,i),1:dt:M)';
        BFint=interp1(1:M,BF_mean(:,i),1:dt:M)';
        % convolution
        temp1=conv(IUH,peffint);
        temp2=conv(NASH(2*gamma,Ab(i),dt,delta_T,1)*dt,BFint);
        % saving outputs
        QQBF(:,i)= temp2(1:round(1/dt):M*round(1/dt));
        QQsim(:,i)= temp1(1:round(1/dt):M*round(1/dt))...
            + temp2(1:round(1/dt):M*round(1/dt)); % in delta_T
    else % if very low area (< 0 km^2) no convolution
        QQBF(:,i) = BF_mean(:,i);
        QQsim(:,i)= peffmean(:,i)+BF_mean(:,i);
    end
    
    % Convolution (Hayami)
    for j=1:Nsez
        if DIST(i+Ninf,j)>0
            g=hayami(dt,DIST(i+Ninf,j),C,Diff,delta_T)*dt;
            % interpolation for dt time step
            QQsimint=interp1(1:M,QQsim(:,i),1:dt:M)';
            QQBFint=interp1(1:M,QQBF(:,i),1:dt:M)';
            % convolution
            temp1=(conv(g',QQsimint));
            temp2=(conv(g',QQBFint));
            % saving outputs
            QQQBF(:,i+Ninf,j)=temp2(1:round(1/dt):M*round(1/dt)); % in delta_T
            QQQsim(:,i+Ninf,j)=temp1(1:round(1/dt):M*round(1/dt));
        elseif DIST(i+Ninf,j)==0
            QQQBF(:,i+Ninf,j)=QQBF(:,i);
            QQQsim(:,i+Ninf,j)=QQsim(:,i);
        else % if negative distance no contribution
            QQQBF(:,i+Ninf,j)=zeros(M,1,1);
            QQQsim(:,i+Ninf,j)=zeros(M,1,1);
        end
    end
end

% Calculation of catchments outlet discharge in m^3/s
for i=1:Nbas
    QQBF(:,i)= QQBF(:,i).*(Ab(i)./delta_T./3.6);
    QQsim(:,i)=QQsim(:,i).*(Ab(i)./delta_T./3.6);
end

% Calculation of sections discharge in m^3/s
for i=1:Nbas
    for j=1:Nsez
        QQQBF(:,i+Ninf,j)=QQQBF(:,i+Ninf,j).*(Ab(i)./delta_T./3.6);
        QQQsim(:,i+Ninf,j)=QQQsim(:,i+Ninf,j).*(Ab(i)./delta_T./3.6);
        
    end
end


% Convolution (Hayami) of upstream discharge
for j=1:Nsez
    for i=1:Ninf
        if DIST(i,j)>0.0
            g=hayami(dt,DIST(i,j),C,Diff,delta_T)*dt;
            QQsimint=interp1(1:M,QM(:,i),1:dt:M)';
            temp1=conv(g',QQsimint);
            QQQsim(:,i,j)=temp1(1:round(1/dt):M*round(1/dt));
        end
    end
end

% Calculation of outlet sections discharge
for j=1:Nsez
    QBF1(:,j)=(nansum(QQQBF(:,:,j)'))';
    Qsim1(:,j)=(nansum(QQQsim(:,:,j)'))';
    KGE_sez(j,1) =klinggupta(Qsim1(:,j),Q(:,ID_bas_app(j)));
    NS_sez(j,1)=1-nansum((Qsim1(:,j)-Q(:,ID_bas_app(j))).^2)./nansum((Q(:,ID_bas_app(j))-nanmean(Q(:,ID_bas_app(j)))).^2);
    
end
Qsim=Qsim1(:,sez_outlet);
QBF_out=[QBF1,QQBF];
Qsim_out=[Qsim1,QQsim];
% Qsim_out(1,:)=Qsim_out(1,:);

% identification Qobs

Qobs =  Q(:,bas_check);

% Calculation of model performance
RMSE=nanmean((Qsim-Qobs).^2).^0.5;
NS=1-nansum((Qsim-Qobs).^2)./nansum((Qobs-nanmean(Qobs)).^2);
ANSE=1-nansum((Qobs+nanmean(Qobs)).*(Qsim-Qobs).^2)./...
    nansum((Qobs+nanmean(Qobs)).*(Qobs-nanmean(Qobs)).^2);
NS_radQ=1-nansum((sqrt(Qsim)-sqrt(Qobs)).^2)./nansum((sqrt(Qobs)-nanmean(sqrt(Qobs))).^2);
NS_lnQ=1-nansum((log(Qsim+0.00001)-log(Qobs+0.00001)).^2)...
    ./nansum((log(Qobs+0.00001)-nanmean(log(Qobs+0.00001))).^2);
X=[Qsim,Qobs]; X(any(isnan(X)'),:) = [];
RRQ=corrcoef(X).^2; RQ=RRQ(2);
KGE_out=klinggupta(Qsim,Qobs);

% -------------------------------------------------------------------------------
% Calculation of Geomorphological Instantaneous Unit Hydrograph
% -------------------------------------------------------------------------------
function IUH=GIUH(gamma,Ab,dt,deltaT)

Lag=(gamma*1.19*Ab^0.33)/deltaT;
hp=0.8/Lag;
data=load('GIUH');
t=data(:,1)*Lag;IUH_0=data(:,2)*hp;
ti=0:dt:max(t);
IUH=interp1(t,IUH_0,ti)';
IUH=IUH./(sum(IUH).*dt);

% -------------------------------------------------------------------------------
% Calculation of Nash Instantaneous Unit Hydrograph
% -------------------------------------------------------------------------------
function IUH=NASH(gamma,Ab,dt,deltaT,n)

K=(gamma*1.19*Ab^0.33)/deltaT;
Tmax=100; % (h)
time=0.00001:dt:Tmax;
IUH=((time/K).^(n-1).*exp(-time/K)/factorial(n-1)/K)';
IUH=IUH./(sum(IUH).*dt);

% -------------------------------------------------------------------------------
% Calculation of Hayami function (diffusive routing)
% -------------------------------------------------------------------------------
function g=hayami(dt,L,C,D,deltaT)

C=C*deltaT;D=D*deltaT;
Tmax=100; % (h)
tt=0.00001:dt:Tmax;
g=(L./sqrt(4*pi*D.*tt.^3)).*exp(-((L-C.*tt).^2)./(4*D.*tt));
g=g./(sum(g).*dt);

%--------------------------------------------------------------------------
% Snow accumulation-melting MODEL
%--------------------------------------------------------------------------

function [rainfall,SWE_melting,SWE_snowpack,snowfall]=snow_model(precipitation, temperature, temp_min, temp_max, Cm)

rainfall = zeros(size(precipitation));
snowfall = zeros(size(precipitation));
SWE_snowpack = zeros(size(precipitation));
SWE_melting = zeros(size(precipitation));

% The precipitation is divided into rainfall and snowfall
% REFERENCES:
% U.S. Army Corps of Engineers (1956)
% Boscarello, L., Ravazzani, G., Pellegrini, M., Dedieu, J. P., & Mancini, M. (2014). Calibration of hydrological model FEST from MODIS images in Alpine Catchments. Politecnico di Milano, Dipartimento di Ingegneria Idraulica, Ambientale, Infrastrutture viarie, Rilevamento.
% Degree Day Method (Mockus, 1964)

% INITIALIZATION
for jj=1:size(precipitation,2)
    if precipitation(1,jj) == NaN || temperature(1,jj) == NaN
        rainfall(1,jj) = NaN;
        snowfall(1,jj) = NaN;
    elseif temperature(1,jj) <= temp_min
        snowfall(1,jj) = precipitation(1,jj);
        rainfall(1,jj) = 0;
        SWE_snowpack(1,jj) = snowfall(1,jj); % [mm]
        SWE_melting(1,jj) = 0; % [mm]
    elseif temperature(1,jj) >= temp_max
        snowfall(1,jj) = 0;
        rainfall(1,jj) = precipitation(1,jj);
        SWE_snowpack(1,jj) = 0; % [mm]
        SWE_melting(1,jj) = 0; % [mm]
    else
        rainfall(1,jj) = precipitation(1,jj) * ((temperature(1,jj)-temp_min)/(temp_max-temp_min));
        snowfall(1,jj) = precipitation(1,jj) - rainfall(1,jj);
        SWE_snowpack(1,jj) = snowfall(1,jj);
        SWE_melting(1,jj) = 0;
    end
    
    
    for i=2:size(precipitation,1)
        if precipitation(i,jj) == NaN || temperature(i,jj) == NaN
            % if the data is missing, it is equal to NaN
            rainfall(i,jj) = NaN;
            snowfall(i,jj) = NaN;
        elseif temperature(i,jj) <= temp_min
            % if the temperature is less than the low threshold,
            % the precipitation is entirely snowfall
            rainfall(i,jj) = 0;
            snowfall(i,jj) = precipitation(i,jj);
            SWE_snowpack(i,jj) = SWE_snowpack(i-1,jj) + snowfall(i,jj);
            SWE_melting(i,jj) = 0;
        elseif temperature(i,jj) > temp_max
            % if the temperature is more than the high threshold,
            % the precipitation is entirely rainfall
            rainfall(i,jj) = precipitation(i,jj);
            snowfall(i,jj) = 0;
            SWE_melting(i,jj) = Cm(jj).* (temperature(i,jj) - temp_max);
            % h_melting(i,1) = rho_water * SWE_melting(i,1) / rho_snow;
            % Check the snowpack SWE
            if SWE_snowpack(i-1,jj) >= SWE_melting(i,jj)
                SWE_snowpack(i,jj) = SWE_snowpack(i-1,jj) - SWE_melting(i,jj);
            else
                SWE_melting(i,jj) = SWE_snowpack(i-1,jj);
                SWE_snowpack(i,jj) = 0;
            end
        else
            rainfall(i,jj) = precipitation(i,jj) * ((temperature(i,jj)-temp_min)/(temp_max-temp_min));
            snowfall(i,jj) = precipitation(i,jj) - rainfall(i,jj);
            SWE_snowpack(i,jj) = SWE_snowpack(i-1,jj) + snowfall(i,jj);
            SWE_melting(i,jj) = 0;
        end
    end
end

