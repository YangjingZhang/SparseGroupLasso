%% processing the climate data NCEP/NCAR Reanalysis.
%% June 13, 2018 Zhang Yangjing.

%% The dataset contains monthly means of climate data measurements 
%% spread across the globe from 1948/1/1 to 2018/05/31 (845 months).
%% There are 7 predictive variables stored in 7 files:
%% air temperature:       air.mon.mean.nc
%% precipitable water:    pr_wtr.eatm.mon.mean.nc
%% relative humidity:     rhum.mon.mean.nc
%% pressure:              pres.mon.mean.nc
%% sea level pressure:    slp.mon.mean.nc
%% horizontal wind speed: uwnd.mon.mean.nc
%% vertical wind speed:   vwnd.mon.mean.nc

%% ncdisp(' ') can see what is in each file
%% ncread(' ') can load what is in each file

function Climatedata_process

%% load data, each is a 3-dimension matrix: 144*73*845
%% first dimension: longtitude \in R^144 = 2.5*[0,1:143] = [0, 2.5, 5, ..., 357.5]
%% second dimension: latitude \in R^73 = [90:-2.5:-90] = [90, 87.5, ..., 0, -2.5, ..., -90]
%% third dimension: months

dataair    = ncread('air.mon.mean.nc','air');
dataprw    = ncread('pr_wtr.eatm.mon.mean.nc','pr_wtr');
datahum    = ncread('rhum.mon.mean.nc','rhum');
datapres   = ncread('pres.mon.mean.nc','pres');
dataslp    = ncread('slp.mon.mean.nc','slp');
datauwind  = ncread('uwnd.mon.mean.nc','uwnd');
datavwind  = ncread('vwnd.mon.mean.nc','vwnd');
[lon,lat,mon] = size(dataair);

%% reshape to 2-dim matrix
%% row refers to places, column refers to months

dataair2   = reshape(dataair,(144*73),845);
dataprw2   = reshape(dataprw,(144*73),845);
datahum2   = reshape(datahum,(144*73),845);
datapres2  = reshape(datapres,(144*73),845);
dataslp2   = reshape(dataslp,(144*73),845);
datauwind2 = reshape(datauwind,(144*73),845);
datavwind2 = reshape(datavwind,(144*73),845);

%% stored in a 3-dim matrix, each page denotes one predictive variables
M1 = zeros(lon*lat,mon,7);
M1(:,:,1) = dataair2;
M1(:,:,2) = dataprw2;
M1(:,:,3) = datahum2;
M1(:,:,4) = datapres2;
M1(:,:,5) = dataslp2;
M1(:,:,6) = datauwind2;
M1(:,:,7) = datavwind2;

%%
A = zeros(mon,lon*lat*7);

for i = 1:mon
    tmp = M1(:,i,:);
    tmp = squeeze(tmp);
    tmp = tmp';
    A(i,:) = tmp(:);
end

%% Dakar's latidude: 14.6937N; longitude: -17.44406W (342.5594E)
%% we can find first dimension  = 138, second dimension = 31.
%% To analyze the air temperature around Dakar, we take the following b
b = squeeze(dataair(138,31,:));

save('Climatedata.mat','A','b');





