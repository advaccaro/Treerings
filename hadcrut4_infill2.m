% THis is a script to produce and global SST analysis from
%  HadCRUT4 data using RegEM
%   History : 
%	J.E.G . GaTech, Dec 4th 2007.
%	 nohup /usr/local/bin/matlab -nojvm -nodisplay < HadSST2_infill.m > HadSST2_infill.log &
%%%%%%%% 

load('hadcrut4.mat'); 

% lil' check on NINO3

%NINO3_x=find(lon >=210 & lon <=270);
%NINO3_y=find(lat>=-5 & lat <=5);
%NINO3=sq(nmean(nmean(HadSST2_raw(:,NINO3_y,NINO3_x),2),3));  
%plot(datenum(t),NINO3)
%datetick('x','yyyy');  big chunk missing in 1860's 

% reshape matrix for RegEM purposes
[nlon, nlat, ntime] = size(H.d);
nloc = nlat * nlon; %number of stations
hadcrut4 = reshape(H.d, [nloc, ntime]);
hadcrut4 = hadcrut4'; % time x space
hadcrut4 = double(hadcrut4);
hadcrut4(hadcrut4<-1000) = NaN;

%  find valid points 
navail = sum(~isnan(hadcrut4));
q10 = quantile(navail,0.1);
station = find(navail >= ntime/5); no = numel(station);
%ocean=find(~isnan(HadSST2_r(end,:,:))); % 1329 out of 2592 pts
X=hadcrut4(:,station);


% apply RegEM

OPTIONS.regress = 'ttls';
OPTIONS.regpar = 30; % this is dataset-dependent. Kaplan uses 80, so not an absurd first guess
OPTIONS.stagtol = 5e-2;
OPTIONS.maxit = 50;
OPTIONS.inflation = 1;
OPTIONS.disp = 1;

%%           RUN REGEM !
disp('Preliminary step : TTLS')
[X0, M0, C0, Xerr0] = regem(X, OPTIONS);
save 'hadcrut4_ttls_guess.mat' X0 M0 C0 Xerr0 station

%load('hadcrut4_ttls_guess.mat');

% Feed that into an iridge regression 
OPTIONS2.Xmis0=X0;
OPTIONS2.C0=C0;
OPTIONS2.regress = 'iridge';
OPTIONS2.stagtol = 5e-3;
OPTIONS2.maxit = 200;
OPTIONS2.inflation = 1; %  
OPTIONS2.disp = 1;

%   Call Multiple Ridge Regression 
disp('Infilling stage : ridge regression')
[Xhf, Mhf, Chf, Xerr] = regem(X,OPTIONS2);

%
save 'hadcrut4_iridge.mat' Xhf Mhf Chf Xerr hadcrut4 station







