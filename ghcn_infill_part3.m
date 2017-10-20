% THis is a script to produce and global SST analysis from
%  HadSST 2 data using RegEM
%   History : 
%	J.E.G . GaTech, Dec 4th 2007.
%	 nohup /usr/local/bin/matlab -nojvm -nodisplay < HadSST2_infill.m > HadSST2_infill.log &
%%%%%%%% 

load('ghcn_db.mat');

% lil' check on NINO3

%NINO3_x=find(lon >=210 & lon <=270);
%NINO3_y=find(lat>=-5 & lat <=5);
%NINO3=sq(nmean(nmean(HadSST2_raw(:,NINO3_y,NINO3_x),2),3));  
%plot(datenum(t),NINO3)
%datetick('x','yyyy');  big chunk missing in 1860's 

%% reshape matrix for RegEM purposes
load('ghcn_db.mat');
ng = length(ghcn);
ntime = 3744; %3732 corresponds to the maximum number of months
G = nan(ntime, ng); %preallocate w/ nan
nyears = floor(ntime/12);
for i = 1:ng
    years = ghcn(i).yr; ryears = years - 1700;
    nyears = length(years);
    for j = 1:nyears
        ind1 = min(ghcn(i).yr(j))-1700;
        G(1+(ind1-1)*12:ind1*12,i) = ghcn(i).anom(1+(j-1)*12:j*12);
    end
end

%  find valid points 
navail = sum(~isnan(G)); 
q10 = quantile(navail,0.2); %%q10 = 0 w/ .1
station = find(navail(4001:6000) >= ntime/10); no = numel(station);
% 1808 out of 2592 pts
X=G(:,station);


%% apply RegEM

OPTIONS.regress = 'ttls';
OPTIONS.regpar = 55; % this is dataset-dependent. Kaplan uses 80, so not an absurd first guess
OPTIONS.stagtol = 5e-2;
OPTIONS.maxit = 50;
%OPTIONS.inflation = 1;
OPTIONS.disp = 1;

%%           RUN REGEM !
disp('Preliminary step : TTLS')
[X0, M0, C0, Xerr0] = regem(X, OPTIONS);
save 'ghcn_part3_ttls_guess' X0 M0 C0 Xerr0 station

%load('RegEM_ttls_guess.mat');

% Feed that into an iridge regression 
OPTIONS2.Xmis0=X0;
OPTIONS2.C0=C0;
OPTIONS2.regress = 'iridge';
OPTIONS2.stagtol = 5e-3;
OPTIONS2.maxit = 200;
%OPTIONS2.inflation = 1; %  
OPTIONS2.disp = 1;

%   Call Multiple Ridge Regression 
disp('Infilling stage : ridge regression')
[Xhf, Mhf, Chf, Xerr] = regem(X,OPTIONS2);

%
save 'ghcn_part3_iridge.mat' Xhf Mhf Chf Xerr hadcrut4 station







