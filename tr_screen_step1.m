%this script prepares treering and station data for screening
% -match treerings to closest station
% -find minimum/maximum years
%(this step takes awhile, which is why i decided to separate it
%from the rest of the screening process)
addpath(genpath('/Users/adam/Desktop/Treerings/'))


%load GHCN data 
load('ghcn_db.mat');


%load treering db (trdb_crn)
load('tr_db_crn.mat');

%set vars
ns = length(ghcn); nt = length(trdb_crn);
r = 6370;

nls = 0; %set locscreen counter = 0




%find closest station

for i = 1:nt
    if trdb_crn(i).locscreen == 1 %preliminary lat/lon check
        tcurve(i).locscreen = 1;
    else
        tcurve(i).locscreen = 0;
        nls = nls + 1;
    end
end

for i = 1:nt %for each site
    for j = 1:ns %for each station
        if tcurve(i).locscreen == 1
            %trdb_crn(i).lon(trdb_crn(i).lon<0) = trdb_crn(i).lon(trdb_crn(i).lon<0) + 360;
            %convert lat/lon to rads
            lat1 = (pi/180)*trdb_crn(i).lat;
            lon1 = (pi/180)*trdb_crn(i).lon;
            lat2 = (pi/180)*ghcn(j).lat;
            lon2 = (pi/180)*ghcn(j).lon;
          
            d(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2, r);
        else
            d(i,j)=NaN;
        end
    end
end

for i = 1:nt
    if tcurve(i).locscreen == 1
        
            [~,ind(i)]=min(d(i,:));%index of closest station
            
            %compile matching station and tree ring data

            tcurve(i).styears = round(ghcn(ind(i)).years);
            tcurve(i).tryears = trdb_crn(i).yr;
            
            tcurve(i).sttemp = ghcn(ind(i)).da;%data from closest station
            tcurve(i).trtemp = trdb_crn(i).x;
            
          
            
            
            
        else
            tcurve(i).tryears = trdb_crn(i).yr;
            na = length(tcurve(i).tryears);
            tcurve(i).tryears = nan(1,na);
            tcurve(i).trtemp = trdb_crn(i).x;
            nb = length(tcurve(i).trtemp);
            tcurve(i).sttemp = nan(1, nb);
    end
    
    
end
    

save('tcurve.mat','tcurve');
    
  
    
    