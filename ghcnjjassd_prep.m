%GHCN prep
%this script prepares treering and station data for screening
% -match treerings to closest station
% -find minimum/maximum years
%(this step takes awhile, which is why i decided to separate it
%from the rest of the screening process)
addpath(genpath('/home/geovault-02/avaccaro/Treerings/'))
nst=1;

%load GHCN data 
load('ghcn_db.mat');


%load treering db (trdb_crn)
load('tr_db_crn.mat');

%set vars
ns = length(ghcn); nt = length(trdb_crn);
r = 6370; %radius of earth in km
%nst = 4; %this is the number of stations that will be used to compute regional mean (1-4)

%set check counters = 0
nls = 0; %loc check
nts1 = 0; %max year check
nts2 = 0; %min year check
nmd = 0; %minimum data check
ns = 0; %final check



%find closest station and create tree-ring/station pair
c = pi/180;
for i = 1:nt %for each site
    for j = 1:ns %for each station
        
        lat1 = c*trdb_crn(i).lat;
        lon1 = c*trdb_crn(i).lon;
        lat2 = c*ghcn(j).lat;
        lon2 = c*ghcn(j).lon;
        
        if trdb_crn(i).locscreen == 1
            d3(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2,r);
        else
            d3(i,j)=NaN;
        end
    end
end


for i = 1:nt
    
    %store MXD meta data
    ghcnjjassd(i).id = trdb_crn(i).id;
    ghcnjjassd(i).lat = trdb_crn(i).lat;
    ghcnjjassd(i).lon = trdb_crn(i).lon;
    ghcnjjassd(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcnjjassd(i).locscreen == 1
        
        %compute regional means

        [smin idmin] = sort(d3(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of lowest 4
        
        %get stations' data and years
        st1 =  ghcn(minIds(1)).da;
        st1yrs = round(ghcn(minIds(1)).years);
        st2 = ghcn(minIds(2)).da;
        st2yrs = round(ghcn(minIds(2)).years);
        st3 = ghcn(minIds(3)).da;
        st3yrs = round(ghcn(minIds(3)).years);
        st4 = ghcn(minIds(4)).da;
        st4yrs = round(ghcn(minIds(4)).years);
        
        %find mins/maxs
        st1min = min(st1yrs); st2min = min(st2yrs);
        st3min = min(st3yrs); st4min = min(st4yrs);
        st1max = max(st1yrs); st2max = max(st2yrs);
        st3max = max(st3yrs); st4max = max(st4yrs);
        stmins = [st1min, st2min, st3min, st4min];
        stmaxs = [st1max, st2max, st3max, st4max];
        
        ghcnjjassd(i).stmin = max(stmins(1:nst));
        ghcnjjassd(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max (overlap)
        ghcnjjassd(i).styears = ghcnjjassd(i).stmin:ghcnjjassd(i).stmax;
        ny = length(ghcnjjassd(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcnjjassd(i).stmean = nan(1,ny);
        
        
        %get indices
        st1ind1 = ismember(ghcnjjassd(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcnjjassd(i).styears);
        st2ind1 = ismember(ghcnjjassd(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcnjjassd(i).styears);
        st3ind1 = ismember(ghcnjjassd(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcnjjassd(i).styears);
        st4ind1 = ismember(ghcnjjassd(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcnjjassd(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcnjjassd(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        %get matching MXD data
        ghcnjjassd(i).tryears = trdb_crn(i).yr;
        ghcnjjassd(i).trtemp = trdb_crn(i).x;
        ghcnjjassd(i).trmin = min(ghcnjjassd(i).tryears);
        ghcnjjassd(i).trmax = max(ghcnjjassd(i).tryears);
        
        %set overlap year range from min/max values
        ghcnjjassd(i).tmin = max(ghcnjjassd(i).trmin, ghcnjjassd(i).stmin);
        ghcnjjassd(i).tmax = min(ghcnjjassd(i).trmax, ghcnjjassd(i).stmax);
        ghcnjjassd(i).years = ghcnjjassd(i).tmin:ghcnjjassd(i).tmax;
  
        
    else
        ghcnjjassd(i).trmin = nan; ghcnjjassd(i).trmax = nan;
        ghcnjjassd(i).tryears = trdb_crn(i).yr;
        ghcnjjassd(i).stmin = nan; ghcnjjassd(i).stmax = nan;
        ghcnjjassd(i).stmean = nan;
        ghcnjjassd(i).years = nan(1,100);
        ghcnjjassd(i).locscreen = 0;
        ghcnjjassd(i).tmin = nan; ghcnjjassd(i).tmax = nan;
        ghcnjjassd(i).years = nan;
    end
    
    
end





%preallocate for decadal smoothing
for i = 1:nt
    if ghcnjjassd(i).locscreen == 1
        ghcnjjassd(i).stsmooth = nan(length(ghcnjjassd(i).styears),1);
        ghcnjjassd(i).trsmooth = nan(length(ghcnjjassd(i).tryears),1);
    end
end

for i = 1:nt
    if ghcnjjassd(i).locscreen == 1
        nntr = sum(~isnan(ghcnjjassd(i).trtemp));
        nnst = sum(~isnan(ghcnjjassd(i).stmean));
        if (nntr > 20 & nnst > 20)
            ghcnjjassd(i).mindatsm = 1;
            
            %decadally smoothed
            indmst = isfinite(ghcnjjassd(i).stmean);
            indmtr = isfinite(ghcnjjassd(i).trtemp);
            ghcnjjassd(i).stsmooth(indmst) = hepta_smooth(ghcnjjassd(i).stmean(indmst), 1/10);
            ghcnjjassd(i).trsmooth(indmtr) = hepta_smooth(ghcnjjassd(i).trtemp(indmtr), 1/10);
            
        else
            ghcnjjassd(i).mindatsm = 0;
        end
    else
        ghcnjjassd(i).mindatsm = 0;
    end
end


%time check using tmin/max values

for i = 1:nt
    if ghcnjjassd(i).locscreen == 0
        nls = nls + 1;
        ghcnjjassd(i).timescreen = 0;
    else
        ghcnjjassd(i).overlap = ghcnjjassd(i).tmax - ghcnjjassd(i).tmin;
        if ghcnjjassd(i).tmax < 1960
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcnjjassd(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcnjjassd(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcnjjassd(i).timescreen == 2
        
        %split years
        ghcnjjassd(i).yearsp1 = [ghcnjjassd(i).tmin:1960];
        ghcnjjassd(i).yearsp2 = [1961:ghcnjjassd(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcnjjassd(i).yearsp1);
        np2 = length(ghcnjjassd(i).yearsp2);
        

        ghcnjjassd(i).stsm1 = nan(np1,1);
        ghcnjjassd(i).stsm2 = nan(np2,1);
        ghcnjjassd(i).trsm1 = nan(np1,1);
        ghcnjjassd(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcnjjassd(i).yearsp1, ghcnjjassd(i).tryears);
        tr2ind = ismember(ghcnjjassd(i).yearsp2, ghcnjjassd(i).tryears);
        st1ind = ismember(ghcnjjassd(i).yearsp1, ghcnjjassd(i).styears);
        st2ind = ismember(ghcnjjassd(i).yearsp2, ghcnjjassd(i).styears);
        
        tr1ind2 = ismember(ghcnjjassd(i).tryears, ghcnjjassd(i).yearsp1);
        tr2ind2 = ismember(ghcnjjassd(i).tryears, ghcnjjassd(i).yearsp2);
        st1ind2 = ismember(ghcnjjassd(i).styears, ghcnjjassd(i).yearsp1);
        st2ind2 = ismember(ghcnjjassd(i).styears, ghcnjjassd(i).yearsp2);
        
        %get matching data
        ghcnjjassd(i).trsm1(tr1ind) = ghcnjjassd(i).trsmooth(tr1ind2);
        ghcnjjassd(i).trsm2(tr2ind) = ghcnjjassd(i).trsmooth(tr2ind2);
        ghcnjjassd(i).stsm1(st1ind) = ghcnjjassd(i).stsmooth(st1ind2);
        ghcnjjassd(i).stsm2(st2ind) = ghcnjjassd(i).stsmooth(st2ind2);
        
        ghcnjjassd(i).trsm = vertcat(ghcnjjassd(i).trsm1, ghcnjjassd(i).trsm2);
        ghcnjjassd(i).stsm = vertcat(ghcnjjassd(i).stsm1, ghcnjjassd(i).stsm2);
    
    else
        ghcnjjassd(i).mindatscreen = 0;
    end
    
end





save('ghcnjjassd.mat', 'ghcnjjassd')
clear all

