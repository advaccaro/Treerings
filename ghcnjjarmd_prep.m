%GHCN prep
%this script prepares treering and station data for screening
% -match treerings to closest station
% -find minimum/maximum years
%(this step takes awhile, which is why i decided to separate it
%from the rest of the screening process)
addpath(genpath('/home/geovault-02/avaccaro/Treerings/'))
nst=4;

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
    ghcnjjarmd(i).id = trdb_crn(i).id;
    ghcnjjarmd(i).lat = trdb_crn(i).lat;
    ghcnjjarmd(i).lon = trdb_crn(i).lon;
    ghcnjjarmd(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcnjjarmd(i).locscreen == 1
        
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
        
        ghcnjjarmd(i).stmin = max(stmins(1:nst));
        ghcnjjarmd(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max (overlap)
        ghcnjjarmd(i).styears = ghcnjjarmd(i).stmin:ghcnjjarmd(i).stmax;
        ny = length(ghcnjjarmd(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcnjjarmd(i).stmean = nan(1,ny);
        
        
        %get indices
        st1ind1 = ismember(ghcnjjarmd(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcnjjarmd(i).styears);
        st2ind1 = ismember(ghcnjjarmd(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcnjjarmd(i).styears);
        st3ind1 = ismember(ghcnjjarmd(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcnjjarmd(i).styears);
        st4ind1 = ismember(ghcnjjarmd(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcnjjarmd(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcnjjarmd(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        %get matching MXD data
        ghcnjjarmd(i).tryears = trdb_crn(i).yr;
        ghcnjjarmd(i).trtemp = trdb_crn(i).x;
        ghcnjjarmd(i).trmin = min(ghcnjjarmd(i).tryears);
        ghcnjjarmd(i).trmax = max(ghcnjjarmd(i).tryears);
        
        %set overlap year range from min/max values
        ghcnjjarmd(i).tmin = max(ghcnjjarmd(i).trmin, ghcnjjarmd(i).stmin);
        ghcnjjarmd(i).tmax = min(ghcnjjarmd(i).trmax, ghcnjjarmd(i).stmax);
        ghcnjjarmd(i).years = ghcnjjarmd(i).tmin:ghcnjjarmd(i).tmax;
  
        
    else
        ghcnjjarmd(i).trmin = nan; ghcnjjarmd(i).trmax = nan;
        ghcnjjarmd(i).tryears = trdb_crn(i).yr;
        ghcnjjarmd(i).stmin = nan; ghcnjjarmd(i).stmax = nan;
        ghcnjjarmd(i).stmean = nan;
        ghcnjjarmd(i).years = nan(1,100);
        ghcnjjarmd(i).locscreen = 0;
        ghcnjjarmd(i).tmin = nan; ghcnjjarmd(i).tmax = nan;
        ghcnjjarmd(i).years = nan;
    end
    
    
end





%preallocate for decadal smoothing
for i = 1:nt
    if ghcnjjarmd(i).locscreen == 1
        ghcnjjarmd(i).stsmooth = nan(length(ghcnjjarmd(i).styears),1);
        ghcnjjarmd(i).trsmooth = nan(length(ghcnjjarmd(i).tryears),1);
    end
end

for i = 1:nt
    if ghcnjjarmd(i).locscreen == 1
        nntr = sum(~isnan(ghcnjjarmd(i).trtemp));
        nnst = sum(~isnan(ghcnjjarmd(i).stmean));
        if (nntr > 20 & nnst > 20)
            ghcnjjarmd(i).mindatsm = 1;
            
            %decadally smoothed
            indmst = isfinite(ghcnjjarmd(i).stmean);
            indmtr = isfinite(ghcnjjarmd(i).trtemp);
            ghcnjjarmd(i).stsmooth(indmst) = hepta_smooth(ghcnjjarmd(i).stmean(indmst), 1/10);
            ghcnjjarmd(i).trsmooth(indmtr) = hepta_smooth(ghcnjjarmd(i).trtemp(indmtr), 1/10);
            
        else
            ghcnjjarmd(i).mindatsm = 0;
        end
    else
        ghcnjjarmd(i).mindatsm = 0;
    end
end


%time check using tmin/max values

for i = 1:nt
    if ghcnjjarmd(i).locscreen == 0
        nls = nls + 1;
        ghcnjjarmd(i).timescreen = 0;
    else
        ghcnjjarmd(i).overlap = ghcnjjarmd(i).tmax - ghcnjjarmd(i).tmin;
        if ghcnjjarmd(i).tmax < 1960
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcnjjarmd(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcnjjarmd(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcnjjarmd(i).timescreen == 2
        
        %split years
        ghcnjjarmd(i).yearsp1 = [ghcnjjarmd(i).tmin:1960];
        ghcnjjarmd(i).yearsp2 = [1961:ghcnjjarmd(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcnjjarmd(i).yearsp1);
        np2 = length(ghcnjjarmd(i).yearsp2);
        

        ghcnjjarmd(i).stsm1 = nan(np1,1);
        ghcnjjarmd(i).stsm2 = nan(np2,1);
        ghcnjjarmd(i).trsm1 = nan(np1,1);
        ghcnjjarmd(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcnjjarmd(i).yearsp1, ghcnjjarmd(i).tryears);
        tr2ind = ismember(ghcnjjarmd(i).yearsp2, ghcnjjarmd(i).tryears);
        st1ind = ismember(ghcnjjarmd(i).yearsp1, ghcnjjarmd(i).styears);
        st2ind = ismember(ghcnjjarmd(i).yearsp2, ghcnjjarmd(i).styears);
        
        tr1ind2 = ismember(ghcnjjarmd(i).tryears, ghcnjjarmd(i).yearsp1);
        tr2ind2 = ismember(ghcnjjarmd(i).tryears, ghcnjjarmd(i).yearsp2);
        st1ind2 = ismember(ghcnjjarmd(i).styears, ghcnjjarmd(i).yearsp1);
        st2ind2 = ismember(ghcnjjarmd(i).styears, ghcnjjarmd(i).yearsp2);
        
        %get matching data
        ghcnjjarmd(i).trsm1(tr1ind) = ghcnjjarmd(i).trsmooth(tr1ind2);
        ghcnjjarmd(i).trsm2(tr2ind) = ghcnjjarmd(i).trsmooth(tr2ind2);
        ghcnjjarmd(i).stsm1(st1ind) = ghcnjjarmd(i).stsmooth(st1ind2);
        ghcnjjarmd(i).stsm2(st2ind) = ghcnjjarmd(i).stsmooth(st2ind2);
        
        ghcnjjarmd(i).trsm = vertcat(ghcnjjarmd(i).trsm1, ghcnjjarmd(i).trsm2);
        ghcnjjarmd(i).stsm = vertcat(ghcnjjarmd(i).stsm1, ghcnjjarmd(i).stsm2);
    
    else
        ghcnjjarmd(i).mindatscreen = 0;
    end
    
end





save('ghcnjjarmd.mat', 'ghcnjjarmd')
clear all

