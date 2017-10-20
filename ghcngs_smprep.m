%GHCN prep
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
    ghcngsd(i).id = trdb_crn(i).id;
    ghcngsd(i).lat = trdb_crn(i).lat;
    ghcngsd(i).lon = trdb_crn(i).lon;
    ghcngsd(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcngsd(i).locscreen == 1
        
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
        
        ghcngsd(i).stmin = max(stmins(1:nst));
        ghcngsd(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max (overlap)
        ghcngsd(i).styears = ghcngsd(i).stmin:ghcngsd(i).stmax;
        ny = length(ghcngsd(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcngsd(i).stmean = nan(1,ny);
        
        
        %get indices
        st1ind1 = ismember(ghcngsd(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcngsd(i).styears);
        st2ind1 = ismember(ghcngsd(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcngsd(i).styears);
        st3ind1 = ismember(ghcngsd(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcngsd(i).styears);
        st4ind1 = ismember(ghcngsd(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcngsd(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcngsd(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        %get matching MXD data
        ghcngsd(i).tryears = trdb_crn(i).yr;
        ghcngsd(i).trtemp = trdb_crn(i).x;
        ghcngsd(i).trmin = min(ghcngsd(i).tryears);
        ghcngsd(i).trmax = max(ghcngsd(i).tryears);
        
        %set overlap year range from min/max values
        ghcngsd(i).tmin = max(ghcngsd(i).trmin, ghcngsd(i).stmin);
        ghcngsd(i).tmax = min(ghcngsd(i).trmax, ghcngsd(i).stmax);
        ghcngsd(i).years = ghcngsd(i).tmin:ghcngsd(i).tmax;
  
        
    else
        ghcngsd(i).trmin = nan; ghcngsd(i).trmax = nan;
        ghcngsd(i).tryears = trdb_crn(i).yr;
        ghcngsd(i).stmin = nan; ghcngsd(i).stmax = nan;
        ghcngsd(i).stmean = nan;
        ghcngsd(i).years = nan(1,100);
        ghcngsd(i).locscreen = 0;
        ghcngsd(i).tmin = nan; ghcngsd(i).tmax = nan;
        ghcngsd(i).years = nan;
    end
    
    
end





%preallocate for decadal smoothing
for i = 1:nt
    if ghcngsd(i).locscreen == 1
        ghcngsd(i).stsmooth = nan(length(ghcngsd(i).styears),1);
        ghcngsd(i).trsmooth = nan(length(ghcngsd(i).tryears),1);
    end
end

for i = 1:nt
    if ghcngsd(i).locscreen == 1
        nntr = sum(~isnan(ghcngsd(i).trtemp));
        nnst = sum(~isnan(ghcngsd(i).stmean));
        if (nntr > 20 & nnst > 20)
            ghcngsd(i).mindatsm = 1;
            
            %decadally smoothed
            indmst = isfinite(ghcngsd(i).stmean);
            indmtr = isfinite(ghcngsd(i).trtemp);
            ghcngsd(i).stsmooth(indmst) = hepta_smooth(ghcngsd(i).stmean(indmst), 1/10);
            ghcngsd(i).trsmooth(indmtr) = hepta_smooth(ghcngsd(i).trtemp(indmtr), 1/10);
            
        else
            ghcngsd(i).mindatsm = 0;
        end
    else
        ghcngsd(i).mindatsm = 0;
    end
end


%time check using tmin/max values

for i = 1:nt
    if ghcngsd(i).locscreen == 0
        nls = nls + 1;
        ghcngsd(i).timescreen = 0;
    else
        ghcngsd(i).overlap = ghcngsd(i).tmax - ghcngsd(i).tmin;
        if ghcngsd(i).tmax < 1960
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcngsd(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcngsd(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcngsd(i).timescreen == 2
        
        %split years
        ghcngsd(i).yearsp1 = [ghcngsd(i).tmin:1960];
        ghcngsd(i).yearsp2 = [1961:ghcngsd(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcngsd(i).yearsp1);
        np2 = length(ghcngsd(i).yearsp2);
        

        ghcngsd(i).stsm1 = nan(np1,1);
        ghcngsd(i).stsm2 = nan(np2,1);
        ghcngsd(i).trsm1 = nan(np1,1);
        ghcngsd(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcngsd(i).yearsp1, ghcngsd(i).tryears);
        tr2ind = ismember(ghcngsd(i).yearsp2, ghcngsd(i).tryears);
        st1ind = ismember(ghcngsd(i).yearsp1, ghcngsd(i).styears);
        st2ind = ismember(ghcngsd(i).yearsp2, ghcngsd(i).styears);
        
        tr1ind2 = ismember(ghcngsd(i).tryears, ghcngsd(i).yearsp1);
        tr2ind2 = ismember(ghcngsd(i).tryears, ghcngsd(i).yearsp2);
        st1ind2 = ismember(ghcngsd(i).styears, ghcngsd(i).yearsp1);
        st2ind2 = ismember(ghcngsd(i).styears, ghcngsd(i).yearsp2);
        
        %get matching data
        ghcngsd(i).trsm1(tr1ind) = ghcngsd(i).trsmooth(tr1ind2);
        ghcngsd(i).trsm2(tr2ind) = ghcngsd(i).trsmooth(tr2ind2);
        ghcngsd(i).stsm1(st1ind) = ghcngsd(i).stsmooth(st1ind2);
        ghcngsd(i).stsm2(st2ind) = ghcngsd(i).stsmooth(st2ind2);
        
        ghcngsd(i).trsm = vertcat(ghcngsd(i).trsm1, ghcngsd(i).trsm2);
        ghcngsd(i).stsm = vertcat(ghcngsd(i).stsm1, ghcngsd(i).stsm2);
    
    else
        ghcngsd(i).mindatscreen = 0;
    end
    
end





save('ghcngsd.mat', 'ghcngsd')

