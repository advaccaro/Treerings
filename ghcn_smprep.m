%ghcn_smprep
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
nloc = length(ghcn); nt = length(trdb_crn);
r = 6370; %radius of earth in km
% nst = 4; %this is the number of stations that will be used to compute regional mean (1-4)

%set check counters = 0
nls = 0; %loc check
nts1 = 0; %max year check
nts2 = 0; %min year check
nmd = 0; %minimum data check
ns = 0; %final check

c = pi/180;

%find closest station and create tree-ring/station pair
for i = 1:nt %for each site
    for j = 1:nloc %for each hadcrut loc
        
        lat1 = c*trdb_crn(i).lat; %tree lat
        lon1 = c*trdb_crn(i).lon; %tree lon
        lat2 = c*ghcn(j).lat; %st lat
        lon2 = c*ghcn(j).lon; %st lon
        
        if trdb_crn(i).locscreen == 1
            d(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2, r);
        else
            d(i,j)=NaN;
        end
    end
end


for i = 1:nt;
    
    %store MXD meta data
    ghcnd(i).id = trdb_crn(i).id;
    ghcnd(i).lat = trdb_crn(i).lat;
    ghcnd(i).lon = trdb_crn(i).lon;
    ghcnd(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcnd(i).locscreen == 1
        
        %compute regional means

        [smin idmin] = sort(d(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of lowest 4
        
        %get stations' data and years
        st1 =  ghcn(minIds(1)).annual;
        st1yrs = round(ghcn(minIds(1)).years);
        st2 = ghcn(minIds(2)).annual;
        st2yrs = round(ghcn(minIds(2)).years);
        st3 = ghcn(minIds(3)).annual;
        st3yrs = round(ghcn(minIds(3)).years);
        st4 = ghcn(minIds(4)).annual;
        st4yrs = round(ghcn(minIds(4)).years);
        
        %find mins/maxs
        st1min = min(st1yrs(isfinite(st1))); st2min = min(st2yrs(isfinite(st2)));
        st3min = min(st3yrs(isfinite(st3))); st4min = min(st4yrs(isfinite(st4)));
        st1max = max(st1yrs(isfinite(st1))); st2max = max(st2yrs(isfinite(st2)));
        st3max = max(st3yrs(isfinite(st3))); st4max = max(st4yrs(isfinite(st4)));
        stmins = [st1min, st2min, st3min, st4min];
        stmaxs = [st1max, st2max, st3max, st4max];
        
        ghcnd(i).stmin = max(stmins(1:nst));
        ghcnd(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max (overlap)
        ghcnd(i).styears = ghcnd(i).stmin:ghcnd(i).stmax;
        ny = length(ghcnd(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcnd(i).stmean = nan(1,ny);
        
        
        %get indices
        st1ind1 = ismember(ghcnd(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcnd(i).styears);
        st2ind1 = ismember(ghcnd(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcnd(i).styears);
        st3ind1 = ismember(ghcnd(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcnd(i).styears);
        st4ind1 = ismember(ghcnd(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcnd(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcnd(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        %get matching MXD data
        ghcnd(i).tryears = trdb_crn(i).yr;
        ghcnd(i).trtemp = trdb_crn(i).x;
        ghcnd(i).trmin = min(ghcnd(i).tryears);
        ghcnd(i).trmax = max(ghcnd(i).tryears);
        
        %set overlap year range from min/max values
        ghcnd(i).tmin = max(ghcnd(i).trmin, ghcnd(i).stmin);
        ghcnd(i).tmax = min(ghcnd(i).trmax, ghcnd(i).stmax);
        ghcnd(i).years = ghcnd(i).tmin:ghcnd(i).tmax;
  
        
    else
        ghcnd(i).trmin = nan; ghcnd(i).trmax = nan;
        ghcnd(i).tryears = trdb_crn(i).yr;
        ghcnd(i).stmin = nan; ghcnd(i).stmax = nan;
        ghcnd(i).stmean = nan;
        ghcnd(i).years = nan(1,100);
        ghcnd(i).locscreen = 0;
        ghcnd(i).tmin = nan; ghcnd(i).tmax = nan;
        ghcnd(i).years = nan;
    end
    
    
end



%preallocate for decadal smoothing
for i = 1:nt
    if ghcnd(i).locscreen == 1
        ghcnd(i).stsmooth = nan(length(ghcnd(i).styears),1);
        ghcnd(i).trsmooth = nan(length(ghcnd(i).tryears),1);
    end
end

for i = 1:nt
    if ghcnd(i).locscreen == 1
        nntr = sum(~isnan(ghcnd(i).trtemp));
        nnst = sum(~isnan(ghcnd(i).stmean));
        if (nntr > 20 & nnst > 20)
            ghcnd(i).mindatsm = 1;
            
            %decadally smoothed
            indmst = isfinite(ghcnd(i).stmean);
            indmtr = isfinite(ghcnd(i).trtemp);
            ghcnd(i).stsmooth(indmst) = hepta_smooth(ghcnd(i).stmean(indmst), 1/10);
            ghcnd(i).trsmooth(indmtr) = hepta_smooth(ghcnd(i).trtemp(indmtr), 1/10);
            
        else
            ghcnd(i).mindatsm = 0;
        end
    else
        ghcnd(i).mindatsm = 0;
    end
end


%time check using tmin/max values

for i = 1:nt
    if ghcnd(i).locscreen == 0
        nls = nls + 1;
        ghcnd(i).timescreen = 0;
    else
        ghcnd(i).overlap = ghcnd(i).tmax - ghcnd(i).tmin;
        if ghcnd(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcnd(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcnd(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcnd(i).timescreen == 2
        
        %split years
        ghcnd(i).yearsp1 = [ghcnd(i).tmin:1960];
        ghcnd(i).yearsp2 = [1961:ghcnd(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcnd(i).yearsp1);
        np2 = length(ghcnd(i).yearsp2);
        

        ghcnd(i).stsm1 = nan(np1,1);
        ghcnd(i).stsm2 = nan(np2,1);
        ghcnd(i).trsm1 = nan(np1,1);
        ghcnd(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcnd(i).yearsp1, ghcnd(i).tryears);
        tr2ind = ismember(ghcnd(i).yearsp2, ghcnd(i).tryears);
        st1ind = ismember(ghcnd(i).yearsp1, ghcnd(i).styears);
        st2ind = ismember(ghcnd(i).yearsp2, ghcnd(i).styears);
        
        tr1ind2 = ismember(ghcnd(i).tryears, ghcnd(i).yearsp1);
        tr2ind2 = ismember(ghcnd(i).tryears, ghcnd(i).yearsp2);
        st1ind2 = ismember(ghcnd(i).styears, ghcnd(i).yearsp1);
        st2ind2 = ismember(ghcnd(i).styears, ghcnd(i).yearsp2);
        
        %get matching data
        ghcnd(i).trsm1(tr1ind) = ghcnd(i).trsmooth(tr1ind2);
        ghcnd(i).trsm2(tr2ind) = ghcnd(i).trsmooth(tr2ind2);
        ghcnd(i).stsm1(st1ind) = ghcnd(i).stsmooth(st1ind2);
        ghcnd(i).stsm2(st2ind) = ghcnd(i).stsmooth(st2ind2);
        
        ghcnd(i).trsm = vertcat(ghcnd(i).trsm1, ghcnd(i).trsm2);
        ghcnd(i).stsm = vertcat(ghcnd(i).stsm1, ghcnd(i).stsm2);
    
    else
        ghcnd(i).mindatscreen = 0;
    end
    
end





save('ghcnd.mat', 'ghcnd')

