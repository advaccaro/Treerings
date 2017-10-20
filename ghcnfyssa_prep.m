%ghcn_prep
%GHCN prep
%this script prepares treering and station data for screening
% -match treerings to closest station
% -find minimum/maximum years
%(this step takes awhile, which is why i decided to separate it
%from the rest of the screening process)
addpath(genpath('/home/geovault-02/avaccaro/Treerings'))


%load GHCN data 
load('ghcn_db.mat');


%load treering db (trdb_crn)
load('tr_db_crn.mat');

%set vars
ng = length(ghcn); nt = length(trdb_crn);
r = 6370; %radius of earth in km
nst = 1; %this is the number of stations that will be used to compute regional mean (1-4)

%set check counters = 0
nls = 0; %loc check
nts1 = 0; %max year check
nts2 = 0; %min year check
nmd = 0; %minimum data check
ns = 0; %final check




%find closest station and create tree-ring/station pair
c = pi/180;
for i = 1:nt %for each site
    for j = 1:ng %for each GHCN station
        
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


for i = 1:nt;
    %store MXD meta data
    ghcnfyssa(i).id = trdb_crn(i).id;
    ghcnfyssa(i).lat = trdb_crn(i).lat;
    ghcnfyssa(i).lon = trdb_crn(i).lon;
    ghcnfyssa(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcnfyssa(i).locscreen == 1
        
        
        %compute regional mean
        
        [smin idmin] = sort(d3(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of closest 4 stations
        %[~,ind(i)]=min(d(i,:));%index of closest station
        
        
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
        st1min = min(st1yrs); st2min = min(st2yrs);
        st3min = min(st3yrs); st4min = min(st4yrs);
        st1max = max(st1yrs); st2max = max(st2yrs);
        st3max = max(st3yrs); st4max = max(st4yrs);
        stmins = [st1min, st2min, st3min, st4min];
        stmaxs = [st1max, st2max, st3max, st4max];
        
        ghcnfyssa(i).stmin = max(stmins(1:nst));
        ghcnfyssa(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max
        ghcnfyssa(i).styears = ghcnfyssa(i).stmin:ghcnfyssa(i).stmax;
        ny = length(ghcnfyssa(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcnfyssa(i).stmean = nan(1,ny);
         
        
        %get indices
        st1ind1 = ismember(ghcnfyssa(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcnfyssa(i).styears);
        st2ind1 = ismember(ghcnfyssa(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcnfyssa(i).styears);
        st3ind1 = ismember(ghcnfyssa(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcnfyssa(i).styears);
        st4ind1 = ismember(ghcnfyssa(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcnfyssa(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcnfyssa(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        
        %get matching MXD data
        ghcnfyssa(i).tryears = trdb_crn(i).yr;
        ghcnfyssa(i).trtemp = trdb_crn(i).x;
        ghcnfyssa(i).trmin = min(ghcnfyssa(i).tryears); %start year
        ghcnfyssa(i).trmax = max(ghcnfyssa(i).tryears); %end year
        
        
        %compute overlap years from min/max values
        ghcnfyssa(i).tmin = max(ghcnfyssa(i).trmin, ghcnfyssa(i).stmin);
        ghcnfyssa(i).tmax = min(ghcnfyssa(i).trmax, ghcnfyssa(i).stmax);
        ghcnfyssa(i).years = ghcnfyssa(i).tmin:ghcnfyssa(i).tmax;

        
        
        
    else
        ghcnfyssa(i).locscreen = 0;
        ghcnfyssa(i).trmin = nan; ghcnfyssa(i).trmax = nan;
        ghcnfyssa(i).stmin = nan; ghcnfyssa(i).stmax = nan;
        ghcnfyssa(i).stmean = nan;
        ghcnfyssa(i).tryears = nan; ghcnfyssa(i).trtemp = nan;
        ghcnfyssa(i).years = nan(1,100);
        
    end
    
    
end


%time check using tmin/max values

for i = 1:nt
    if ghcnfyssa(i).locscreen == 0
        nls = nls + 1;
        ghcnfyssa(i).timescreen = 0;
    else
        ghcnfyssa(i).overlap = ghcnfyssa(i).tmax - ghcnfyssa(i).tmin;
        if ghcnfyssa(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcnfyssa(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcnfyssa(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcnfyssa(i).timescreen == 2
        
        %split years
        ghcnfyssa(i).yearsp1 = [ghcnfyssa(i).tmin:1960];
        ghcnfyssa(i).yearsp2 = [1961:ghcnfyssa(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcnfyssa(i).yearsp1);
        np2 = length(ghcnfyssa(i).yearsp2);
        
        ghcnfyssa(i).stdatp1 = nan(np1,1);
        ghcnfyssa(i).stdatp2 = nan(np2,1);
        ghcnfyssa(i).trdatp1 = nan(np1,1);
        ghcnfyssa(i).trdatp2 = nan(np2,1);
        ghcnfyssa(i).stsm1 = nan(np1,1);
        ghcnfyssa(i).stsm2 = nan(np2,1);
        ghcnfyssa(i).trsm1 = nan(np1,1);
        ghcnfyssa(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcnfyssa(i).yearsp1, ghcnfyssa(i).tryears);
        tr2ind = ismember(ghcnfyssa(i).yearsp2, ghcnfyssa(i).tryears);
        st1ind = ismember(ghcnfyssa(i).yearsp1, ghcnfyssa(i).styears);
        st2ind = ismember(ghcnfyssa(i).yearsp2, ghcnfyssa(i).styears);
        
        tr1ind2 = ismember(ghcnfyssa(i).tryears, ghcnfyssa(i).yearsp1);
        tr2ind2 = ismember(ghcnfyssa(i).tryears, ghcnfyssa(i).yearsp2);
        st1ind2 = ismember(ghcnfyssa(i).styears, ghcnfyssa(i).yearsp1);
        st2ind2 = ismember(ghcnfyssa(i).styears, ghcnfyssa(i).yearsp2);
        
        %get matching data
        ghcnfyssa(i).trdatp1(tr1ind) = ghcnfyssa(i).trtemp(tr1ind2);
        ghcnfyssa(i).trdatp2(tr2ind) = ghcnfyssa(i).trtemp(tr2ind2);
        ghcnfyssa(i).stdatp1(st1ind) = ghcnfyssa(i).stmean(st1ind2);
        ghcnfyssa(i).stdatp2(st2ind) = ghcnfyssa(i).stmean(st2ind2);
        
        ghcnfyssa(i).trdat = vertcat(ghcnfyssa(i).trdatp1, ghcnfyssa(i).trdatp2);
        ghcnfyssa(i).stdat = vertcat(ghcnfyssa(i).stdatp1, ghcnfyssa(i).stdatp2);
        
        
        
        
        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(ghcnfyssa(i).trdatp1));
        ntrp2 = sum(~isnan(ghcnfyssa(i).trdatp2));
        nstp1 = sum(~isnan(ghcnfyssa(i).stdatp1));
        nstp2 = sum(~isnan(ghcnfyssa(i).stdatp2));
        
        if ntrp1>15
            p1screen = p1screen + 1;
        end
        
        if nstp1>15
            p1screen = p1screen + 1;
        end
        
        if ntrp2>10
            p2screen = p2screen + 1;
        end
        
        if nstp2>10
            p2screen = p2screen + 1;
        end
        
        ghcnfyssa(i).mindatscreen = p1screen + p2screen;
        if ghcnfyssa(i).mindatscreen < 4
            nmd = nmd + 1;
        end
    else
        ghcnfyssa(i).mindatscreen = 0;
    end
    
    ghcnfyssa(i).finalscreen = ...
        ghcnfyssa(i).timescreen + ghcnfyssa(i).locscreen + ghcnfyssa(i).mindatscreen;
    
    if ghcnfyssa(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(ghcnfyssa(i).id)
        ns = ns + 1;
    end
end




np = nt - ns;

loccheck = sprintf('%u record(s) failed loc check', nls);
timecheck1 = sprintf('%u record(s) failed max year check', nts1);
timecheck2 = sprintf('%u record(s) failed min year check', nts2);
mdcheck = sprintf('%u record(s) failed minimum data check', nmd);
finalcheck = sprintf('%u record(s) omitted', ns);
passed = sprintf('%u record(s) passed screening', np);


display(loccheck)
display(timecheck1)
display(timecheck2)
display(mdcheck)
display(finalcheck)
display(passed)



save('ghcnfyssa.mat', 'ghcnfyssa')
clear all