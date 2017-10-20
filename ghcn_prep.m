%ghcn_prep

%GHCN prep
%this script prepares treering and station data for screening
% -match treerings to closest station
% -find minimum/maximum years
%(this step takes awhile, which is why i decided to separate it
%from the rest of the screening process)
%addpath(genpath('/Users/adam/Desktop/Treerings/'))


%load GHCN data 
load('ghcn_db.mat');


%load treering db (trdb_crn)
load('tr_db_crn.mat');

%set vars
ng = length(ghcn); nt = length(trdb_crn);
r = 6370; %radius of earth in km
% nst = 4; %this is the number of stations that will be used to compute regional mean (1-4)

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
    ghcna(i).id = trdb_crn(i).id;
    ghcna(i).lat = trdb_crn(i).lat;
    ghcna(i).lon = trdb_crn(i).lon;
    ghcna(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcna(i).locscreen == 1
        
        
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
        
        ghcna(i).stmin = max(stmins(1:nst));
        ghcna(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max
        ghcna(i).styears = ghcna(i).stmin:ghcna(i).stmax;
        ny = length(ghcna(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcna(i).stmean = nan(1,ny);
         
        
        %get indices
        st1ind1 = ismember(ghcna(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcna(i).styears);
        st2ind1 = ismember(ghcna(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcna(i).styears);
        st3ind1 = ismember(ghcna(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcna(i).styears);
        st4ind1 = ismember(ghcna(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcna(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcna(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        
        %get matching MXD data
        ghcna(i).tryears = trdb_crn(i).yr;
        ghcna(i).trtemp = trdb_crn(i).x;
        ghcna(i).trmin = min(ghcna(i).tryears); %start year
        ghcna(i).trmax = max(ghcna(i).tryears); %end year
        
        
        %compute overlap years from min/max values
        ghcna(i).tmin = max(ghcna(i).trmin, ghcna(i).stmin);
        ghcna(i).tmax = min(ghcna(i).trmax, ghcna(i).stmax);
        ghcna(i).years = ghcna(i).tmin:ghcna(i).tmax;

        
        
        
    else
        ghcna(i).locscreen = 0;
        ghcna(i).trmin = nan; ghcna(i).trmax = nan;
        ghcna(i).stmin = nan; ghcna(i).stmax = nan;
        ghcna(i).stmean = nan;
        ghcna(i).tryears = nan; ghcna(i).trtemp = nan;
        ghcna(i).years = nan(1,100);
        
    end
    
    
end


%time check using tmin/max values

for i = 1:nt
    if ghcna(i).locscreen == 0
        nls = nls + 1;
        ghcna(i).timescreen = 0;
    else
        ghcna(i).overlap = ghcna(i).tmax - ghcna(i).tmin;
        if ghcna(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcna(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcna(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcna(i).timescreen == 2
        
        %split years
        ghcna(i).yearsp1 = [ghcna(i).tmin:1960];
        ghcna(i).yearsp2 = [1961:ghcna(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcna(i).yearsp1);
        np2 = length(ghcna(i).yearsp2);
        
        ghcna(i).stdatp1 = nan(np1,1);
        ghcna(i).stdatp2 = nan(np2,1);
        ghcna(i).trdatp1 = nan(np1,1);
        ghcna(i).trdatp2 = nan(np2,1);
        ghcna(i).stsm1 = nan(np1,1);
        ghcna(i).stsm2 = nan(np2,1);
        ghcna(i).trsm1 = nan(np1,1);
        ghcna(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcna(i).yearsp1, ghcna(i).tryears);
        tr2ind = ismember(ghcna(i).yearsp2, ghcna(i).tryears);
        st1ind = ismember(ghcna(i).yearsp1, ghcna(i).styears);
        st2ind = ismember(ghcna(i).yearsp2, ghcna(i).styears);
        
        tr1ind2 = ismember(ghcna(i).tryears, ghcna(i).yearsp1);
        tr2ind2 = ismember(ghcna(i).tryears, ghcna(i).yearsp2);
        st1ind2 = ismember(ghcna(i).styears, ghcna(i).yearsp1);
        st2ind2 = ismember(ghcna(i).styears, ghcna(i).yearsp2);
        
        %get matching data
        ghcna(i).trdatp1(tr1ind) = ghcna(i).trtemp(tr1ind2);
        ghcna(i).trdatp2(tr2ind) = ghcna(i).trtemp(tr2ind2);
        ghcna(i).stdatp1(st1ind) = ghcna(i).stmean(st1ind2);
        ghcna(i).stdatp2(st2ind) = ghcna(i).stmean(st2ind2);
        
        ghcna(i).trdat = vertcat(ghcna(i).trdatp1, ghcna(i).trdatp2);
        ghcna(i).stdat = vertcat(ghcna(i).stdatp1, ghcna(i).stdatp2);
        
        
        
        
        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(ghcna(i).trdatp1));
        ntrp2 = sum(~isnan(ghcna(i).trdatp2));
        nstp1 = sum(~isnan(ghcna(i).stdatp1));
        nstp2 = sum(~isnan(ghcna(i).stdatp2));
        
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
        
        ghcna(i).mindatscreen = p1screen + p2screen;
        if ghcna(i).mindatscreen < 4
            nmd = nmd + 1;
        end
    else
        ghcna(i).mindatscreen = 0;
    end
    
    ghcna(i).finalscreen = ...
        ghcna(i).timescreen + ghcna(i).locscreen + ghcna(i).mindatscreen;
    
    if ghcna(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(ghcna(i).id)
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



save('ghcna.mat', 'ghcna')
