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
nst=1;
%set vars
ng = length(ghcn); nt = length(trdb_crn);
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
    ghcnjjassa(i).id = trdb_crn(i).id;
    ghcnjjassa(i).lat = trdb_crn(i).lat;
    ghcnjjassa(i).lon = trdb_crn(i).lon;
    ghcnjjassa(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcnjjassa(i).locscreen == 1
        
        
        %compute regional mean
        
        [smin idmin] = sort(d3(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of closest 4 stations
        %[~,ind(i)]=min(d(i,:));%index of closest station
        
        
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
        
        ghcnjjassa(i).stmin = max(stmins(1:nst));
        ghcnjjassa(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max
        ghcnjjassa(i).styears = ghcnjjassa(i).stmin:ghcnjjassa(i).stmax;
        ny = length(ghcnjjassa(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcnjjassa(i).stmean = nan(1,ny);
         
        
        %get indices
        st1ind1 = ismember(ghcnjjassa(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcnjjassa(i).styears);
        st2ind1 = ismember(ghcnjjassa(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcnjjassa(i).styears);
        st3ind1 = ismember(ghcnjjassa(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcnjjassa(i).styears);
        st4ind1 = ismember(ghcnjjassa(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcnjjassa(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcnjjassa(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        
        %get matching MXD data
        ghcnjjassa(i).tryears = trdb_crn(i).yr;
        ghcnjjassa(i).trtemp = trdb_crn(i).x;
        ghcnjjassa(i).trmin = min(ghcnjjassa(i).tryears); %start year
        ghcnjjassa(i).trmax = max(ghcnjjassa(i).tryears); %end year
        
        
        %compute overlap years from min/max values
        ghcnjjassa(i).tmin = max(ghcnjjassa(i).trmin, ghcnjjassa(i).stmin);
        ghcnjjassa(i).tmax = min(ghcnjjassa(i).trmax, ghcnjjassa(i).stmax);
        ghcnjjassa(i).years = ghcnjjassa(i).tmin:ghcnjjassa(i).tmax;

        
        
        
    else
        ghcnjjassa(i).locscreen = 0;
        ghcnjjassa(i).trmin = nan; ghcnjjassa(i).trmax = nan;
        ghcnjjassa(i).stmin = nan; ghcnjjassa(i).stmax = nan;
        ghcnjjassa(i).stmean = nan;
        ghcnjjassa(i).tryears = nan; ghcnjjassa(i).trtemp = nan;
        ghcnjjassa(i).years = nan(1,100);
        
    end
    
    
end


%time check using tmin/max values

for i = 1:nt
    if ghcnjjassa(i).locscreen == 0
        nls = nls + 1;
        ghcnjjassa(i).timescreen = 0;
    else
        ghcnjjassa(i).overlap = ghcnjjassa(i).tmax - ghcnjjassa(i).tmin;
        if ghcnjjassa(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcnjjassa(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcnjjassa(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcnjjassa(i).timescreen == 2
        
        %split years
        ghcnjjassa(i).yearsp1 = [ghcnjjassa(i).tmin:1960];
        ghcnjjassa(i).yearsp2 = [1961:ghcnjjassa(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcnjjassa(i).yearsp1);
        np2 = length(ghcnjjassa(i).yearsp2);
        
        ghcnjjassa(i).stdatp1 = nan(np1,1);
        ghcnjjassa(i).stdatp2 = nan(np2,1);
        ghcnjjassa(i).trdatp1 = nan(np1,1);
        ghcnjjassa(i).trdatp2 = nan(np2,1);
        ghcnjjassa(i).stsm1 = nan(np1,1);
        ghcnjjassa(i).stsm2 = nan(np2,1);
        ghcnjjassa(i).trsm1 = nan(np1,1);
        ghcnjjassa(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcnjjassa(i).yearsp1, ghcnjjassa(i).tryears);
        tr2ind = ismember(ghcnjjassa(i).yearsp2, ghcnjjassa(i).tryears);
        st1ind = ismember(ghcnjjassa(i).yearsp1, ghcnjjassa(i).styears);
        st2ind = ismember(ghcnjjassa(i).yearsp2, ghcnjjassa(i).styears);
        
        tr1ind2 = ismember(ghcnjjassa(i).tryears, ghcnjjassa(i).yearsp1);
        tr2ind2 = ismember(ghcnjjassa(i).tryears, ghcnjjassa(i).yearsp2);
        st1ind2 = ismember(ghcnjjassa(i).styears, ghcnjjassa(i).yearsp1);
        st2ind2 = ismember(ghcnjjassa(i).styears, ghcnjjassa(i).yearsp2);
        
        %get matching data
        ghcnjjassa(i).trdatp1(tr1ind) = ghcnjjassa(i).trtemp(tr1ind2);
        ghcnjjassa(i).trdatp2(tr2ind) = ghcnjjassa(i).trtemp(tr2ind2);
        ghcnjjassa(i).stdatp1(st1ind) = ghcnjjassa(i).stmean(st1ind2);
        ghcnjjassa(i).stdatp2(st2ind) = ghcnjjassa(i).stmean(st2ind2);
        
        ghcnjjassa(i).trdat = vertcat(ghcnjjassa(i).trdatp1, ghcnjjassa(i).trdatp2);
        ghcnjjassa(i).stdat = vertcat(ghcnjjassa(i).stdatp1, ghcnjjassa(i).stdatp2);
        
        
        
        
        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(ghcnjjassa(i).trdatp1));
        ntrp2 = sum(~isnan(ghcnjjassa(i).trdatp2));
        nstp1 = sum(~isnan(ghcnjjassa(i).stdatp1));
        nstp2 = sum(~isnan(ghcnjjassa(i).stdatp2));
        
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
        
        ghcnjjassa(i).mindatscreen = p1screen + p2screen;
        if ghcnjjassa(i).mindatscreen < 4
            nmd = nmd + 1;
        end
    else
        ghcnjjassa(i).mindatscreen = 0;
    end
    
    ghcnjjassa(i).finalscreen = ...
        ghcnjjassa(i).timescreen + ghcnjjassa(i).locscreen + ghcnjjassa(i).mindatscreen;
    
    if ghcnjjassa(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(ghcnjjassa(i).id)
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



save('ghcnjjassa.mat', 'ghcnjjassa')
clear all