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
    ghcngsa(i).id = trdb_crn(i).id;
    ghcngsa(i).lat = trdb_crn(i).lat;
    ghcngsa(i).lon = trdb_crn(i).lon;
    ghcngsa(i).locscreen = trdb_crn(i).locscreen;
    
    if ghcngsa(i).locscreen == 1
        
        
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
        
        ghcngsa(i).stmin = max(stmins(1:nst));
        ghcngsa(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max
        ghcngsa(i).styears = ghcngsa(i).stmin:ghcngsa(i).stmax;
        ny = length(ghcngsa(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        ghcngsa(i).stmean = nan(1,ny);
         
        
        %get indices
        st1ind1 = ismember(ghcngsa(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, ghcngsa(i).styears);
        st2ind1 = ismember(ghcngsa(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, ghcngsa(i).styears);
        st3ind1 = ismember(ghcngsa(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, ghcngsa(i).styears);
        st4ind1 = ismember(ghcngsa(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, ghcngsa(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            ghcngsa(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        
        %get matching MXD data
        ghcngsa(i).tryears = trdb_crn(i).yr;
        ghcngsa(i).trtemp = trdb_crn(i).x;
        ghcngsa(i).trmin = min(ghcngsa(i).tryears); %start year
        ghcngsa(i).trmax = max(ghcngsa(i).tryears); %end year
        
        
        %compute overlap years from min/max values
        ghcngsa(i).tmin = max(ghcngsa(i).trmin, ghcngsa(i).stmin);
        ghcngsa(i).tmax = min(ghcngsa(i).trmax, ghcngsa(i).stmax);
        ghcngsa(i).years = ghcngsa(i).tmin:ghcngsa(i).tmax;

        
        
        
    else
        ghcngsa(i).locscreen = 0;
        ghcngsa(i).trmin = nan; ghcngsa(i).trmax = nan;
        ghcngsa(i).stmin = nan; ghcngsa(i).stmax = nan;
        ghcngsa(i).stmean = nan;
        ghcngsa(i).tryears = nan; ghcngsa(i).trtemp = nan;
        ghcngsa(i).years = nan(1,100);
        
    end
    
    
end


%time check using tmin/max values

for i = 1:nt
    if ghcngsa(i).locscreen == 0
        nls = nls + 1;
        ghcngsa(i).timescreen = 0;
    else
        ghcngsa(i).overlap = ghcngsa(i).tmax - ghcngsa(i).tmin;
        if ghcngsa(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if ghcngsa(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    ghcngsa(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if ghcngsa(i).timescreen == 2
        
        %split years
        ghcngsa(i).yearsp1 = [ghcngsa(i).tmin:1960];
        ghcngsa(i).yearsp2 = [1961:ghcngsa(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(ghcngsa(i).yearsp1);
        np2 = length(ghcngsa(i).yearsp2);
        
        ghcngsa(i).stdatp1 = nan(np1,1);
        ghcngsa(i).stdatp2 = nan(np2,1);
        ghcngsa(i).trdatp1 = nan(np1,1);
        ghcngsa(i).trdatp2 = nan(np2,1);
        ghcngsa(i).stsm1 = nan(np1,1);
        ghcngsa(i).stsm2 = nan(np2,1);
        ghcngsa(i).trsm1 = nan(np1,1);
        ghcngsa(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(ghcngsa(i).yearsp1, ghcngsa(i).tryears);
        tr2ind = ismember(ghcngsa(i).yearsp2, ghcngsa(i).tryears);
        st1ind = ismember(ghcngsa(i).yearsp1, ghcngsa(i).styears);
        st2ind = ismember(ghcngsa(i).yearsp2, ghcngsa(i).styears);
        
        tr1ind2 = ismember(ghcngsa(i).tryears, ghcngsa(i).yearsp1);
        tr2ind2 = ismember(ghcngsa(i).tryears, ghcngsa(i).yearsp2);
        st1ind2 = ismember(ghcngsa(i).styears, ghcngsa(i).yearsp1);
        st2ind2 = ismember(ghcngsa(i).styears, ghcngsa(i).yearsp2);
        
        %get matching data
        ghcngsa(i).trdatp1(tr1ind) = ghcngsa(i).trtemp(tr1ind2);
        ghcngsa(i).trdatp2(tr2ind) = ghcngsa(i).trtemp(tr2ind2);
        ghcngsa(i).stdatp1(st1ind) = ghcngsa(i).stmean(st1ind2);
        ghcngsa(i).stdatp2(st2ind) = ghcngsa(i).stmean(st2ind2);
        
        ghcngsa(i).trdat = vertcat(ghcngsa(i).trdatp1, ghcngsa(i).trdatp2);
        ghcngsa(i).stdat = vertcat(ghcngsa(i).stdatp1, ghcngsa(i).stdatp2);
        
        
        
        
        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(ghcngsa(i).trdatp1));
        ntrp2 = sum(~isnan(ghcngsa(i).trdatp2));
        nstp1 = sum(~isnan(ghcngsa(i).stdatp1));
        nstp2 = sum(~isnan(ghcngsa(i).stdatp2));
        
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
        
        ghcngsa(i).mindatscreen = p1screen + p2screen;
        if ghcngsa(i).mindatscreen < 4
            nmd = nmd + 1;
        end
    else
        ghcngsa(i).mindatscreen = 0;
    end
    
    ghcngsa(i).finalscreen = ...
        ghcngsa(i).timescreen + ghcngsa(i).locscreen + ghcngsa(i).mindatscreen;
    
    if ghcngsa(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(ghcngsa(i).id)
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



save('ghcngsa.mat', 'ghcngsa')
