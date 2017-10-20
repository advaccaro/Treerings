%resize hadcrut4 dataset and prepare for intra_annual_avg (change to
%fractional time-axis)

load('inst.mat'); %load hadcrut3 data
load('latlon_temp.mat');
load('tr_db_crn.mat'); %load MXD data

hadcrut3 = flipud(inst(:,2:end));
nt = length(trdb_crn); %number of tree sites
nloc = length(t_lat); %number of stations
ntime = length(inst(:,1));

r = 6370; %radius of Earth in km

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
    for j = 1:nloc %for each hadcrut loc        
        lat1 = c*trdb_crn(i).lat; %tree lat
        lon1 = c*trdb_crn(i).lon; %tree lon
        lat2 = c*t_lat(j); %st lat
        lon2 = c*t_lon(j); %st lon
        
        if trdb_crn(i).locscreen == 1
            d2(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2,r);
        else
            d2(i,j)=NaN;
        end
    end
end



for i = 1:nt;
    %store MXD meta data
    had3a(i).id = trdb_crn(i).id;
    had3a(i).lat = trdb_crn(i).lat;
    had3a(i).lon = trdb_crn(i).lon;
    had3a(i).locscreen = trdb_crn(i).locscreen;
    
    if had3a(i).locscreen == 1
        
        
        %compute regional mean
        
        [smin idmin] = sort(d2(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of closest 4 stations
        %[~,ind(i)]=min(d(i,:));%index of closest station
        
        
        %get stations' data and years
        st1 =  hadcrut3(:,minIds(1));
        st1yrs = inst(:,1);
        st2 = hadcrut3(:,minIds(2));
        st2yrs = inst(:,1);
        st3 = hadcrut3(:,minIds(3));
        st3yrs = inst(:,1);
        st4 = hadcrut3(:,minIds(4));
        st4yrs = inst(:,1);
        
        %find mins/maxs
        st1min = min(st1yrs(isfinite(st1))); st2min = min(st2yrs(isfinite(st2)));
        st3min = min(st3yrs(isfinite(st3))); st4min = min(st4yrs(isfinite(st4)));
        st1max = max(st1yrs(isfinite(st1))); st2max = max(st2yrs(isfinite(st2)));
        st3max = max(st3yrs(isfinite(st3))); st4max = max(st4yrs(isfinite(st4)));
        stmins = [st1min, st2min, st3min, st4min];
        stmaxs = [st1max, st2max, st3max, st4max];
        
        had3a(i).stmin = max(stmins(1:nst));
        had3a(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max
        had3a(i).styears = had3a(i).stmin:had3a(i).stmax;
        ny = length(had3a(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        had3a(i).stmean = nan(1,ny);
         
        
        %get indices
        st1ind1 = ismember(had3a(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, had3a(i).styears);
        st2ind1 = ismember(had3a(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, had3a(i).styears);
        st3ind1 = ismember(had3a(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, had3a(i).styears);
        st4ind1 = ismember(had3a(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, had3a(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            had3a(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        
        %get matching MXD data
        had3a(i).tryears = trdb_crn(i).yr;
        had3a(i).trtemp = trdb_crn(i).x;
        had3a(i).trmin = min(had3a(i).tryears); %start year
        had3a(i).trmax = max(had3a(i).tryears); %end year
        
        
        %compute overlap years from min/max values
        had3a(i).tmin = max(had3a(i).trmin, had3a(i).stmin);
        had3a(i).tmax = min(had3a(i).trmax, had3a(i).stmax);
        had3a(i).years = had3a(i).tmin:had3a(i).tmax;

        
        
        
    else
        had3a(i).locscreen = 0;
        had3a(i).trmin = nan; had3a(i).trmax = nan;
        had3a(i).stmin = nan; had3a(i).stmax = nan;
        had3a(i).stmean = nan;
        had3a(i).tryears = nan; had3a(i).trtemp = nan;
        had3a(i).years = nan(1,100);
        
    end
    
    
end


%time check using tmin/max values

for i = 1:nt
    if had3a(i).locscreen == 0
        nls = nls + 1;
        had3a(i).timescreen = 0;
    else
        had3a(i).overlap = had3a(i).tmax - had3a(i).tmin;
        if had3a(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if had3a(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    had3a(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if had3a(i).timescreen == 2
        
        %split years
        had3a(i).yearsp1 = [had3a(i).tmin:1960];
        had3a(i).yearsp2 = [1961:had3a(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(had3a(i).yearsp1);
        np2 = length(had3a(i).yearsp2);
        
        had3a(i).stdatp1 = nan(np1,1);
        had3a(i).stdatp2 = nan(np2,1);
        had3a(i).trdatp1 = nan(np1,1);
        had3a(i).trdatp2 = nan(np2,1);
        had3a(i).stsm1 = nan(np1,1);
        had3a(i).stsm2 = nan(np2,1);
        had3a(i).trsm1 = nan(np1,1);
        had3a(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(had3a(i).yearsp1, had3a(i).tryears);
        tr2ind = ismember(had3a(i).yearsp2, had3a(i).tryears);
        st1ind = ismember(had3a(i).yearsp1, had3a(i).styears);
        st2ind = ismember(had3a(i).yearsp2, had3a(i).styears);
        
        tr1ind2 = ismember(had3a(i).tryears, had3a(i).yearsp1);
        tr2ind2 = ismember(had3a(i).tryears, had3a(i).yearsp2);
        st1ind2 = ismember(had3a(i).styears, had3a(i).yearsp1);
        st2ind2 = ismember(had3a(i).styears, had3a(i).yearsp2);
        
        %get matching data
        had3a(i).trdatp1(tr1ind) = had3a(i).trtemp(tr1ind2);
        had3a(i).trdatp2(tr2ind) = had3a(i).trtemp(tr2ind2);
        had3a(i).stdatp1(st1ind) = had3a(i).stmean(st1ind2);
        had3a(i).stdatp2(st2ind) = had3a(i).stmean(st2ind2);
        
        had3a(i).trdat = vertcat(had3a(i).trdatp1, had3a(i).trdatp2);
        had3a(i).stdat = vertcat(had3a(i).stdatp1, had3a(i).stdatp2);
        
        
        
        
        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(had3a(i).trdatp1));
        ntrp2 = sum(~isnan(had3a(i).trdatp2));
        nstp1 = sum(~isnan(had3a(i).stdatp1));
        nstp2 = sum(~isnan(had3a(i).stdatp2));
        
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
        
        had3a(i).mindatscreen = p1screen + p2screen;
        if had3a(i).mindatscreen < 4
            nmd = nmd + 1;
        end
    else
        had3a(i).mindatscreen = 0;
    end
    
    had3a(i).finalscreen = ...
        had3a(i).timescreen + had3a(i).locscreen + had3a(i).mindatscreen;
    
    if had3a(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(had3a(i).id)
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



save('had3a.mat', 'had3a')
