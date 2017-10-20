%hadcrut4inf prep
%resize hadcrut4 dataset and prepare for intra_annual_avg (change to
%fractional time-axis)
addpath(genpath('/home/geovault-02/avaccaro/Treerings/'))

load('hadcrut4_iridge_v1.mat'); %load hadcrut4 infilled data
load('hadcrut4.mat'); %load raw hadcrut4 data
load('tr_db_crn.mat'); %load MXD data

nt = length(trdb_crn); %number of tree sites

r = 6370; %radius of Earth in km
nst = 1;
% nst = 4; %this is the number of stations that will be used to compute regional mean (1-4)

[ntime, nloc] = size(Xhf);
%nloc = nlat * nlon; %number of stations
nyears = floor(ntime/12);
hyears = 1850:2012; %hadcrut4a years

%hadcrut4 = reshape(H.d, [nloc, ntime]);
%hadcrut4 = hadcrut4'; % time x space
%hadcrut4 = reshape(H.d, [ntime nloc]);
%hadcrut4(hadcrut4<-1000) = NaN;

[y,x] = meshgrid(H.lat,H.lon);
tloc = [x(:),y(:)];
loc = tloc(station,:);
%weights = cosd(loc);


%preallocate w/ nan
hadcrut4infa = nan(nyears,nloc);

%annualize hadcrut4 (via 12month mean, replaces intra_annual_avg)
for i = 1:nloc
    for j = 1:nyears
        hadcrut4infa(j,i) = nmean(Xhf(1+(j-1)*12:12*j,i));
    end
end


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
        
        lat1 = c*trdb_crn(i).lat;
        lon1 = c*trdb_crn(i).lon;
        lat2 = c*loc(j,2);
        lon2 = c*loc(j,1);
        
        if trdb_crn(i).locscreen == 1
            d(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2,r);
        else
            d(i,j)=NaN;
        end
    end
end


for i = 1:nt;
    %store MXD meta data
    had4fyssa(i).id = trdb_crn(i).id;
    had4fyssa(i).lat = trdb_crn(i).lat;
    had4fyssa(i).lon = trdb_crn(i).lon;
    had4fyssa(i).locscreen = trdb_crn(i).locscreen;
    
    if had4fyssa(i).locscreen == 1
        
        
        %compute regional mean
        
        [smin idmin] = sort(d(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of closest 4 stations
        %[~,ind(i)]=min(d(i,:));%index of closest station
        
        
        %get stations' data and years
        st1 =  hadcrut4infa(:,minIds(1));
        st1yrs = hyears;
        st2 = hadcrut4infa(:,minIds(2));
        st2yrs = hyears;
        st3 = hadcrut4infa(:,minIds(3));
        st3yrs = hyears;
        st4 = hadcrut4infa(:,minIds(4));
        st4yrs = hyears;
        
        %find mins/maxs
        st1min = min(st1yrs(isfinite(st1))); st2min = min(st2yrs(isfinite(st2)));
        st3min = min(st3yrs(isfinite(st3))); st4min = min(st4yrs(isfinite(st4)));
        st1max = max(st1yrs(isfinite(st1))); st2max = max(st2yrs(isfinite(st2)));
        st3max = max(st3yrs(isfinite(st3))); st4max = max(st4yrs(isfinite(st4)));
        stmins = [st1min, st2min, st3min, st4min];
        stmaxs = [st1max, st2max, st3max, st4max];
        
        had4fyssa(i).stmin = max(stmins(1:nst));
        had4fyssa(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max
        had4fyssa(i).styears = had4fyssa(i).stmin:had4fyssa(i).stmax;
        ny = length(had4fyssa(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        had4fyssa(i).stmean = nan(1,ny);
         
        
        %get indices
        st1ind1 = ismember(had4fyssa(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, had4fyssa(i).styears);
        st2ind1 = ismember(had4fyssa(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, had4fyssa(i).styears);
        st3ind1 = ismember(had4fyssa(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, had4fyssa(i).styears);
        st4ind1 = ismember(had4fyssa(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, had4fyssa(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            had4fyssa(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        
        %get matching MXD data
        had4fyssa(i).tryears = trdb_crn(i).yr;
        had4fyssa(i).trtemp = trdb_crn(i).x;
        had4fyssa(i).trmin = min(had4fyssa(i).tryears); %start year
        had4fyssa(i).trmax = max(had4fyssa(i).tryears); %end year
        
        
        %compute overlap years from min/max values
        had4fyssa(i).tmin = max(had4fyssa(i).trmin, had4fyssa(i).stmin);
        had4fyssa(i).tmax = min(had4fyssa(i).trmax, had4fyssa(i).stmax);
        had4fyssa(i).years = had4fyssa(i).tmin:had4fyssa(i).tmax;

        
        
        
    else
        had4fyssa(i).locscreen = 0;
        had4fyssa(i).trmin = nan; had4fyssa(i).trmax = nan;
        had4fyssa(i).stmin = nan; had4fyssa(i).stmax = nan;
        had4fyssa(i).stmean = nan;
        had4fyssa(i).tryears = nan; had4fyssa(i).trtemp = nan;
        had4fyssa(i).years = nan(1,100);
        
    end
    
    
end


%time check using tmin/max values

for i = 1:nt
    if had4fyssa(i).locscreen == 0
        nls = nls + 1;
        had4fyssa(i).timescreen = 0;
    else
        had4fyssa(i).overlap = had4fyssa(i).tmax - had4fyssa(i).tmin;
        if had4fyssa(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if had4fyssa(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    had4fyssa(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if had4fyssa(i).timescreen == 2
        
        %split years
        had4fyssa(i).yearsp1 = [had4fyssa(i).tmin:1960];
        had4fyssa(i).yearsp2 = [1961:had4fyssa(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(had4fyssa(i).yearsp1);
        np2 = length(had4fyssa(i).yearsp2);
        
        had4fyssa(i).stdatp1 = nan(np1,1);
        had4fyssa(i).stdatp2 = nan(np2,1);
        had4fyssa(i).trdatp1 = nan(np1,1);
        had4fyssa(i).trdatp2 = nan(np2,1);
        had4fyssa(i).stsm1 = nan(np1,1);
        had4fyssa(i).stsm2 = nan(np2,1);
        had4fyssa(i).trsm1 = nan(np1,1);
        had4fyssa(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(had4fyssa(i).yearsp1, had4fyssa(i).tryears);
        tr2ind = ismember(had4fyssa(i).yearsp2, had4fyssa(i).tryears);
        st1ind = ismember(had4fyssa(i).yearsp1, had4fyssa(i).styears);
        st2ind = ismember(had4fyssa(i).yearsp2, had4fyssa(i).styears);
        
        tr1ind2 = ismember(had4fyssa(i).tryears, had4fyssa(i).yearsp1);
        tr2ind2 = ismember(had4fyssa(i).tryears, had4fyssa(i).yearsp2);
        st1ind2 = ismember(had4fyssa(i).styears, had4fyssa(i).yearsp1);
        st2ind2 = ismember(had4fyssa(i).styears, had4fyssa(i).yearsp2);
        
        %get matching data
        had4fyssa(i).trdatp1(tr1ind) = had4fyssa(i).trtemp(tr1ind2);
        had4fyssa(i).trdatp2(tr2ind) = had4fyssa(i).trtemp(tr2ind2);
        had4fyssa(i).stdatp1(st1ind) = had4fyssa(i).stmean(st1ind2);
        had4fyssa(i).stdatp2(st2ind) = had4fyssa(i).stmean(st2ind2);
        
        had4fyssa(i).trdat = vertcat(had4fyssa(i).trdatp1, had4fyssa(i).trdatp2);
        had4fyssa(i).stdat = vertcat(had4fyssa(i).stdatp1, had4fyssa(i).stdatp2);
        
        
        
        
        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(had4fyssa(i).trdatp1));
        ntrp2 = sum(~isnan(had4fyssa(i).trdatp2));
        nstp1 = sum(~isnan(had4fyssa(i).stdatp1));
        nstp2 = sum(~isnan(had4fyssa(i).stdatp2));
        
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
        
        had4fyssa(i).mindatscreen = p1screen + p2screen;
        if had4fyssa(i).mindatscreen < 4
            nmd = nmd + 1;
        end
    else
        had4fyssa(i).mindatscreen = 0;
    end
    
    had4fyssa(i).finalscreen = ...
        had4fyssa(i).timescreen + had4fyssa(i).locscreen + had4fyssa(i).mindatscreen;
    
    if had4fyssa(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(had4fyssa(i).id)
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



save('had4fyssa.mat', 'had4fyssa')
clear all
