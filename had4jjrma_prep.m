%had4jjrma_prep
addpath(genpath('/home/geovault-02/avaccaro/Treerings'))
%resize hadcrut4 dataset and prepare for intra_annual_avg (change to
%fractional time-axis)

load('hadcrut4.mat'); %load hadcrut4 data
load('hadcrut4_iridge_v1.mat'); %load hadcrut4 infilled data
load('tr_db_crn.mat'); %load MXD data
nst = 4;
nt = length(trdb_crn); %number of tree sites

r = 6370; %radius of Earth in km

%nst = 4; %this is the number of stations that will be used to compute regional mean (1-4)

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


%set check counters = 0
nls = 0; %loc check
nts1 = 0; %max year check
nts2 = 0; %min year check
nmd = 0; %minimum data check
ns = 0; %final check



%create annual time axis
dtime = double(675700 + floor(H.time)); %days since Jan 1 0 AD 00:00
dtime2 = datevec(dtime); %convert to datenum format
for j = 1:size(dtime2,1)
    dtime3(j) = dtime2(j,1) + dtime2(j,2)/12; %fractional t-axis
end


% %intra_annual_avg
for i = 1:nloc
    [da(:,i),ta(:,i)] = intra_annual_avg(Xhf(:,i),dtime3',6,7);
end







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
    had4jjrma(i).id = trdb_crn(i).id;
    had4jjrma(i).lat = trdb_crn(i).lat;
    had4jjrma(i).lon = trdb_crn(i).lon;
    had4jjrma(i).locscreen = trdb_crn(i).locscreen;
    
    if had4jjrma(i).locscreen == 1
        
        
        %compute regional mean
        
        [smin idmin] = sort(d(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of closest 4 stations
        %[~,ind(i)]=min(d(i,:));%index of closest station
        
        
        %get stations' data and years
        st1 =  da(:,minIds(1));
        st1yrs = ta(:,minIds(1));
        st2 = da(:,minIds(2));
        st2yrs = ta(:,minIds(2));
        st3 = da(:,minIds(3));
        st3yrs = ta(:,minIds(3));
        st4 = da(:,minIds(4));
        st4yrs = ta(:,minIds(4));
        
        %find mins/maxs
        st1min = min(st1yrs(isfinite(st1))); st2min = min(st2yrs(isfinite(st2)));
        st3min = min(st3yrs(isfinite(st3))); st4min = min(st4yrs(isfinite(st4)));
        st1max = max(st1yrs(isfinite(st1))); st2max = max(st2yrs(isfinite(st2)));
        st3max = max(st3yrs(isfinite(st3))); st4max = max(st4yrs(isfinite(st4)));
        stmins = [st1min, st2min, st3min, st4min];
        stmaxs = [st1max, st2max, st3max, st4max];
        
        had4jjrma(i).stmin = max(stmins(1:nst));
        had4jjrma(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max
        had4jjrma(i).styears = had4jjrma(i).stmin:had4jjrma(i).stmax;
        ny = length(had4jjrma(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        had4jjrma(i).stmean = nan(1,ny);
         
        
        %get indices
        st1ind1 = ismember(had4jjrma(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, had4jjrma(i).styears);
        st2ind1 = ismember(had4jjrma(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, had4jjrma(i).styears);
        st3ind1 = ismember(had4jjrma(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, had4jjrma(i).styears);
        st4ind1 = ismember(had4jjrma(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, had4jjrma(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            had4jjrma(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        
        %get matching MXD data
        had4jjrma(i).tryears = trdb_crn(i).yr;
        had4jjrma(i).trtemp = trdb_crn(i).x;
        had4jjrma(i).trmin = min(had4jjrma(i).tryears); %start year
        had4jjrma(i).trmax = max(had4jjrma(i).tryears); %end year
        
        
        %compute overlap years from min/max values
        had4jjrma(i).tmin = max(had4jjrma(i).trmin, had4jjrma(i).stmin);
        had4jjrma(i).tmax = min(had4jjrma(i).trmax, had4jjrma(i).stmax);
        had4jjrma(i).years = had4jjrma(i).tmin:had4jjrma(i).tmax;

        
        
        
    else
        had4jjrma(i).locscreen = 0;
        had4jjrma(i).trmin = nan; had4jjrma(i).trmax = nan;
        had4jjrma(i).stmin = nan; had4jjrma(i).stmax = nan;
        had4jjrma(i).stmean = nan;
        had4jjrma(i).tryears = nan; had4jjrma(i).trtemp = nan;
        had4jjrma(i).years = nan(1,100);
        
    end
    
    
end


%time check using tmin/max values

for i = 1:nt
    if had4jjrma(i).locscreen == 0
        nls = nls + 1;
        had4jjrma(i).timescreen = 0;
    else
        had4jjrma(i).overlap = had4jjrma(i).tmax - had4jjrma(i).tmin;
        if had4jjrma(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if had4jjrma(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    had4jjrma(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if had4jjrma(i).timescreen == 2
        
        %split years
        had4jjrma(i).yearsp1 = [had4jjrma(i).tmin:1960];
        had4jjrma(i).yearsp2 = [1961:had4jjrma(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(had4jjrma(i).yearsp1);
        np2 = length(had4jjrma(i).yearsp2);
        
        had4jjrma(i).stdatp1 = nan(np1,1);
        had4jjrma(i).stdatp2 = nan(np2,1);
        had4jjrma(i).trdatp1 = nan(np1,1);
        had4jjrma(i).trdatp2 = nan(np2,1);
        had4jjrma(i).stsm1 = nan(np1,1);
        had4jjrma(i).stsm2 = nan(np2,1);
        had4jjrma(i).trsm1 = nan(np1,1);
        had4jjrma(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(had4jjrma(i).yearsp1, had4jjrma(i).tryears);
        tr2ind = ismember(had4jjrma(i).yearsp2, had4jjrma(i).tryears);
        st1ind = ismember(had4jjrma(i).yearsp1, had4jjrma(i).styears);
        st2ind = ismember(had4jjrma(i).yearsp2, had4jjrma(i).styears);
        
        tr1ind2 = ismember(had4jjrma(i).tryears, had4jjrma(i).yearsp1);
        tr2ind2 = ismember(had4jjrma(i).tryears, had4jjrma(i).yearsp2);
        st1ind2 = ismember(had4jjrma(i).styears, had4jjrma(i).yearsp1);
        st2ind2 = ismember(had4jjrma(i).styears, had4jjrma(i).yearsp2);
        
        %get matching data
        had4jjrma(i).trdatp1(tr1ind) = had4jjrma(i).trtemp(tr1ind2);
        had4jjrma(i).trdatp2(tr2ind) = had4jjrma(i).trtemp(tr2ind2);
        had4jjrma(i).stdatp1(st1ind) = had4jjrma(i).stmean(st1ind2);
        had4jjrma(i).stdatp2(st2ind) = had4jjrma(i).stmean(st2ind2);
        
        had4jjrma(i).trdat = vertcat(had4jjrma(i).trdatp1, had4jjrma(i).trdatp2);
        had4jjrma(i).stdat = vertcat(had4jjrma(i).stdatp1, had4jjrma(i).stdatp2);
        
        
        
        
        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(had4jjrma(i).trdatp1));
        ntrp2 = sum(~isnan(had4jjrma(i).trdatp2));
        nstp1 = sum(~isnan(had4jjrma(i).stdatp1));
        nstp2 = sum(~isnan(had4jjrma(i).stdatp2));
        
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
        
        had4jjrma(i).mindatscreen = p1screen + p2screen;
        if had4jjrma(i).mindatscreen < 4
            nmd = nmd + 1;
        end
    else
        had4jjrma(i).mindatscreen = 0;
    end
    
    had4jjrma(i).finalscreen = ...
        had4jjrma(i).timescreen + had4jjrma(i).locscreen + had4jjrma(i).mindatscreen;
    
    if had4jjrma(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(had4jjrma(i).id)
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



save('had4jjrma.mat', 'had4jjrma')
clear all

