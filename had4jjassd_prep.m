%hadcrut4infgs_smprep
%resize hadcrut4 dataset and prepare for intra_annual_avg (change to
%fractional time-axis)
addpath(genpath('/home/geovault-02/avaccaro/Treerings'))
load('hadcrut4_iridge_v1.mat'); %load hadcrut4 infilled data
load('hadcrut4.mat'); %load hadcrut4 data
load('tr_db_crn.mat'); %load MXD data

nt = length(trdb_crn); %number of tree sites

r = 6370; %radius of Earth in km
nst = 1;
%nst = 4; %this is the number of stations that will be used to compute regional mean (1-4)


[ntime, nloc] = size(Xhf);
%nloc = nlat * nlon; %number of stations
nyears = floor(ntime/12);
hyears = 1850:2012; %hadcrut4a years

%hadcrut4 = reshape(H.d, [nloc, ntime]);
%hadcrut4 = hadcrut4'; % time x space
%hadcrut4 = reshape(H.d, [ntime nloc]);
%hadcrut4(hadcrut4<-1000) = NaN;

[x,y] = meshgrid(H.lon,H.lat);
tloc = [x(:),y(:)];
loc = tloc(station,:);



%create annual time axis
dtime = double(675700 + floor(H.time)); %days since 00:00 Jan 1 AD 0 
dtime2 = datevec(dtime); %convert to datenum format
for j = 1:size(dtime2,1)
    dtime3(j) = dtime2(j,1) + dtime2(j,2)/12; %fractional t-axis
end


%%intra_annual_avg
for i = 1:nloc
    [da2(:,i),ta2(:,i)]=intra_annual_avg(Xhf(:,i),dtime3',6,8); 
end


c = pi/180;

%find closest station and create tree-ring/station pair
for i = 1:nt %for each site
    for j = 1:nloc %for each hadcrut loc
        
        lat1 = c*trdb_crn(i).lat; %tree lat
        lon1 = c*trdb_crn(i).lon; %tree lon
        lat2 = c*loc(j,2); %st lat
        lon2 = c*loc(j,1); %st lon
        
        if trdb_crn(i).locscreen == 1
            d(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2, r);
        else
            d(i,j)=NaN;
        end
    end
end


for i = 1:nt;
    
    %store MXD meta data
    had4jjassd(i).id = trdb_crn(i).id;
    had4jjassd(i).lat = trdb_crn(i).lat;
    had4jjassd(i).lon = trdb_crn(i).lon;
    had4jjassd(i).locscreen = trdb_crn(i).locscreen;
    
    if had4jjassd(i).locscreen == 1
        
        %compute regional means

        [smin idmin] = sort(d(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of lowest 4
        
        %get stations' data and years
        st1 =  da2(:,minIds(1));
        st1yrs = ta2(:,minIds(1));
        st2 = da2(:,minIds(2));
        st2yrs = ta2(:,minIds(2));
        st3 = da2(:,minIds(3));
        st3yrs = ta2(:,minIds(3));
        st4 = da2(:,minIds(4));
        st4yrs = ta2(:,minIds(4));
        
        %find mins/maxs
        st1min = min(st1yrs(isfinite(st1))); st2min = min(st2yrs(isfinite(st2)));
        st3min = min(st3yrs(isfinite(st3))); st4min = min(st4yrs(isfinite(st4)));
        st1max = max(st1yrs(isfinite(st1))); st2max = max(st2yrs(isfinite(st2)));
        st3max = max(st3yrs(isfinite(st3))); st4max = max(st4yrs(isfinite(st4)));
        stmins = [st1min, st2min, st3min, st4min];
        stmaxs = [st1max, st2max, st3max, st4max];
        
        had4jjassd(i).stmin = max(stmins(1:nst));
        had4jjassd(i).stmax = min(stmaxs(1:nst));
        
        %st years from min/max (overlap)
        had4jjassd(i).styears = had4jjassd(i).stmin:had4jjassd(i).stmax;
        ny = length(had4jjassd(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        had4jjassd(i).stmean = nan(1,ny);
        
        
        %get indices
        st1ind1 = ismember(had4jjassd(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, had4jjassd(i).styears);
        st2ind1 = ismember(had4jjassd(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, had4jjassd(i).styears);
        st3ind1 = ismember(had4jjassd(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, had4jjassd(i).styears);
        st4ind1 = ismember(had4jjassd(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, had4jjassd(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            had4jjassd(i).stmean(j) = nmean(A(1:nst,j));
        end
        
        %get matching MXD data
        had4jjassd(i).tryears = trdb_crn(i).yr;
        had4jjassd(i).trtemp = trdb_crn(i).x;
        had4jjassd(i).trmin = min(had4jjassd(i).tryears);
        had4jjassd(i).trmax = max(had4jjassd(i).tryears);
        
        %set overlap year range from min/max values
        had4jjassd(i).tmin = max(had4jjassd(i).trmin, had4jjassd(i).stmin);
        had4jjassd(i).tmax = min(had4jjassd(i).trmax, had4jjassd(i).stmax);
        had4jjassd(i).years = had4jjassd(i).tmin:had4jjassd(i).tmax;
  
        
    else
        had4jjassd(i).trmin = nan; had4jjassd(i).trmax = nan;
        had4jjassd(i).tryears = trdb_crn(i).yr;
        had4jjassd(i).stmin = nan; had4jjassd(i).stmax = nan;
        had4jjassd(i).stmean = nan;
        had4jjassd(i).years = nan(1,100);
        had4jjassd(i).locscreen = 0;
        had4jjassd(i).tmin = nan; had4jjassd(i).tmax = nan;
        had4jjassd(i).years = nan;
    end
    
    
end





%preallocate for decadal smoothing
for i = 1:nt
    if had4jjassd(i).locscreen == 1
        had4jjassd(i).stsmooth = nan(length(had4jjassd(i).styears),1);
        had4jjassd(i).trsmooth = nan(length(had4jjassd(i).tryears),1);
    end
end

for i = 1:nt
    if had4jjassd(i).locscreen == 1
        nntr = sum(~isnan(had4jjassd(i).trtemp));
        nnst = sum(~isnan(had4jjassd(i).stmean));
        if (nntr > 20 & nnst > 20)
            had4jjassd(i).mindatsm = 1;
            
            %decadally smoothed
            indmst = isfinite(had4jjassd(i).stmean);
            indmtr = isfinite(had4jjassd(i).trtemp);
            had4jjassd(i).stsmooth(indmst) = hepta_smooth(had4jjassd(i).stmean(indmst), 1/10);
            had4jjassd(i).trsmooth(indmtr) = hepta_smooth(had4jjassd(i).trtemp(indmtr), 1/10);
            
        else
            had4jjassd(i).mindatsm = 0;
        end
    else
        had4jjassd(i).mindatsm = 0;
    end
end


%time check using tmin/max values
nls = 0;
nts1 = 0;
nts2 = 0;
for i = 1:nt
    if had4jjassd(i).locscreen == 0
        nls = nls + 1;
        had4jjassd(i).timescreen = 0;
    else
        had4jjassd(i).overlap = had4jjassd(i).tmax - had4jjassd(i).tmin;
        if had4jjassd(i).tmax < 1960
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if had4jjassd(i).tmin > 1950
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    had4jjassd(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1960
    if had4jjassd(i).timescreen == 2
        
        %split years
        had4jjassd(i).yearsp1 = [had4jjassd(i).tmin:1960];
        had4jjassd(i).yearsp2 = [1961:had4jjassd(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(had4jjassd(i).yearsp1);
        np2 = length(had4jjassd(i).yearsp2);
        

        had4jjassd(i).stsm1 = nan(np1,1);
        had4jjassd(i).stsm2 = nan(np2,1);
        had4jjassd(i).trsm1 = nan(np1,1);
        had4jjassd(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(had4jjassd(i).yearsp1, had4jjassd(i).tryears);
        tr2ind = ismember(had4jjassd(i).yearsp2, had4jjassd(i).tryears);
        st1ind = ismember(had4jjassd(i).yearsp1, had4jjassd(i).styears);
        st2ind = ismember(had4jjassd(i).yearsp2, had4jjassd(i).styears);
        
        tr1ind2 = ismember(had4jjassd(i).tryears, had4jjassd(i).yearsp1);
        tr2ind2 = ismember(had4jjassd(i).tryears, had4jjassd(i).yearsp2);
        st1ind2 = ismember(had4jjassd(i).styears, had4jjassd(i).yearsp1);
        st2ind2 = ismember(had4jjassd(i).styears, had4jjassd(i).yearsp2);
        
        %get matching data
        had4jjassd(i).trsm1(tr1ind) = had4jjassd(i).trsmooth(tr1ind2);
        had4jjassd(i).trsm2(tr2ind) = had4jjassd(i).trsmooth(tr2ind2);
        had4jjassd(i).stsm1(st1ind) = had4jjassd(i).stsmooth(st1ind2);
        had4jjassd(i).stsm2(st2ind) = had4jjassd(i).stsmooth(st2ind2);
        
        had4jjassd(i).trsm = vertcat(had4jjassd(i).trsm1, had4jjassd(i).trsm2);
        had4jjassd(i).stsm = vertcat(had4jjassd(i).stsm1, had4jjassd(i).stsm2);
    
    else
        had4jjassd(i).mindatscreen = 0;
    end
    
end





save('had4jjassd.mat', 'had4jjassd')
clear all