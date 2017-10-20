%resize hadcrut4 dataset and prepare for intra_annual_avg (change to
%fractional time-axis)

load('hadcrut4.mat'); %load hadcrut4 data
load('tr_db_crn.mat'); %load MXD data

nt = length(trdb_crn); %number of tree sites

r = 6370; %radius of Earth in km


[nlon, nlat, ntime] = size(H.d);
nloc = nlat * nlon; %number of stations

hadcrut4 = reshape(H.d, [nloc, ntime]);
hadcrut4(hadcrut4<-1000) = NaN;

%SHOULD PROBABLY INTRODUCE LOC HERE

nyears = round(ntime/12);

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

    [da(:,i),ta(:,i)]=intra_annual_avg(hadcrut4(i,:),dtime3',6,8);

end



% H.tfrac(1) = 1850; %or does it start at 1850+1/12?
% 
% for i = 2:ntime
%     H.tfrac(i) = H.tfrac(i-1)+(1/12);
% end

[x,y] = meshgrid(H.lon,H.lat);
loc = [x(:),y(:)];




c = pi/180;

%tic;
%find closest station and create tree-ring/station pair
for i = 1:nt %for each site
    for j = 1:nloc %for each hadcrut loc
            
            lat1 = c*trdb_crn(i).lat;
            lon1 = c*trdb_crn(i).lon;
            lat2 = c*loc(j,2); %%THIS IS WRONG
            lon2 = c*loc(j,1); %%WHERE DID LOC COME FROM!?
            
            if trdb_crn(i).locscreen == 1
                d(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2);
            else
                d(i,j)=NaN;
            end
        end
    end

%toc;

%compute regional means
for i = 1:nt;
    if trdb_crn(i).locscreen == 1


        [smin idmin] = sort(d(i,:)); %sort d values
        minVals = smin(1:4); %lowest 4 d values
        minIds = idmin(1:4); %indices of lowest 4
        
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
        
        hcurve(i).stmin = max(stmins);
        hcurve(i).stmax = min(stmaxs);
        
        %st years from min/max
        hcurve(i).styears = hcurve(i).stmin:hcurve(i).stmax;
        ny = length(hcurve(i).styears);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);st4da = nan(1,ny);
        stdas = nan(4,ny);
        hcurve(i).stmean = nan(1,ny);
        
        
        %get indices
        st1ind1 = ismember(hcurve(i).styears, st1yrs);
        st1ind2 = ismember(st1yrs, hcurve(i).styears);
        st2ind1 = ismember(hcurve(i).styears, st2yrs);
        st2ind2 = ismember(st2yrs, hcurve(i).styears);
        st3ind1 = ismember(hcurve(i).styears, st3yrs);
        st3ind2 = ismember(st3yrs, hcurve(i).styears);
        st4ind1 = ismember(hcurve(i).styears, st4yrs);
        st4ind2 = ismember(st4yrs, hcurve(i).styears);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        st4da(st4ind1) = st4(st4ind2);
        
        A = vertcat(st1da, st2da, st3da, st4da);
        
        %compute stmean
        for j = 1:ny
            hcurve(i).stmean(j) = nmean(A(:,j));
        end
        
        else
        hcurve(i).trmin = nan; hcurve(i).trmax = nan;
        hcurve(i).stmin = nan; hcurve(i).stmax = nan;
        hcurve(i).stmean = nan;
        hcurve(i).years = nan(1,100);
    end


end



%find indices for MXD/station pair
for i = 1:nt
    if trdb_crn(i).locscreen == 1
        [~,ind(i)]=min(d(i,:));%index of closest station
        
        %compile matching MXD/station data
        %hcurve(i).styears = ta(:,1); %1850-2012
        hcurve(i).tryears = trdb_crn(i).yr;
        
        %hcurve(i).sttemp = da(:,ind(i));
        hcurve(i).trtemp = trdb_crn(i).x;
        
    
    else
        hcurve(i).tryears = trdb_crn(i).yr;
        na = length(hcurve(i).tryears);
        hcurve(i).tryears = nan(1,na);
        hcurve(i).trtemp = trdb_crn(i).x;
        nb = length(hcurve(i).trtemp);
        hcurve(i).sttemp = nan(1, nb);
    end
    
    
end

%store tr site metadata in hcurve
for i = 1:nt
    hcurve(i).id = trdb_crn(i).id;
    hcurve(i).lat = trdb_crn(i).lat;
    hcurve(i).lon = trdb_crn(i).lon;
    hcurve(i).locscreen = trdb_crn(i).locscreen;
    
    if hcurve(i).locscreen == 1
        
        %set hcurve year range by finding min/max values
        %hcurve(i).stmin = min(hcurve(i).styears);
        %hcurve(i).stmax = max(hcurve(i).styears);
        hcurve(i).trmin = min(hcurve(i).tryears);
        hcurve(i).trmax = max(hcurve(i).tryears);
        
        hcurve(i).tmin = max(hcurve(i).stmin, hcurve(i).trmin);
        hcurve(i).tmax = min(hcurve(i).stmax, hcurve(i).trmax);
        
        hcurve(i).years = hcurve(i).tmin:hcurve(i).tmax;
        
    else
        hcurve(i).trmin = nan; hcurve(i).trmax = nan;
        hcurve(i).stmin = nan; hcurve(i).stmax = nan;
        hcurve(i).years = nan(1,100);
    end
end
    

%time check using tmin/max values

for i = 1:nt
    if hcurve(i).locscreen == 0
        nls = nls + 1;
        hcurve(i).timescreen = 0;
    else
        hcurve(i).overlap = hcurve(i).tmax - hcurve(i).tmin;
        if hcurve(i).tmax < 1970
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if hcurve(i).tmin > 1935
            screen2 = 0;
            nts2 = nts2 +1;
        else
            screen2 = 1;
        end
    end
    
    hcurve(i).timescreen = screen1 + screen2;
    
    %split suitable records into pre-/post-1950
    if hcurve(i).timescreen == 2
        
        %split years
        hcurve(i).yearsp1 = [hcurve(i).tmin:1959];
        hcurve(i).yearsp2 = [1960:hcurve(i).tmax];
        
        %preallocate for st/tr years
        np1 = length(hcurve(i).yearsp1);
        np2 = length(hcurve(i).yearsp2);
        
        hcurve(i).stdatp1 = nan(np1,1);
        hcurve(i).stdatp2 = nan(np2,1);
        hcurve(i).trdatp1 = nan(np1,1);
        hcurve(i).trdatp2 = nan(np2,1);
        hcurve(i).stsm1 = nan(np1,1);
        hcurve(i).stsm2 = nan(np2,1);
        hcurve(i).trsm1 = nan(np1,1);
        hcurve(i).trsm2 = nan(np2,1);
        
        %find indices for each time period
        tr1ind = ismember(hcurve(i).yearsp1, hcurve(i).tryears);
        tr2ind = ismember(hcurve(i).yearsp2, hcurve(i).tryears);
        st1ind = ismember(hcurve(i).yearsp1, hcurve(i).styears);
        st2ind = ismember(hcurve(i).yearsp2, hcurve(i).styears);
        
        tr1ind2 = ismember(hcurve(i).tryears, hcurve(i).yearsp1);
        tr2ind2 = ismember(hcurve(i).tryears, hcurve(i).yearsp2);
        st1ind2 = ismember(hcurve(i).styears, hcurve(i).yearsp1);
        st2ind2 = ismember(hcurve(i).styears, hcurve(i).yearsp2);
        
        %get matching data
        hcurve(i).trdatp1(tr1ind) = hcurve(i).trtemp(tr1ind2);
        hcurve(i).trdatp2(tr2ind) = hcurve(i).trtemp(tr2ind2);
%         hcurve(i).stdatp1(st1ind) = hcurve(i).sttemp(st1ind2);
        hcurve(i).stdatp1(st1ind) = hcurve(i).stmean(st1ind2);
%         hcurve(i).stdatp2(st2ind) = hcurve(i).sttemp(st2ind2);
        hcurve(i).stdatp2(st2ind) = hcurve(i).stmean(st2ind2);
        
        hcurve(i).trdat = vertcat(hcurve(i).trdatp1, hcurve(i).trdatp2);
        hcurve(i).stdat = vertcat(hcurve(i).stdatp1, hcurve(i).stdatp2);

        
        %decadally smoothed
%         indmst = isfinite(hcurve(i).stdat);
%         indmtr = isfinite(hcurve(i).trdat);
%         hcurve(i).stsm = hepta_smooth(hcurve(i).stdat(indmst), 10);
%         hcurve(i).trsm = hepta_smooth(hcurve(i).trdat(indmtr), 10);
%         hcurve(i).stsm1(st1ind) = hcurve(i).stsm(st1ind2);
%         hcurve(i).stsm2(st2ind) = hcurve(i).stsm(st2ind2);
%         hcurve(i).trsm1(tr1ind) = hcurve(i).trsm(tr1ind2);
%         hcurve(i).trsm2(tr2ind) = hcurve(i).trsm(tr2ind2);


        %minimum data check p1
        p1screen = 0;
        p2screen = 0;
        
        ntrp1 = sum(~isnan(hcurve(i).trdatp1));
        ntrp2 = sum(~isnan(hcurve(i).trdatp2));
        nstp1 = sum(~isnan(hcurve(i).stdatp1));
        nstp2 = sum(~isnan(hcurve(i).stdatp2));
        
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
        
        hcurve(i).mindatscreen = p1screen + p2screen;
        if hcurve(i).mindatscreen < 4
            nmd = nmd + 1;
        end
     else
         hcurve(i).mindatscreen = 0;
     end
    
    hcurve(i).finalscreen = ...
        hcurve(i).timescreen + hcurve(i).locscreen + hcurve(i).mindatscreen;
    
    if hcurve(i).finalscreen < 7
        display('record screened out')
        display(i)
        display(hcurve(i).id)
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


    
save('hcurve.mat', 'hcurve')
