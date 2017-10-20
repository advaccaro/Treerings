%compare hadcrut4 to hadcrut3
load('latlon_temp.mat') %t_lat, t_lon
load('inst.mat') %inst
load('hadcrut4.mat') %hadcrut4

years = 1850:2006;

r = 6370; %radius of Earth in km


[nlon, nlat, ntime] = size(H.d);
nloc = nlat * nlon; %number of stations hadcrut4
nh = length(t_lat); %number of stations hadcrut3

hadcrut4 = reshape(H.d, [nloc, ntime]);
hadcrut4(hadcrut4<-1000) = NaN;

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



[x,y] = meshgrid(H.lon,H.lat);
loc = [x(:),y(:)];




c = pi/180;

%tic;
%find closest station and create tree-ring/station pair
for i = 1:nh %for each hadcrut3 loc
    for j = 1:nloc %for each hadcrut4 loc
        
        lat1 = c*t_lat(i); %had3 lat
        lon1 = c*t_lon(i); %had3 lon
        lat2 = c*loc(j,2); %had4 lat
        lon2 = c*loc(j,1); %had4 lon
        
        dh(i,j)=greatCircleDistance(lat1, lon1, lat2, lon2);
        
        [~,ind(i)]=min(d(i,:));%index of closest station
        
        %compile matching data
        htest(i).had3 = trdb_crn(i).yr;
        
       
        htest(i).had3 = trdb_crn(i).x;
        
    


        
        
    end
end