%prepare hadcrut4 data for MXD availability map

%prepares MXD data for plotting by splitting the data in three groups: raw,
%quality controlled, and screened

addpath(genpath('/Users/adam/Desktop/Treerings/'));

load('tr_db_crn.mat'); load('hcurve.mat'); nt = length(hcurve);
%separate the data into three groups: 3=screened, 2=quality controlled, 1=
%failed quality control
for i = 1:nt
    hcurve(i).group = 0; %preset all groups = 0
end

for i = 1:nt
    if hcurve(i).smsignif1 == 1 & hcurve(i).smsignif2 == 1
        hcurve(i).group = 3;
    else
        if hcurve(i).smsignif1 == 1
        hcurve(i).group = 2;
        else
        if hcurve(i).locscreen == 1
        hcurve(i).group = 1;
        end
        end
    end
end


%CALCULATE
%set min/max years (first, create array, get index, fill array, ->
%min/max)


%create index to match data to group
for i =1:nt
    for j = 1:3
        %indgroup(j,i) = (isequal(tcurve(i).group, j));
       indhgroup(j,i) = (hcurve(i).group>=j);
    end
end

%preallocate/predefine
hcurve(223).trmin = NaN; hcurve(223).trmax = NaN;

for j = 1:3
    hgroup(j).yearmin = zeros(nt,1);
    %group(j).yearmin(223) = NaN; 
    hgroup(j).yearmax = zeros(nt,1);
    %group(j).yearmax(223) = NaN; 
    
    hgroup(j).minyears = cell2mat({hcurve.trmin}); %vector containing min values
    hgroup(j).maxyears = cell2mat({hcurve.trmax});%vec containing max values

end





for j = 1:3 %change this to 3 when group 3 is done!            
            
%find year_i/year_f
hgroup(j).yearmin(indhgroup(j,:)) = hgroup(j).minyears(indhgroup(j,:)); %min year vector
indmin = hgroup(j).yearmin(:)>1; %ind for min year values (>1)
hgroup(j).year_i = min(hgroup(j).yearmin(indmin));%min year value

hgroup(j).yearmax(indhgroup(j,:)) = hgroup(j).maxyears(indhgroup(j,:)); %max year vector
indmax = hgroup(j).yearmax(:)<2014; %ind for max year values (<2014)
hgroup(j).year_f = max(hgroup(j).yearmax(indmax));%max year value


%create year_vec
hgroup(j).year_vec = hgroup(j).year_i:hgroup(j).year_f; %vector containing all relevant years
hgroup(j).ny = length(hgroup(j).year_vec); % # of years
hgroup(j).tm = hgroup(j).year_vec;

hcurve(223).trtemp = repmat(-9999, hgroup(j).ny);


%preallocate for mxd&nproxy

hgroup(j).mxd = nan(nt, hgroup(1).ny);

hgroup(j).nproxy = nan(hgroup(1).ny,1);
end

for j = 1:3


for i = 1:nt
    if hcurve(i).locscreen == 0
        hgroup(j).mxd(i,:) = nan(1,hgroup(1).ny);
    elseif hcurve(i).group >= j
%     group(j).indyr = ismember(group(1).year_vec,tcurve(i).tryears);
%     group(j).indyr2 = ismember(tcurve(i).tryears,group(1).year_vec);
    indyr = ismember(hgroup(1).year_vec, hcurve(i).tryears);
    indyr2 = ismember(hcurve(i).tryears, hgroup(1).year_vec);
%     group(j).mxd(i,group(1).indyr) = tcurve(i).trtemp(group(j).indyr2);
    hgroup(j).mxd(i,indyr) = hcurve(i).trtemp(indyr2);
    end
end



%create nproxy array (for bargraph)
for i = 1:hgroup(1).ny
    hgroup(j).nproxy(i) = sum(~isnan(hgroup(j).mxd(:,i)));
end



% ptype = b.ptype; np = size(ptype,2);
hgroup(j).proxy = hgroup(j).mxd; 

end


%screen lat/lons
lon1 = cell2mat({hcurve.lon}); %lon of unscreened trees
lat1 = cell2mat({hcurve.lat}); %lat of unscreened trees


for j = 1:3
    n = 1; %set counter = 1
for i = 1:nt
    if hcurve(i).group >= j
        hgroup(j).lon(n) = lon1(i);
        hgroup(j).lat(n) = lat1(i);
        n = n+1;
    end
end
end


%counters of # members/group
ng1 = nt;
ng2 = sum(indhgroup(2,:));
ng3 = sum(indhgroup(3,:));

save('hcurve.mat','hcurve')
save('hgroup.mat', 'hgroup')