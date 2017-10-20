%prepares MXD data for plotting by splitting the data in three groups: raw,
%quality controlled, and screened

addpath(genpath('/Users/adam/Desktop/Treerings/'));

load('tr_db_crn.mat'); load('tcurve.mat'); nt = length(tcurve);
%separate the data into three groups: 3=screened, 2=quality controlled, 1=
%failed quality control
for i = 1:nt
    tcurve(i).group = 0; %preset all groups = 0
end

for i = 1:nt
    if tcurve(i).corrscreen == 1
        tcurve(i).group = 3;
    else
        if tcurve(i).finalscreen == 7
        tcurve(i).group = 2;
        else
        if tcurve(i).locscreen == 1
        tcurve(i).group = 1;
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
       indgroup(j,i) = (tcurve(i).group>=j);
    end
end

%preallocate/predefine
tcurve(223).trmin = NaN; tcurve(223).trmax = NaN;

for j = 1:3
    group(j).yearmin = zeros(nt,1);
    %group(j).yearmin(223) = NaN; 
    group(j).yearmax = zeros(nt,1);
    %group(j).yearmax(223) = NaN; 
    
    group(j).minyears = cell2mat({tcurve.trmin}); %vector containing min values
    group(j).maxyears = cell2mat({tcurve.trmax});%vec containing max values

end





for j = 1:3 %change this to 3 when group 3 is done!            
            
%find year_i/year_f
group(j).yearmin(indgroup(j,:)) = group(j).minyears(indgroup(j,:)); %min year vector
indmin = group(j).yearmin(:)>1; %ind for min year values (>1)
group(j).year_i = min(group(j).yearmin(indmin));%min year value

group(j).yearmax(indgroup(j,:)) = group(j).maxyears(indgroup(j,:)); %max year vector
indmax = group(j).yearmax(:)<2014; %ind for max year values (<2014)
group(j).year_f = max(group(j).yearmax(indmax));%max year value


%create year_vec
group(j).year_vec = group(j).year_i:group(j).year_f; %vector containing all relevant years
group(j).ny = length(group(j).year_vec); % # of years
group(j).tm = group(j).year_vec;

tcurve(223).trtemp = repmat(-9999, group(j).ny);


%preallocate for mxd&nproxy

group(j).mxd = nan(nt, group(1).ny);

group(j).nproxy = nan(group(1).ny,1);
end

for j = 1:3


for i = 1:nt
    if tcurve(i).locscreen == 0
        group(j).mxd(i,:) = nan(1,group(1).ny);
    elseif tcurve(i).group >= j
%     group(j).indyr = ismember(group(1).year_vec,tcurve(i).tryears);
%     group(j).indyr2 = ismember(tcurve(i).tryears,group(1).year_vec);
    indyr = ismember(group(1).year_vec, tcurve(i).tryears);
    indyr2 = ismember(tcurve(i).tryears, group(1).year_vec);
%     group(j).mxd(i,group(1).indyr) = tcurve(i).trtemp(group(j).indyr2);
    group(j).mxd(i,indyr) = tcurve(i).trtemp(indyr2);
    end
end



%create nproxy array (for bargraph)
for i = 1:group(1).ny
    group(j).nproxy(i) = sum(~isnan(group(j).mxd(:,i)));
end



% ptype = b.ptype; np = size(ptype,2);
group(j).proxy = group(j).mxd; 

end


%screen lat/lons
lon1 = cell2mat({tcurve.lon}); %lon of unscreened trees
lat1 = cell2mat({tcurve.lat}); %lat of unscreened trees


for j = 1:3
    n = 1; %set counter = 1
for i = 1:nt
    if tcurve(i).group >= j
        group(j).lon(n) = lon1(i);
        group(j).lat(n) = lat1(i);
        n = n+1;
    end
end
end


%counters of # members/group
ng1 = nt;
ng2 = sum(indgroup(2,:));
ng3 = sum(indgroup(3,:));

save('tcurve.mat','tcurve')
save('group.mat', 'group')