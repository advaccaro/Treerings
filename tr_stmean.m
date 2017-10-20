%create average from closet n stations

addpath(genpath('/Users/adam/Desktop/Treerings/'));
%load('tcurve.mat');
load('tr_db_crn.mat');
load('ghcn_db.mat');

nt = length(tcurve);

for i = 1:nt;
    if tcurve(i).locscreen == 1


        [smin idmin] = sort(d(i,:)); %sort d values
        minVals = smin(1:3); %lowest 3 d values
        minIds = idmin(1:3); %indices of lowest 3
        
        %get stations' data and years
        st1 =  ghcn(minIds(1)).da;
        st1yrs = round(ghcn(minIds(1)).years); 
        st2 = ghcn(minIds(2)).da;
        st2yrs = round(ghcn(minIds(2)).years);
        st3 = ghcn(minIds(3)).da;
        st3yrs = round(ghcn(minIds(3)).years);
        
        %find mins/maxs
        st1min = min(st1yrs); st2min = min(st2yrs); st3min = min(st3yrs);
        st1max = max(st1yrs); st2max = max(st2yrs); st3max = max(st3yrs);
        stmins = [st1min, st2min, st3min];
        stmaxs = [st1max, st2max, st3max];
        
        tcurve(i).stmin = min(stmins);
        tcurve(i).stmax = max(stmaxs);
        
        %st years from min/max
        tcurve(i).styrs = tcurve(i).stmin:tcurve(i).stmax;
        ny = length(tcurve(i).styrs);
        
        %preallocate w/ nan
        st1da = nan(1,ny); st2da = nan(1,ny); st3da = nan(1,ny);
        stdas = nan(3,ny);
        tcurve(i).stmean = nan(1,ny);
        
        
        %get indices
        st1ind1 = ismember(tcurve(i).styrs, st1yrs);
        st1ind2 = ismember(st1yrs, tcurve(i).styrs);
        st2ind1 = ismember(tcurve(i).styrs, st2yrs);
        st2ind2 = ismember(st2yrs, tcurve(i).styrs);
        st3ind1 = ismember(tcurve(i).styrs, st3yrs);
        st3ind2 = ismember(st3yrs, tcurve(i).styrs);
        
        %fill data
        st1da(st1ind1) = st1(st1ind2);
        st2da(st2ind1) = st2(st2ind2);
        st3da(st3ind1) = st3(st3ind2);
        
        A = vertcat(st1da, st2da, st3da);
        
        %compute stmean
        for j = 1:ny
            tcurve(i).stmean(j) = nmean(A(:,j));
        end
        
        else
        tcurve(i).trmin = nan; tcurve(i).trmax = nan;
        tcurve(i).stmin = nan; tcurve(i).stmax = nan;
        tcurve(i).stmean = nan;
        tcurve(i).years = nan(1,100);
    end


end
            