addpath(genpath('/Users/adam/Desktop/Treerings/'))
load('tr_db_crn.mat');
load('tcurve.mat');

%STEP 2: trim years and data


nt = length(trdb_crn); %number of tree sites

%set counters for screened records = 0
nls = 0; %loc check
nts1 = 0;  %max year check
nts2 = 0; %min year check
nmd = 0; %minimum data check
ns = 0; %finalscreen
%load('tcurve.mat');





%Set tcurve year range by finding min/max values
for h = 1:nt
    
    %store tr site metadata in tcurve 
    tcurve(h).id=trdb_crn(h).id; 
    tcurve(h).lat=trdb_crn(h).lat;
    tcurve(h).lon=trdb_crn(h).lon;
        
    if tcurve(h).locscreen == 1

    
    


%find time bounds for tr & st
    tcurve(h).stmin = min(tcurve(h).styears); 
    tcurve(h).stmax = max(tcurve(h).styears);
    tcurve(h).trmin = min(tcurve(h).tryears); 
    tcurve(h).trmax = max(tcurve(h).tryears);
    
    
%use tr & st min/max values to calculate tmin/max
    tcurve(h).tmin = max(tcurve(h).stmin, tcurve(h).trmin);
    tcurve(h).tmax = min(tcurve(h).stmax, tcurve(h).trmax);
    
%set tcurve year range from tmin/max values
    tcurve(h).years = tcurve(h).tmin:tcurve(h).tmax;

    
%     (set this aside for now, probably best later in script)
%     tcurve(h).stdat = tcurve(h).sttemp(ismember(
%     tcurve(h).trdat = tcurve(h).trtemp(ismember(
%     
    else
        tcurve(h).trmin = nan; tcurve(h).trmax = nan;
        tcurve(h).stmin = nan; tcurve(h).stmax = nan;
        tcurve(h).years = nan(1,100);
    end
    
end



%timescreen using tmin/max values
%(can rewrite this section to include a screen for overlap w/ overlap as
%variable)

for j=1:nt
%     if j == 223 %skip missing record 223
%         tcurve(j).timescreen = 0;
    if tcurve(j).locscreen == 0
        nls=nls+1;
        tcurve(j).timescreen = 0;
    else
       
        tcurve(j).overlap = tcurve(j).tmax - tcurve(j).tmin;
        if tcurve(j).tmax < 1960
            screen1 = 0;
            nts1 = nts1 + 1;
        else
            screen1 = 1;
        end
        
        if tcurve(j).tmin > 1935
            screen2 = 0;
            nts2 = nts2 + 1;
            
        else
            screen2 = 1;
        end
        
        end    
        tcurve(j).timescreen = screen1 + screen2;
    
        
        %split suitable records into pre/post 1960
        
        if tcurve(j).timescreen == 2
            
            %split years
            tcurve(j).yearsp1 = [tcurve(j).tmin:1950];
            tcurve(j).yearsp2 = [1951:tcurve(j).tmax];
            
            
         
            
            %preallocate for st/tr years
            np1 = length(tcurve(j).yearsp1);
            np2 = length(tcurve(j).yearsp2);
            
            tcurve(j).stdatp1 = nan(np1,1);
            tcurve(j).stdatp2 = nan(np2,1);
            tcurve(j).trdatp1 = nan(np1,1);
            tcurve(j).trdatp2 = nan(np2,1);
            
            
            
            
            %find indices corresponding to those time intervals
            tr1ind = ismember(tcurve(j).yearsp1, tcurve(j).tryears);
            tr2ind = ismember(tcurve(j).yearsp2, tcurve(j).tryears);
            st1ind = ismember(tcurve(j).yearsp1, tcurve(j).styears);
            st2ind = ismember(tcurve(j).yearsp2, tcurve(j).styears);
            
            tr1ind2 = ismember(tcurve(j).tryears, tcurve(j).yearsp1);
            tr2ind2 = ismember(tcurve(j).tryears, tcurve(j).yearsp2);
            st1ind2 = ismember(tcurve(j).styears, tcurve(j).yearsp1);
            st2ind2 = ismember(tcurve(j).styears, tcurve(j).yearsp2);
            
            

            
            %get matching data
            tcurve(j).trdatp1(tr1ind) = tcurve(j).trtemp(tr1ind2);
            tcurve(j).trdatp2(tr2ind) = tcurve(j).trtemp(tr2ind2);
            tcurve(j).stdatp1(st1ind) = tcurve(j).sttemp(st1ind2);
            tcurve(j).stdatp2(st2ind) = tcurve(j).sttemp(st2ind2);
            
            tcurve(j).trdat = vertcat(tcurve(j).trdatp1, tcurve(j).trdatp2);
            tcurve(j).stdat = vertcat(tcurve(j).stdatp1, tcurve(j).stdatp2);
        
    
 
                        
            %min data screen p1
            p1screen = 0;
            p2screen = 0;
    
            ntrp1 = sum(~isnan(tcurve(j).trdatp1));
            ntrp2 = sum(~isnan(tcurve(j).trdatp2));
            nstp1 = sum(~isnan(tcurve(j).stdatp1));
            nstp2 = sum(~isnan(tcurve(j).stdatp2));
   
            if ntrp1>15
                p1screen=p1screen+1;
            end
    
            if nstp1>15
                p1screen=p1screen+1;
            end
    
            if ntrp2>10
                p2screen=p2screen+1;
            end
    
            if nstp2>10
                p2screen=p2screen+1;
            end
    
            tcurve(j).mindatscreen = p1screen + p2screen;
            if tcurve(j).mindatscreen < 4
                nmd = nmd +1;
                display(j)
                display('Min dat screened')
            end
        else
            tcurve(j).mindatscreen = 0;
            
    end
    
    tcurve(j).finalscreen = ...
        tcurve(j).timescreen + tcurve(j).locscreen + tcurve(j).mindatscreen;
    
    
    if tcurve(j).finalscreen < 7
        display('record screened out')
        display(j)
        display(tcurve(j).id)
        ns = ns +1;
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

save('tcurve.mat','tcurve');
