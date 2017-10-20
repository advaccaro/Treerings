addpath(genpath('/Users/adam/Desktop/Treerings/'))

%load GHCN data
load('ghcn_db.mat');


%process data
    
       
 for j = 1:length(ghcn)
           
     %remove_season
     [ghcn(j).anom,ghcn(j).clim] = remove_season(ghcn(j).arrayt,ghcn(j).tn);
        
 


	%intra_annual_avg
        
   

	 [ghcn(j).da,ghcn(j).ta,ghcn(j).ts]=intra_annual_avg(ghcn(j).anom,ghcn(j).tfrac,6,8);
            %this is the final processed version
     %end
     %end
 end

%save
save('ghcn_db.mat', 'ghcn')


