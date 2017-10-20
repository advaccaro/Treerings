addpath(genpath('/Users/adam/Desktop/Treerings/'))

%load tree ring data
load('tr_db_crn.mat');


%create fractional time axis
    nt = length(trdb_crn)  
    for j = 1:nt
        ny=length(trdb_crn(j).yr); %number of years
        %nt=ny*12; %number of t-steps
        
        %fmon = 5/12;
        %fmon2 = [repmat(fmon,ny,1)];
        
        %fyr = trdb_crn(j).yr;
        %fyr2 = [repmat(fyr,1,12)];
        
        %fmat=fyr2+fmon2;
        
      
        
       trdb_crn(j).tfrac=trdb_crn(j).yr(:);  
    
    
        
   %convert to datenum
      
      year = floor(trdb_crn(j).tfrac);
      month = floor((trdb_crn(j).tfrac - year)*12+1);
      
      %define datevec
      ts = [year, month, repmat(15, [ny 1])];
      trdb_crn(j).tn=datenum(ts);
      
     
      
  %convert data to time axis format
    
    trdb_crn(j).arrayt = trdb_crn(j).x;
       

    
    

%process data
    
  
           
     %remove_season
     [trdb_crn(j).anom,trdb_crn(j).clim] = remove_season(trdb_crn(j).arrayt,trdb_crn(j).tn);
        
 


	%intra_annual_avg
        
   

	 [trdb_crn(j).da,trdb_crn(j).ta,trdb_crn(j).ts]=intra_annual_avg(trdb_crn(j).anom,trdb_crn(j).tfrac,6,8);
            %this is the final processed version
     %end
     %end
 end

%save
save('tr_db_crn.mat', 'trdb_crn')