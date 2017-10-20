addpath(genpath('/Users/adam/Desktop/Treerings/'))

%load GHCN monthly data
[A,TXT,RAW] = xlsread('ghcnm.tavg.v3.2.0.20121119.qca.inv.xls');
B = dlmread('ghcnm.tavg.v3.2.0.20121119.qca.dat.adj.txt' ,'\t');
B(find(B == -9999)) = NaN;

ns = size(A,1);%number of stations

%assemble GHCN raw data
for j = 1:ns
    ghcn(j).id=[A(j,1)]; ghcn(j).lat=[A(j,2)]; ghcn(j).lon=[A(j,3)];
    
    ind1 = (A(j,1)==floor(B(:,1)/10000));
    
    ghcn(j).years =((B(ind1,1)/10000)-floor(B(ind1,1)/10000))*10000;
    ghcn(j).yr = round(ghcn(j).years);
    
    %data in YxM
    ghcn(j).array1 = (B(ind1,2:13));
    
    
    
    
    %create fractional time axis
    
        ny=length(ghcn(j).years); %number of years
        nt=ny*12; %number of t-steps
        
        fmon = [0:11]/12;
        fmon2 = [repmat(fmon,ny,1)];
        
        fyr = ghcn(j).yr;
        fyr2 = [repmat(fyr,1,12)];
        
        fmat=fyr2+fmon2;
        
       ghcn(j).tfrac=sort(fmat(:));  
    
    
        
   %convert to datenum
      
      year = floor(ghcn(j).tfrac);
      month = floor((ghcn(j).tfrac - year)*12+1);
      
      %define datevec
      ts = [year, month, repmat(15, [nt 1])];
      ghcn(j).tn=datenum(ts);
      
     
      
  %convert data to time axis format
    
    
    array2 = ghcn(j).array1';
    ghcn(j).arrayt = reshape(array2,nt,1);
       
end



%save
save('ghcn_db.mat', 'ghcn')


