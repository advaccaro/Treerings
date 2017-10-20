%uses Meko's crn2vec2 to assemble .crn's in database
%record 223 does not convert (this number may change if more
%records are added to the database)

addpath('/Users/adam/Desktop/Treerings/trunk/');

%folder w/ crn's
crnpath = '/Users/adam/Desktop/Treerings/data/rawdata/';



%select .crn files
list = dir(sprintf('%s*.crn',crnpath));


nc = length(list);





%convert and store in db
for j=1:nc
    
    
    nm1 = list(j).name; %.crn filename
    idleng = length(nm1); %length of filename
    
    if j == 223 %skip record 223; unable to open with crn2vec
        nm2 = sprintf('%s corrupt', nm1);
        screen(j) = 0;
   
    
    elseif idleng < 13
        nm2 = nm1(1:end-5); %trimmed site id for long id
        screen(j) =1;
    elseif idleng > 14
        nm2 = nm1(1:end-10); %trimmed site id for short id
        screen(j) =1;
    else
        display(j)
        display(nm1)
        nm2 = nm1(1:end-6);
        screen(j) =1;
    end
    
    trdb_crn(j).id=nm2; 




    if screen(j) == 1

    pfin = sprintf('%s',crnpath, nm1);
    
    
    [trdb_crn(j).x, trdb_crn(j).s, trdb_crn(j).yr] = crn2vec2(pfin);
    trdb_crn(j).locscreen = 1;
    else
        sprintf('unable to open %s', nm1)
        trdb_crn(j).locscreen = 0;
    end
    
end
    




 %preallocate space for lat/lon (need to add separately)
    
%  for m = 1:nc
%     trdb_crn(m).lat=0;trdb_crn(m).lon=0;
%  end
%  
 

save('tr_db_crn.mat','trdb_crn')