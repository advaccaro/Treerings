%add lat/lon to existing treering db
addpath(genpath('/Users/adam/Desktop/Treerings/'))

% load('tr_db_crn.mat')

newid = importdata('id.txt');
newid2 = cellstr(lower(newid(:)));
newid3 = deblank(newid2);

nls = 0; %set locscreen counter = 0




newlat = load('north.txt');

newlon = load('east.txt');

nt = length(trdb_crn);

nd = length(newid);

for i = 1:nt
    id{i}=trdb_crn(i).id;
    
    
    ind = strcmp(id{i},newid3);
    
    if sum(ind) == 0
        display(i)
        display(trdb_crn(i).id)
        trdb_crn(i).locscreen = 0;
        nls = nls + 1;
    else
        trdb_crn(i).locscreen = 1;
        trdb_crn(i).lat=newlat(ind==1);
        trdb_crn(i).lon=newlon(ind==1);
        
    end
    

end

%manually add missing lat/lon info
trdb_crn(23).lat = 57;
trdb_crn(23).lon = -3.58;
trdb_crn(23).locscreen = 1;


trdb_crn(60).lat = 61.68;
trdb_crn(60).lon = -120.72;
trdb_crn(60).locscreen = 1;

trdb_crn(110).lat = 51.42;
trdb_crn(110).lon = -117.33;
trdb_crn(110).locscreen = 1;

trdb_crn(114).lat = 53.95;
trdb_crn(114).lon = -105.15;
trdb_crn(114).locscreen = 1;

trdb_crn(141).lat = 50.75;
trdb_crn(141).lon = 15.55;
trdb_crn(141).locscreen = 1;

trdb_crn(142).lat = 50.75;
trdb_crn(142).lon = 15.55;
trdb_crn(142).locscreen =1;


nls=nls-5;

display(nls)


save('tr_db_crn.mat', 'trdb_crn')