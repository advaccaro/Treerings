%fix tree ring ids to match standardized versions
addpath(genpath('/Users/adam/Desktop/Treerings/'))

load('tr_db_crn.mat')

newid = importdata('id.txt');
newid2 = cellstr(lower(newid(:)));
newid3 = deblank(newid2);




newlat = load('north.txt');

newlon = load('east.txt');

nt = length(trdb_crn);

nd = length(newid);

for i = 1:nt
    if trdb_crn(i).locscreen == 1
        id{i}=trdb_crn(i).id;
        
        
        ind = strcmp(id{i},newid3);
        if sum(ind) == 0
            display(i)
            display(trdb_crn(i).id)
            ncid = length(trdb_crn(i).id);
            if ncid > 10
                display('trimming 5 characters from %s', trdb_crn(i).id)
                badid = trdb_crn(i).id;
                trdb_crn(i).id = badid(1:end-5);
            end
            
        end
    end
end
%= regexprep('str', 'expr', 'repstr', options