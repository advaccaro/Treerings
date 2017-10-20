%check for missing t-tick in hadcrut4 dataset

for i = 1:floor(1955/12)
    startval = (i-1).*12+1;
    stopval = i.*12;
    
    testvec = floor(H.time(startval:stopval)/365);
    testval = floor(H.time(startval)/365);
    valvec = repmat (testval,12,1);
    
    if ~isequal(valvec, testvec)
        display(i)
    end
end