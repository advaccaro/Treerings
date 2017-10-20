%tr_stmean_smooth
load('tcurve.mat')

nt = length(tcurve);

for i = 1:nt
    if tcurve(i).finalscreen == 7
        tcurve(i).trmovavg = smooth(tcurve(i).trdat,7);
        tcurve(i).stmovavg = smooth(tcurve(i).stdat,7);
        
%         tcurve(i).trtrend = filter(1/5, [1 (1/5)-1], tcurve(i).trmovavg);
%         tcurve(i).sttrend = filter(1/5, [1 (1/5)-1], tcurve(i).stmovavg);
        tcurve(i).trtrend = filter(1/2, [1 (1/2)-1], tcurve(i).trdat);
        tcurve(i).sttrend = filter(1/2, [1 (1/2)-1], tcurve(i).stdat);

    end
end


save('tcurve.mat', 'tcurve')