%this script smooths tree-ring and station data using a moving average.

load('tcurve.mat')

nt = length(tcurve);

for i = 1:nt
    if tcurve(i).finalscreen == 7
        %ind1 = isfinite(tcurve(i).trdat);
        %ind2 = isfinite(tcurve(i).stdat);
        tcurve(i).trmovavg = smooth(tcurve(i).trdat,10);
        tcurve(i).stmovavg = smooth(tcurve(i).stdat,10);
        
%         tcurve(i).trtrend = filter(1/5, [1 (1/5)-1], tcurve(i).trmovavg);
%         tcurve(i).sttrend = filter(1/5, [1 (1/5)-1], tcurve(i).stmovavg);
        tcurve(i).trttrend = filter(1/2, [1 (1/2)-1], tcurve(i).trdat);
        tcurve(i).sttrend = filter(1/2, [1 (1/2)-1], tcurve(i).stdat);

    end
end


save('tcurve.mat', 'tcurve')