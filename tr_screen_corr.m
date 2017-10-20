%compute correlation between tr proxy and closest station

%this script only works after tr_screen_step1 and step2 have been used
load('tcurve.mat');
nt = length(tcurve);
n= 0;
for i = 1:nt
    
    if tcurve(i).finalscreen == 7
    ind=~isnan(tcurve(i).stdat) & ~isnan(tcurve(i).trdat);
    %[tcurve(i).R,tcurve(i).sig,tcurve(i).pval,tcurve(i).g]= corr_signif(tcurve(i).sttemp(ind),tcurve(i).trtemp(ind));
    ind1=~isnan(tcurve(i).stdatp1)&~isnan(tcurve(i).trdatp1);
    [tcurve(i).Rp1,tcurve(i).sigp1,tcurve(i).pvalp1,tcurve(i).gp1]= corr_signif(tcurve(i).stdatp1(ind1),tcurve(i).trdatp1(ind1));
    ind2=~isnan(tcurve(i).stdatp2)&~isnan(tcurve(i).trdatp2);
    [tcurve(i).Rp2,tcurve(i).sigp2,tcurve(i).pvalp2,tcurve(i).gp2]= corr_signif(tcurve(i).stdatp2(ind2),tcurve(i).trdatp2(ind2));
    
    [tcurve(i).r, tcurve(i).signif, tcurve(i).p] = corr_sig(tcurve(i).sttemp(ind), tcurve(i).trtemp(ind));
    %add 2 more here
    
    else
        display(i)
        n=n+1;   
    end
end
sprintf('%d records screened out', n)

save('tcurve.mat', 'tcurve')
