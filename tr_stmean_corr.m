%tr_stmean_corr
%compute correlation between tr proxy and regional trend

%this script only works after tr_stmean and step2 have been used
load('tcurve.mat');
nt = length(tcurve);
n= 0;
for i = 1:nt
    %if i == 460
    %    tcurve(i).finalscreen = 0; %460 has DP but cannot compute corr
    %end
    if tcurve(i).finalscreen == 7
    ind=~isnan(tcurve(i).stdat) & ~isnan(tcurve(i).trdat);
    [tcurve(i).R, tcurve(i).pval]= corrcoef(tcurve(i).stmean(ind),tcurve(i).trtemp(ind));
    ind1=~isnan(tcurve(i).stdatp1)&~isnan(tcurve(i).trdatp1);
    [tcurve(i).Rp1,tcurve(i).pvalp1]= corrcoef(tcurve(i).stdatp1(ind1),tcurve(i).trdatp1(ind1));
%     ind2=~isnan(tcurve(i).stdatp2)&~isnan(tcurve(i).trdatp2);
%     [tcurve(i).Rp2,tcurve(i).sigp2,tcurve(i).pvalp2,tcurve(i).gp2]= corr_signif(tcurve(i).stdatp2(ind2),tcurve(i).trdatp2(ind2));
    else
        display(i)
        n=n+1;   
    end
end
sprintf('%d records screened out', n)

save('tcurve.mat', 'tcurve')
