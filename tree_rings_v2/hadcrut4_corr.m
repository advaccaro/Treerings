%compute correlations for MXD/hadcrut4 pairs

%load('hcurve.mat')
nt = length(hcurve);
n = 0;
m = 0;
for i = 1:nt
    hcurve(i).smtestcheck = hcurve(i).locscreen + hcurve(i).mindatsm + hcurve(i).timescreen;
end
for i = 1:nt
    %if hcurve(i).finalscreen == 7
    if hcurve(i).smtestcheck == 4
        ind = ~isnan(hcurve(i).stdat) & ~isnan(hcurve(i).trdat);
        ind1 = ~isnan(hcurve(i).stdatp1) & ~isnan(hcurve(i).trdatp1);
        ind2 = ~isnan(hcurve(i).stdatp2) & ~isnan(hcurve(i).stdatp2);
        indsm = ~isnan(hcurve(i).stsm) & ~isnan(hcurve(i).trsm);
        indsm1 = ~isnan(hcurve(i).stsm1) & ~isnan(hcurve(i).trsm1);
        indsm2 = ~isnan(hcurve(i).stsm2) & ~isnan(hcurve(i).stsm2);
        
%         [hcurve(i).r, hcurve(i).signif, hcurve(i).p] = ...
%             corr_sig(hcurve(i).stdat(ind)*100, hcurve(i).trdat(ind));
%         [hcurve(i).r1, hcurve(i).signif1, hcurve(i).p1] = ...
%             corr_sig(hcurve(i).stdatp1(ind1)*100, hcurve(i).trdatp1(ind1));
%         [hcurve(i).r2, hcurve(i).signif2, hcurve(i).p2] = ...
%             corr_sig(hcurve(i).stdatp2(ind2)*100, hcurve(i).trdatp2(ind2));
         [hcurve(i).smr, hcurve(i).smsignif, hcurve(i).smp] = ...
             corr_sig(hcurve(i).stsm(indsm), hcurve(i).trsm(indsm));
         [hcurve(i).smr1, hcurve(i).smsignif1, hcurve(i).smp1] = ...
             corr_sig(hcurve(i).stsm1(indsm1), hcurve(i).trsm1(indsm1));
         [hcurve(i).smr2, hcurve(i).smsignif2, hcurve(i).smp2] = ...
             corr_sig(hcurve(i).stsm2(indsm2), hcurve(i).trsm2(indsm2));
    else
        display(i)
        n=n+1;
    end
end
sprintf('%d records screened out (quality control)', n)



% for i = 1:nt
%     if hcurve(i).finalscreen == 7
%         hcurve(i).dcorr = abs(hcurve(i).r1 - hcurve(i).r2);
%     else
%         hcurve(i).dcorr = NaN;
%     end
% end
% 
% for i = 1:nt
%     if hcurve(i).dcorr < .25
%         hcurve(i).corrscreen = 1;
%         m = m +1;
%     else
%         hcurve(i).corrscreen = 0;
%     end
% end

sprintf('%d records passed screen', m)


save('hcurve.mat', 'hcurve')