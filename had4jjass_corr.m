%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('had4jjassa.mat'); load('had4jjassd.mat');
na = length(had4jjassa); nd = length(had4jjassd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if had4jjassa(i).finalscreen == 7
        ind = ~isnan(had4jjassa(i).stdat) & ~isnan(had4jjassa(i).trdat);
        ind1 = ~isnan(had4jjassa(i).stdatp1) & ~isnan(had4jjassa(i).trdatp1);
        ind2 = ~isnan(had4jjassa(i).stdatp2) & ~isnan(had4jjassa(i).stdatp2);
        [had4jjassa(i).r, had4jjassa(i).signif, had4jjassa(i).p] = ...
            corr_sig(had4jjassa(i).stdat(ind), had4jjassa(i).trdat(ind));
        [had4jjassa(i).r1, had4jjassa(i).signif1, had4jjassa(i).p1] = ...
            corr_sig(had4jjassa(i).stdatp1(ind1), had4jjassa(i).trdatp1(ind1));
        [had4jjassa(i).r2, had4jjassa(i).signif2, had4jjassa(i).p2] = ...
            corr_sig(had4jjassa(i).stdatp2(ind2), had4jjassa(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    %had4jjassd(i).smtestcheck = had4jjassd(i).locscreen + had4jjassd(i).mindatsm + had4jjassd(i).timescreen;
if sum(~isnan(had4jjassd(i).stsm1)) > 15 & sum(~isnan(had4jjassd(i).stsm2)) > 10
	had4jjassd(i).smtestcheck = 4;
	else
had4jjassd(i).smtestcheck =0;
end
    
    if had4jjassd(i).smtestcheck == 4
        
        indsm = ~isnan(had4jjassd(i).stsm) & ~isnan(had4jjassd(i).trsm);
        indsm1 = ~isnan(had4jjassd(i).stsm1) & ~isnan(had4jjassd(i).trsm1);
        indsm2 = ~isnan(had4jjassd(i).stsm2) & ~isnan(had4jjassd(i).stsm2);
        
        
        [had4jjassd(i).smr, had4jjassd(i).smsignif, had4jjassd(i).smp] = ...
            corr_sig(had4jjassd(i).stsm(indsm), had4jjassd(i).trsm(indsm));
        [had4jjassd(i).smr1, had4jjassd(i).smsignif1, had4jjassd(i).smp1] = ...
            corr_sig(had4jjassd(i).stsm1(indsm1), had4jjassd(i).trsm1(indsm1));
        [had4jjassd(i).smr2, had4jjassd(i).smsignif2, had4jjassd(i).smp2] = ...
            corr_sig(had4jjassd(i).stsm2(indsm2), had4jjassd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('had4jjassa.mat', 'had4jjassa')
save('had4jjassd.mat', 'had4jjassd')

for i = 1:na
    if had4jjassa(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if had4jjassa(i).signif == 1
        naall = naall + 1;
    end
    if had4jjassa(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if had4jjassa(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if had4jjassa(i).signif1 == 1 & had4jjassa(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if had4jjassa(i).signif1 == 1 & had4jjassa(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if had4jjassa(i).signif1 == 0 & had4jjassa(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if had4jjassa(i).signif1 == 0 & had4jjassa(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if had4jjassd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if had4jjassd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if had4jjassd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if had4jjassd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if had4jjassd(i).smsignif1 == 1 & had4jjassd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if had4jjassd(i).smsignif1 == 1 & had4jjassd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if had4jjassd(i).smsignif1 == 0 & had4jjassd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if had4jjassd(i).smsignif1 == 0 & had4jjassd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save had4jjassa_results napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save had4jjassd_results ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2



clear all