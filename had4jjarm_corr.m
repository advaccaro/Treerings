%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('had4jjarma.mat'); load('had4jjarmd.mat');
na = length(had4jjarma); nd = length(had4jjarmd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if had4jjarma(i).finalscreen == 7
        ind = ~isnan(had4jjarma(i).stdat) & ~isnan(had4jjarma(i).trdat);
        ind1 = ~isnan(had4jjarma(i).stdatp1) & ~isnan(had4jjarma(i).trdatp1);
        ind2 = ~isnan(had4jjarma(i).stdatp2) & ~isnan(had4jjarma(i).stdatp2);
        [had4jjarma(i).r, had4jjarma(i).signif, had4jjarma(i).p] = ...
            corr_sig(had4jjarma(i).stdat(ind), had4jjarma(i).trdat(ind));
        [had4jjarma(i).r1, had4jjarma(i).signif1, had4jjarma(i).p1] = ...
            corr_sig(had4jjarma(i).stdatp1(ind1), had4jjarma(i).trdatp1(ind1));
        [had4jjarma(i).r2, had4jjarma(i).signif2, had4jjarma(i).p2] = ...
            corr_sig(had4jjarma(i).stdatp2(ind2), had4jjarma(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    %had4jjarmd(i).smtestcheck = had4jjarmd(i).locscreen + had4jjarmd(i).mindatsm + had4jjarmd(i).timescreen;
if sum(~isnan(had4jjarmd(i).stsm1)) > 15 & sum(~isnan(had4jjarmd(i).stsm2)) > 10
	had4jjarmd(i).smtestcheck = 4;
	else
had4jjarmd(i).smtestcheck =0;
end
    
    if had4jjarmd(i).smtestcheck == 4
        
        indsm = ~isnan(had4jjarmd(i).stsm) & ~isnan(had4jjarmd(i).trsm);
        indsm1 = ~isnan(had4jjarmd(i).stsm1) & ~isnan(had4jjarmd(i).trsm1);
        indsm2 = ~isnan(had4jjarmd(i).stsm2) & ~isnan(had4jjarmd(i).stsm2);
        
        
        [had4jjarmd(i).smr, had4jjarmd(i).smsignif, had4jjarmd(i).smp] = ...
            corr_sig(had4jjarmd(i).stsm(indsm), had4jjarmd(i).trsm(indsm));
        [had4jjarmd(i).smr1, had4jjarmd(i).smsignif1, had4jjarmd(i).smp1] = ...
            corr_sig(had4jjarmd(i).stsm1(indsm1), had4jjarmd(i).trsm1(indsm1));
        [had4jjarmd(i).smr2, had4jjarmd(i).smsignif2, had4jjarmd(i).smp2] = ...
            corr_sig(had4jjarmd(i).stsm2(indsm2), had4jjarmd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('had4jjarma.mat', 'had4jjarma')
save('had4jjarmd.mat', 'had4jjarmd')

for i = 1:na
    if had4jjarma(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if had4jjarma(i).signif == 1
        naall = naall + 1;
    end
    if had4jjarma(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if had4jjarma(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if had4jjarma(i).signif1 == 1 & had4jjarma(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if had4jjarma(i).signif1 == 1 & had4jjarma(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if had4jjarma(i).signif1 == 0 & had4jjarma(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if had4jjarma(i).signif1 == 0 & had4jjarma(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if had4jjarmd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if had4jjarmd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if had4jjarmd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if had4jjarmd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if had4jjarmd(i).smsignif1 == 1 & had4jjarmd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if had4jjarmd(i).smsignif1 == 1 & had4jjarmd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if had4jjarmd(i).smsignif1 == 0 & had4jjarmd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if had4jjarmd(i).smsignif1 == 0 & had4jjarmd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save had4jjarma_results napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save had4jjarmd_results ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2



clear all