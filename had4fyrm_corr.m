%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('had4fyrma.mat'); load('had4fyrmd.mat');
na = length(had4fyrma); nd = length(had4fyrmd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if had4fyrma(i).finalscreen == 7
        ind = ~isnan(had4fyrma(i).stdat) & ~isnan(had4fyrma(i).trdat);
        ind1 = ~isnan(had4fyrma(i).stdatp1) & ~isnan(had4fyrma(i).trdatp1);
        ind2 = ~isnan(had4fyrma(i).stdatp2) & ~isnan(had4fyrma(i).stdatp2);
        [had4fyrma(i).r, had4fyrma(i).signif, had4fyrma(i).p] = ...
            corr_sig(had4fyrma(i).stdat(ind), had4fyrma(i).trdat(ind));
        [had4fyrma(i).r1, had4fyrma(i).signif1, had4fyrma(i).p1] = ...
            corr_sig(had4fyrma(i).stdatp1(ind1), had4fyrma(i).trdatp1(ind1));
        [had4fyrma(i).r2, had4fyrma(i).signif2, had4fyrma(i).p2] = ...
            corr_sig(had4fyrma(i).stdatp2(ind2), had4fyrma(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    %had4fyrmd(i).smtestcheck = had4fyrmd(i).locscreen + had4fyrmd(i).mindatsm + had4fyrmd(i).timescreen;
if sum(~isnan(had4fyrmd(i).stsm1)) > 15 & sum(~isnan(had4fyrmd(i).stsm2)) > 10
	had4fyrmd(i).smtestcheck = 4;
	else
had4fyrmd(i).smtestcheck =0;
end
    
    if had4fyrmd(i).smtestcheck == 4
        
        indsm = ~isnan(had4fyrmd(i).stsm) & ~isnan(had4fyrmd(i).trsm);
        indsm1 = ~isnan(had4fyrmd(i).stsm1) & ~isnan(had4fyrmd(i).trsm1);
        indsm2 = ~isnan(had4fyrmd(i).stsm2) & ~isnan(had4fyrmd(i).stsm2);
        
        
        [had4fyrmd(i).smr, had4fyrmd(i).smsignif, had4fyrmd(i).smp] = ...
            corr_sig(had4fyrmd(i).stsm(indsm), had4fyrmd(i).trsm(indsm));
        [had4fyrmd(i).smr1, had4fyrmd(i).smsignif1, had4fyrmd(i).smp1] = ...
            corr_sig(had4fyrmd(i).stsm1(indsm1), had4fyrmd(i).trsm1(indsm1));
        [had4fyrmd(i).smr2, had4fyrmd(i).smsignif2, had4fyrmd(i).smp2] = ...
            corr_sig(had4fyrmd(i).stsm2(indsm2), had4fyrmd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('had4fyrma.mat', 'had4fyrma')
save('had4fyrmd.mat', 'had4fyrmd')

for i = 1:na
    if had4fyrma(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if had4fyrma(i).signif == 1
        naall = naall + 1;
    end
    if had4fyrma(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if had4fyrma(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if had4fyrma(i).signif1 == 1 & had4fyrma(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if had4fyrma(i).signif1 == 1 & had4fyrma(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if had4fyrma(i).signif1 == 0 & had4fyrma(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if had4fyrma(i).signif1 == 0 & had4fyrma(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if had4fyrmd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if had4fyrmd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if had4fyrmd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if had4fyrmd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if had4fyrmd(i).smsignif1 == 1 & had4fyrmd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if had4fyrmd(i).smsignif1 == 1 & had4fyrmd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if had4fyrmd(i).smsignif1 == 0 & had4fyrmd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if had4fyrmd(i).smsignif1 == 0 & had4fyrmd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save had4fyrma_results napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save had4fyrmd_results ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2


clear all
