%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('had3a.mat'); load('had3d.mat');
na = length(had3a); nd = length(had3d);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if had3a(i).finalscreen == 7
        ind = ~isnan(had3a(i).stdat) & ~isnan(had3a(i).trdat);
        ind1 = ~isnan(had3a(i).stdatp1) & ~isnan(had3a(i).trdatp1);
        ind2 = ~isnan(had3a(i).stdatp2) & ~isnan(had3a(i).stdatp2);
        [had3a(i).r, had3a(i).signif, had3a(i).p] = ...
            corr_sig(had3a(i).stdat(ind), had3a(i).trdat(ind));
        [had3a(i).r1, had3a(i).signif1, had3a(i).p1] = ...
            corr_sig(had3a(i).stdatp1(ind1), had3a(i).trdatp1(ind1));
        [had3a(i).r2, had3a(i).signif2, had3a(i).p2] = ...
            corr_sig(had3a(i).stdatp2(ind2), had3a(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    %had3d(i).smtestcheck = had3d(i).locscreen + had3d(i).mindatsm + had3d(i).timescreen;
if sum(~isnan(had3d(i).stsm1)) > 15 & sum(~isnan(had3d(i).stsm2)) > 10
	had3d(i).smtestcheck = 4;
	else
had3d(i).smtestcheck =0;
end
    
    if had3d(i).smtestcheck == 4
        
        indsm = ~isnan(had3d(i).stsm) & ~isnan(had3d(i).trsm);
        indsm1 = ~isnan(had3d(i).stsm1) & ~isnan(had3d(i).trsm1);
        indsm2 = ~isnan(had3d(i).stsm2) & ~isnan(had3d(i).stsm2);
        
        
        [had3d(i).smr, had3d(i).smsignif, had3d(i).smp] = ...
            corr_sig(had3d(i).stsm(indsm), had3d(i).trsm(indsm));
        [had3d(i).smr1, had3d(i).smsignif1, had3d(i).smp1] = ...
            corr_sig(had3d(i).stsm1(indsm1), had3d(i).trsm1(indsm1));
        [had3d(i).smr2, had3d(i).smsignif2, had3d(i).smp2] = ...
            corr_sig(had3d(i).stsm2(indsm2), had3d(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('had3a.mat', 'had3a')
save('had3d.mat', 'had3d')

for i = 1:na
    if had3a(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if had3a(i).signif == 1
        naall = naall + 1;
    end
    if had3a(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if had3a(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if had3a(i).signif1 == 1 & had3a(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if had3a(i).signif1 == 1 & had3a(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if had3a(i).signif1 == 0 & had3a(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if had3a(i).signif1 == 0 & had3a(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if had3d(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if had3d(i).smsignif == 1
        ndall = ndall + 1;
    end
    if had3d(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if had3d(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if had3d(i).smsignif1 == 1 & had3d(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if had3d(i).smsignif1 == 1 & had3d(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if had3d(i).smsignif1 == 0 & had3d(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if had3d(i).smsignif1 == 0 & had3d(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save had3a_results2 napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save had3d_results2 ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2



