%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('ghcnjjarma.mat'); load('ghcnjjarmd.mat');
na = length(ghcnjjarma); nd = length(ghcnjjarmd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if ghcnjjarma(i).finalscreen == 7
        ind = ~isnan(ghcnjjarma(i).stdat) & ~isnan(ghcnjjarma(i).trdat);
        ind1 = ~isnan(ghcnjjarma(i).stdatp1) & ~isnan(ghcnjjarma(i).trdatp1);
        ind2 = ~isnan(ghcnjjarma(i).stdatp2) & ~isnan(ghcnjjarma(i).stdatp2);
        [ghcnjjarma(i).r, ghcnjjarma(i).signif, ghcnjjarma(i).p] = ...
            corr_sig(ghcnjjarma(i).stdat(ind), ghcnjjarma(i).trdat(ind));
        [ghcnjjarma(i).r1, ghcnjjarma(i).signif1, ghcnjjarma(i).p1] = ...
            corr_sig(ghcnjjarma(i).stdatp1(ind1), ghcnjjarma(i).trdatp1(ind1));
        [ghcnjjarma(i).r2, ghcnjjarma(i).signif2, ghcnjjarma(i).p2] = ...
            corr_sig(ghcnjjarma(i).stdatp2(ind2), ghcnjjarma(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    %ghcnjjarmd(i).smtestcheck = ghcnjjarmd(i).locscreen + ghcnjjarmd(i).mindatsm + ghcnjjarmd(i).timescreen;
if sum(~isnan(ghcnjjarmd(i).stsm1)) > 15 & sum(~isnan(ghcnjjarmd(i).stsm2)) > 10
	ghcnjjarmd(i).smtestcheck = 4;
	else
ghcnjjarmd(i).smtestcheck =0;
end
    
    if ghcnjjarmd(i).smtestcheck == 4
        
        indsm = ~isnan(ghcnjjarmd(i).stsm) & ~isnan(ghcnjjarmd(i).trsm);
        indsm1 = ~isnan(ghcnjjarmd(i).stsm1) & ~isnan(ghcnjjarmd(i).trsm1);
        indsm2 = ~isnan(ghcnjjarmd(i).stsm2) & ~isnan(ghcnjjarmd(i).stsm2);
        
        
        [ghcnjjarmd(i).smr, ghcnjjarmd(i).smsignif, ghcnjjarmd(i).smp] = ...
            corr_sig(ghcnjjarmd(i).stsm(indsm), ghcnjjarmd(i).trsm(indsm));
        [ghcnjjarmd(i).smr1, ghcnjjarmd(i).smsignif1, ghcnjjarmd(i).smp1] = ...
            corr_sig(ghcnjjarmd(i).stsm1(indsm1), ghcnjjarmd(i).trsm1(indsm1));
        [ghcnjjarmd(i).smr2, ghcnjjarmd(i).smsignif2, ghcnjjarmd(i).smp2] = ...
            corr_sig(ghcnjjarmd(i).stsm2(indsm2), ghcnjjarmd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('ghcnjjarma.mat', 'ghcnjjarma')
save('ghcnjjarmd.mat', 'ghcnjjarmd')

for i = 1:na
    if ghcnjjarma(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if ghcnjjarma(i).signif == 1
        naall = naall + 1;
    end
    if ghcnjjarma(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if ghcnjjarma(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if ghcnjjarma(i).signif1 == 1 & ghcnjjarma(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if ghcnjjarma(i).signif1 == 1 & ghcnjjarma(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if ghcnjjarma(i).signif1 == 0 & ghcnjjarma(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if ghcnjjarma(i).signif1 == 0 & ghcnjjarma(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if ghcnjjarmd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if ghcnjjarmd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if ghcnjjarmd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if ghcnjjarmd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if ghcnjjarmd(i).smsignif1 == 1 & ghcnjjarmd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if ghcnjjarmd(i).smsignif1 == 1 & ghcnjjarmd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if ghcnjjarmd(i).smsignif1 == 0 & ghcnjjarmd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if ghcnjjarmd(i).smsignif1 == 0 & ghcnjjarmd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save ghcnjjarma_results napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save ghcnjjarmd_results ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2

clear all

