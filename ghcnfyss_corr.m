%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('ghcnfyssa.mat'); load('ghcnfyssd.mat');
na = length(ghcnfyssa); nd = length(ghcnfyssd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if ghcnfyssa(i).finalscreen == 7
        ind = ~isnan(ghcnfyssa(i).stdat) & ~isnan(ghcnfyssa(i).trdat);
        ind1 = ~isnan(ghcnfyssa(i).stdatp1) & ~isnan(ghcnfyssa(i).trdatp1);
        ind2 = ~isnan(ghcnfyssa(i).stdatp2) & ~isnan(ghcnfyssa(i).stdatp2);
        [ghcnfyssa(i).r, ghcnfyssa(i).signif, ghcnfyssa(i).p] = ...
            corr_sig(ghcnfyssa(i).stdat(ind), ghcnfyssa(i).trdat(ind));
        [ghcnfyssa(i).r1, ghcnfyssa(i).signif1, ghcnfyssa(i).p1] = ...
            corr_sig(ghcnfyssa(i).stdatp1(ind1), ghcnfyssa(i).trdatp1(ind1));
        [ghcnfyssa(i).r2, ghcnfyssa(i).signif2, ghcnfyssa(i).p2] = ...
            corr_sig(ghcnfyssa(i).stdatp2(ind2), ghcnfyssa(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    if sum(~isnan(ghcnfyssd(i).stsm1)) > 15 & sum(~isnan(ghcnfyssd(i).stsm2)) > 10
	ghcnfyssd(i).smtestcheck = 4;
	else
ghcnfyssd(i).smtestcheck =0;
end

    
    if ghcnfyssd(i).smtestcheck == 4
        
        indsm = ~isnan(ghcnfyssd(i).stsm) & ~isnan(ghcnfyssd(i).trsm);
        indsm1 = ~isnan(ghcnfyssd(i).stsm1) & ~isnan(ghcnfyssd(i).trsm1);
        indsm2 = ~isnan(ghcnfyssd(i).stsm2) & ~isnan(ghcnfyssd(i).stsm2);
        
        
        [ghcnfyssd(i).smr, ghcnfyssd(i).smsignif, ghcnfyssd(i).smp] = ...
            corr_sig(ghcnfyssd(i).stsm(indsm), ghcnfyssd(i).trsm(indsm));
        [ghcnfyssd(i).smr1, ghcnfyssd(i).smsignif1, ghcnfyssd(i).smp1] = ...
            corr_sig(ghcnfyssd(i).stsm1(indsm1), ghcnfyssd(i).trsm1(indsm1));
        [ghcnfyssd(i).smr2, ghcnfyssd(i).smsignif2, ghcnfyssd(i).smp2] = ...
            corr_sig(ghcnfyssd(i).stsm2(indsm2), ghcnfyssd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('ghcnfyssa.mat', 'ghcnfyssa')
save('ghcnfyssd.mat', 'ghcnfyssd')

for i = 1:na
    if ghcnfyssa(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if ghcnfyssa(i).signif == 1
        naall = naall + 1;
    end
    if ghcnfyssa(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if ghcnfyssa(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if ghcnfyssa(i).signif1 == 1 & ghcnfyssa(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if ghcnfyssa(i).signif1 == 1 & ghcnfyssa(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if ghcnfyssa(i).signif1 == 0 & ghcnfyssa(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if ghcnfyssa(i).signif1 == 0 & ghcnfyssa(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if ghcnfyssd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if ghcnfyssd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if ghcnfyssd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if ghcnfyssd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if ghcnfyssd(i).smsignif1 == 1 & ghcnfyssd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if ghcnfyssd(i).smsignif1 == 1 & ghcnfyssd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if ghcnfyssd(i).smsignif1 == 0 & ghcnfyssd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if ghcnfyssd(i).smsignif1 == 0 & ghcnfyssd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save ghcnfyssa_results napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save ghcnfyssd_results ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2

clear all



