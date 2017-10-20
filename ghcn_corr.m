%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('ghcna.mat'); load('ghcnd.mat');
na = length(ghcna); nd = length(ghcnd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if ghcna(i).finalscreen == 7
        ind = ~isnan(ghcna(i).stdat) & ~isnan(ghcna(i).trdat);
        ind1 = ~isnan(ghcna(i).stdatp1) & ~isnan(ghcna(i).trdatp1);
        ind2 = ~isnan(ghcna(i).stdatp2) & ~isnan(ghcna(i).stdatp2);
        [ghcna(i).r, ghcna(i).signif, ghcna(i).p] = ...
            corr_sig(ghcna(i).stdat(ind), ghcna(i).trdat(ind));
        [ghcna(i).r1, ghcna(i).signif1, ghcna(i).p1] = ...
            corr_sig(ghcna(i).stdatp1(ind1), ghcna(i).trdatp1(ind1));
        [ghcna(i).r2, ghcna(i).signif2, ghcna(i).p2] = ...
            corr_sig(ghcna(i).stdatp2(ind2), ghcna(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    if sum(~isnan(ghcnd(i).stsm1)) > 15 & sum(~isnan(ghcnd(i).stsm2)) > 10
	ghcnd(i).smtestcheck = 4;
	else
ghcnd(i).smtestcheck =0;
end

    
    if ghcnd(i).smtestcheck == 4
        
        indsm = ~isnan(ghcnd(i).stsm) & ~isnan(ghcnd(i).trsm);
        indsm1 = ~isnan(ghcnd(i).stsm1) & ~isnan(ghcnd(i).trsm1);
        indsm2 = ~isnan(ghcnd(i).stsm2) & ~isnan(ghcnd(i).stsm2);
        
        
        [ghcnd(i).smr, ghcnd(i).smsignif, ghcnd(i).smp] = ...
            corr_sig(ghcnd(i).stsm(indsm), ghcnd(i).trsm(indsm));
        [ghcnd(i).smr1, ghcnd(i).smsignif1, ghcnd(i).smp1] = ...
            corr_sig(ghcnd(i).stsm1(indsm1), ghcnd(i).trsm1(indsm1));
        [ghcnd(i).smr2, ghcnd(i).smsignif2, ghcnd(i).smp2] = ...
            corr_sig(ghcnd(i).stsm2(indsm2), ghcnd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('ghcna.mat', 'ghcna')
save('ghcnd.mat', 'ghcnd')

for i = 1:na
    if ghcna(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if ghcna(i).signif == 1
        naall = naall + 1;
    end
    if ghcna(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if ghcna(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if ghcna(i).signif1 == 1 & ghcna(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if ghcna(i).signif1 == 1 & ghcna(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if ghcna(i).signif1 == 0 & ghcna(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if ghcna(i).signif1 == 0 & ghcna(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if ghcnd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if ghcnd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if ghcnd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if ghcnd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if ghcnd(i).smsignif1 == 1 & ghcnd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if ghcnd(i).smsignif1 == 1 & ghcnd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if ghcnd(i).smsignif1 == 0 & ghcnd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if ghcnd(i).smsignif1 == 0 & ghcnd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save ghcna_results2 napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save ghcnd_results2 ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2



