%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('ghcngsa.mat'); load('ghcngsd.mat');
na = length(ghcngsa); nd = length(ghcngsd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if ghcngsa(i).finalscreen == 7
        ind = ~isnan(ghcngsa(i).stdat) & ~isnan(ghcngsa(i).trdat);
        ind1 = ~isnan(ghcngsa(i).stdatp1) & ~isnan(ghcngsa(i).trdatp1);
        ind2 = ~isnan(ghcngsa(i).stdatp2) & ~isnan(ghcngsa(i).stdatp2);
        [ghcngsa(i).r, ghcngsa(i).signif, ghcngsa(i).p] = ...
            corr_sig(ghcngsa(i).stdat(ind), ghcngsa(i).trdat(ind));
        [ghcngsa(i).r1, ghcngsa(i).signif1, ghcngsa(i).p1] = ...
            corr_sig(ghcngsa(i).stdatp1(ind1), ghcngsa(i).trdatp1(ind1));
        [ghcngsa(i).r2, ghcngsa(i).signif2, ghcngsa(i).p2] = ...
            corr_sig(ghcngsa(i).stdatp2(ind2), ghcngsa(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    %ghcngsd(i).smtestcheck = ghcngsd(i).locscreen + ghcngsd(i).mindatsm + ghcngsd(i).timescreen;
if sum(~isnan(ghcngsd(i).stsm1)) > 15 & sum(~isnan(ghcngsd(i).stsm2)) > 10
	ghcngsd(i).smtestcheck = 4;
	else
ghcngsd(i).smtestcheck =0;
end
    
    if ghcngsd(i).smtestcheck == 4
        
        indsm = ~isnan(ghcngsd(i).stsm) & ~isnan(ghcngsd(i).trsm);
        indsm1 = ~isnan(ghcngsd(i).stsm1) & ~isnan(ghcngsd(i).trsm1);
        indsm2 = ~isnan(ghcngsd(i).stsm2) & ~isnan(ghcngsd(i).stsm2);
        
        
        [ghcngsd(i).smr, ghcngsd(i).smsignif, ghcngsd(i).smp] = ...
            corr_sig(ghcngsd(i).stsm(indsm), ghcngsd(i).trsm(indsm));
        [ghcngsd(i).smr1, ghcngsd(i).smsignif1, ghcngsd(i).smp1] = ...
            corr_sig(ghcngsd(i).stsm1(indsm1), ghcngsd(i).trsm1(indsm1));
        [ghcngsd(i).smr2, ghcngsd(i).smsignif2, ghcngsd(i).smp2] = ...
            corr_sig(ghcngsd(i).stsm2(indsm2), ghcngsd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('ghcngsa.mat', 'ghcngsa')
save('ghcngsd.mat', 'ghcngsd')

for i = 1:na
    if ghcngsa(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if ghcngsa(i).signif == 1
        naall = naall + 1;
    end
    if ghcngsa(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if ghcngsa(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if ghcngsa(i).signif1 == 1 & ghcngsa(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if ghcngsa(i).signif1 == 1 & ghcngsa(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if ghcngsa(i).signif1 == 0 & ghcngsa(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if ghcngsa(i).signif1 == 0 & ghcngsa(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if ghcngsd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if ghcngsd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if ghcngsd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if ghcngsd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if ghcngsd(i).smsignif1 == 1 & ghcngsd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if ghcngsd(i).smsignif1 == 1 & ghcngsd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if ghcngsd(i).smsignif1 == 0 & ghcngsd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if ghcngsd(i).smsignif1 == 0 & ghcngsd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save ghcngsa_results2 napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save ghcngsd_results2 ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2



