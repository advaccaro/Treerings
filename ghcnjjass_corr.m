%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('ghcnjjassa.mat'); load('ghcnjjassd.mat');
na = length(ghcnjjassa); nd = length(ghcnjjassd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if ghcnjjassa(i).finalscreen == 7
        ind = ~isnan(ghcnjjassa(i).stdat) & ~isnan(ghcnjjassa(i).trdat);
        ind1 = ~isnan(ghcnjjassa(i).stdatp1) & ~isnan(ghcnjjassa(i).trdatp1);
        ind2 = ~isnan(ghcnjjassa(i).stdatp2) & ~isnan(ghcnjjassa(i).stdatp2);
        [ghcnjjassa(i).r, ghcnjjassa(i).signif, ghcnjjassa(i).p] = ...
            corr_sig(ghcnjjassa(i).stdat(ind), ghcnjjassa(i).trdat(ind));
        [ghcnjjassa(i).r1, ghcnjjassa(i).signif1, ghcnjjassa(i).p1] = ...
            corr_sig(ghcnjjassa(i).stdatp1(ind1), ghcnjjassa(i).trdatp1(ind1));
        [ghcnjjassa(i).r2, ghcnjjassa(i).signif2, ghcnjjassa(i).p2] = ...
            corr_sig(ghcnjjassa(i).stdatp2(ind2), ghcnjjassa(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    %ghcnjjassd(i).smtestcheck = ghcnjjassd(i).locscreen + ghcnjjassd(i).mindatsm + ghcnjjassd(i).timescreen;
if sum(~isnan(ghcnjjassd(i).stsm1)) > 15 & sum(~isnan(ghcnjjassd(i).stsm2)) > 10
	ghcnjjassd(i).smtestcheck = 4;
	else
ghcnjjassd(i).smtestcheck =0;
end
    
    if ghcnjjassd(i).smtestcheck == 4
        
        indsm = ~isnan(ghcnjjassd(i).stsm) & ~isnan(ghcnjjassd(i).trsm);
        indsm1 = ~isnan(ghcnjjassd(i).stsm1) & ~isnan(ghcnjjassd(i).trsm1);
        indsm2 = ~isnan(ghcnjjassd(i).stsm2) & ~isnan(ghcnjjassd(i).stsm2);
        
        
        [ghcnjjassd(i).smr, ghcnjjassd(i).smsignif, ghcnjjassd(i).smp] = ...
            corr_sig(ghcnjjassd(i).stsm(indsm), ghcnjjassd(i).trsm(indsm));
        [ghcnjjassd(i).smr1, ghcnjjassd(i).smsignif1, ghcnjjassd(i).smp1] = ...
            corr_sig(ghcnjjassd(i).stsm1(indsm1), ghcnjjassd(i).trsm1(indsm1));
        [ghcnjjassd(i).smr2, ghcnjjassd(i).smsignif2, ghcnjjassd(i).smp2] = ...
            corr_sig(ghcnjjassd(i).stsm2(indsm2), ghcnjjassd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('ghcnjjassa.mat', 'ghcnjjassa')
save('ghcnjjassd.mat', 'ghcnjjassd')

for i = 1:na
    if ghcnjjassa(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if ghcnjjassa(i).signif == 1
        naall = naall + 1;
    end
    if ghcnjjassa(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if ghcnjjassa(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if ghcnjjassa(i).signif1 == 1 & ghcnjjassa(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if ghcnjjassa(i).signif1 == 1 & ghcnjjassa(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if ghcnjjassa(i).signif1 == 0 & ghcnjjassa(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if ghcnjjassa(i).signif1 == 0 & ghcnjjassa(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if ghcnjjassd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if ghcnjjassd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if ghcnjjassd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if ghcnjjassd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if ghcnjjassd(i).smsignif1 == 1 & ghcnjjassd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if ghcnjjassd(i).smsignif1 == 1 & ghcnjjassd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if ghcnjjassd(i).smsignif1 == 0 & ghcnjjassd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if ghcnjjassd(i).smsignif1 == 0 & ghcnjjassd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save ghcnjjassa_results napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save ghcnjjassd_results ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2

clear all

