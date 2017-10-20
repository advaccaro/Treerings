%hadcrut4infgs_corr
%ghcngs_corr
%compute correlations for MXD/ghcngs pairs
load('ghcnfyrma.mat'); load('ghcnfyrmd.mat');
na = length(ghcnfyrma); nd = length(ghcnfyrmd);


n = 0;
m = 0;
napairs = 0; naall = 0; nap1 = 0; nap2 = 0; 
nap12 = 0; nap1x2 = 0; napx12 = 0; napx1x2 = 0;
ndpairs = 0; ndall = 0; ndp1 = 0; ndp2 = 0; 
ndp12 = 0; ndp1x2 = 0; ndpx12 = 0; ndpx1x2 = 0;


for i = 1:na
    if ghcnfyrma(i).finalscreen == 7
        ind = ~isnan(ghcnfyrma(i).stdat) & ~isnan(ghcnfyrma(i).trdat);
        ind1 = ~isnan(ghcnfyrma(i).stdatp1) & ~isnan(ghcnfyrma(i).trdatp1);
        ind2 = ~isnan(ghcnfyrma(i).stdatp2) & ~isnan(ghcnfyrma(i).stdatp2);
        [ghcnfyrma(i).r, ghcnfyrma(i).signif, ghcnfyrma(i).p] = ...
            corr_sig(ghcnfyrma(i).stdat(ind), ghcnfyrma(i).trdat(ind));
        [ghcnfyrma(i).r1, ghcnfyrma(i).signif1, ghcnfyrma(i).p1] = ...
            corr_sig(ghcnfyrma(i).stdatp1(ind1), ghcnfyrma(i).trdatp1(ind1));
        [ghcnfyrma(i).r2, ghcnfyrma(i).signif2, ghcnfyrma(i).p2] = ...
            corr_sig(ghcnfyrma(i).stdatp2(ind2), ghcnfyrma(i).trdatp2(ind2));
    else
        m = m + 1;
    end
end
for i = 1:nd
    if sum(~isnan(ghcnfyrmd(i).stsm1)) > 15 & sum(~isnan(ghcnfyrmd(i).stsm2)) > 10
	ghcnfyrmd(i).smtestcheck = 4;
	else
ghcnfyrmd(i).smtestcheck =0;
end

    
    if ghcnfyrmd(i).smtestcheck == 4
        
        indsm = ~isnan(ghcnfyrmd(i).stsm) & ~isnan(ghcnfyrmd(i).trsm);
        indsm1 = ~isnan(ghcnfyrmd(i).stsm1) & ~isnan(ghcnfyrmd(i).trsm1);
        indsm2 = ~isnan(ghcnfyrmd(i).stsm2) & ~isnan(ghcnfyrmd(i).stsm2);
        
        
        [ghcnfyrmd(i).smr, ghcnfyrmd(i).smsignif, ghcnfyrmd(i).smp] = ...
            corr_sig(ghcnfyrmd(i).stsm(indsm), ghcnfyrmd(i).trsm(indsm));
        [ghcnfyrmd(i).smr1, ghcnfyrmd(i).smsignif1, ghcnfyrmd(i).smp1] = ...
            corr_sig(ghcnfyrmd(i).stsm1(indsm1), ghcnfyrmd(i).trsm1(indsm1));
        [ghcnfyrmd(i).smr2, ghcnfyrmd(i).smsignif2, ghcnfyrmd(i).smp2] = ...
            corr_sig(ghcnfyrmd(i).stsm2(indsm2), ghcnfyrmd(i).trsm2(indsm2));
        
    else
        n=n+1;
    end
end

sprintf('%d annual records removed(quality control)', m)
sprintf('%d decadal records removed (quality control)', n)


save('ghcnfyrma.mat', 'ghcnfyrma')
save('ghcnfyrmd.mat', 'ghcnfyrmd')

for i = 1:na
    if ghcnfyrma(i).finalscreen == 7
        napairs = napairs + 1;
    end
    if ghcnfyrma(i).signif == 1
        naall = naall + 1;
    end
    if ghcnfyrma(i).signif1 == 1
        nap1 = nap1 + 1;
    end
    if ghcnfyrma(i).signif2 == 1
        nap2 = nap2 + 1;
    end
    if ghcnfyrma(i).signif1 == 1 & ghcnfyrma(i).signif2 == 1
        nap12 = nap12 + 1;
    end
    if ghcnfyrma(i).signif1 == 1 & ghcnfyrma(i).signif2 == 0
        nap1x2 = nap1x2 + 1;
    end
    if ghcnfyrma(i).signif1 == 0 & ghcnfyrma(i).signif2 == 1
        napx12 = napx12 + 1;
    end
    if ghcnfyrma(i).signif1 == 0 & ghcnfyrma(i).signif2 == 0
        napx1x2 = napx1x2 + 1;
    end
end

for i = 1:nd
    if ghcnfyrmd(i).smtestcheck == 4
        ndpairs = ndpairs + 1;
    end
    if ghcnfyrmd(i).smsignif == 1
        ndall = ndall + 1;
    end
    if ghcnfyrmd(i).smsignif1 == 1
        ndp1 = ndp1 + 1;
    end
    if ghcnfyrmd(i).smsignif2 == 1
        ndp2 = ndp2 + 1;
    end
    if ghcnfyrmd(i).smsignif1 == 1 & ghcnfyrmd(i).smsignif2 == 1
        ndp12 = ndp12 + 1;
    end
    if ghcnfyrmd(i).smsignif1 == 1 & ghcnfyrmd(i).smsignif2 == 0
        ndp1x2 = ndp1x2 + 1;
    end
    if ghcnfyrmd(i).smsignif1 == 0 & ghcnfyrmd(i).smsignif2 == 1
        ndpx12 = ndpx12 + 1;
    end
    if ghcnfyrmd(i).smsignif1 == 0 & ghcnfyrmd(i).smsignif2 == 0
        ndpx1x2 = ndpx1x2 + 1;
    end
end

save ghcnfyrma_results napairs naall nap1 nap2 nap12 nap1x2 napx12 napx1x2
save ghcnfyrmd_results ndpairs ndall ndp1 ndp2 ndp12 ndp1x2 ndpx12 ndpx1x2

clear all



