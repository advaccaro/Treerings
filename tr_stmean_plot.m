%tr_stmean_plot
function [trplot3] = tr_stmean_plot(data,n)
%this script only works after tr_screen_step1 and step2 have been used
load('tcurve.mat');
tcurve = data;

if tcurve(n).finalscreen < 7
    display('record failed screening')
    display(n)
    display(tcurve(n).id)
else

trmean = nmean(tcurve(n).trdat);
trmean2 = nmean(tcurve(n).trdat);


t1 = tcurve(n).yearsp1;
t2 = tcurve(n).yearsp2;
% trt1 = tcurve(n).tryearsp1;
% trt2 = tcurve(n).tryearsp2;
% stt1 = tcurve(n).styearsp1;
% stt2 = tcurve(n).styearsp2;
tr1 = tcurve(n).trdatp1;
tr2 = tcurve(n).trdatp2;
st1 = tcurve(n).stdatp1;
st2 = tcurve(n).stdatp2;

%create plot (divergence test)

    figure;
    %axis([tcurve(i).min, tcurve(i).max,
    %min(tcurve(i).temp)-5,max(tcurve(i).temp)+5]);
    trplot3 = plot(t1,tr1-trmean,'g--',t1,st1,'r--',t2,tr2-trmean,'g--',t2,st2,'r--',tcurve(n).styrs,tcurve(n).stmean,'c');
    hold on
    line([1950 1950],[-250 250]);
    hold off
    
end

save('tcurve.mat','tcurve')