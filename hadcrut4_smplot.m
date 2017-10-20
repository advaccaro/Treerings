function [trplot] = hadcrut4_plot(data,n)
%this script only works after tr_screen_step1 and step2 have been used
%load('hcurve.mat');
hcurve = data;

if hcurve(n).smtestcheck < 4
    display('record failed screening')
    display(n)
    display(hcurve(n).id)
else

trmean = nmean(hcurve(n).trsm)*100;

t1 = hcurve(n).yearsp1;
t2 = hcurve(n).yearsp2;
% trt1 = tcurve(n).tryearsp1;
% trt2 = tcurve(n).tryearsp2;
% stt1 = tcurve(n).styearsp1;
% stt2 = tcurve(n).styearsp2;
tr1 = hcurve(n).trsm1*100;
tr2 = hcurve(n).trsm2*100;
st1 = hcurve(n).stsm1*100;
st2 = hcurve(n).stsm2*100;

%create plot

    figure;
    %axis([tcurve(i).min, tcurve(i).max,
    %min(tcurve(i).temp)-5,max(tcurve(i).temp)+5]);
    trplot = plot(t1,tr1-trmean,'g',t1,st1,'r',t2,tr2-trmean,'g--',t2,st2,'r--');
    hold on
    line([1960 1960],[-250 250]);
    hold off
    
end
    
