%correlation screen
load('tcurve.mat');
nt = length(tcurve);
n=0; m = 0;%set counters =0 

for i = 1:nt
    if (tcurve(i).finalscreen == 7)
        tcurve(i).dcorr = abs(tcurve(i).Rp1 - tcurve(i).Rp2);
    else
        tcurve(i).dcorr = NaN;
        %display(i)
        n=n+1;
    end
end
 for i = 1:nt
     %if all(1/2 < tcurve(i).ratiocorr & tcurve(i).ratiocorr < 2/1)
     if tcurve(i).dcorr <.1
         tcurve(i).corrscreen = 1;
         display(i)
         m = m+1;
     else
         tcurve(i).corrscreen = 0;
     end
 end
%sprintf('%d records failed quality controll' ,n)
sprintf('%d records passed screen', m)

save('tcurve.mat', 'tcurve')