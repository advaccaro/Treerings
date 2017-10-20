load('hgroup.mat');
load JEG_graphics
% load proxy files
% ==========================
% Make the color figures for MXD/station plots, 6 to a page
nt = length(hcurve);
m = 1;
for i = 1:nt
    if hcurve(i).smtestcheck == 4
        hc(m).tr = hcurve(i).trsmooth;
        hc(m).tryrs = hcurve(i).tryears;
        hc(m).st = hcurve(i).stsmooth;
        hc(m).styrs = hcurve(i).styears;
        
        m = m + 1;
    end
end

nh = length(hc);


for n=1:6:nh
    if n+5 < nh
        fig(['Proxy Records ' num2str(n:n+5)]), clf
    else
        fig(['Proxy Records ' num2str(n:nh)]), clf
    end
    
    %orient landscape
    %k=p;
    for k=n:n+5
        if k<=nh
            % S(k).Site
            subplot(3,2,k-n+1)
            trmean = nmean(hc(k).tr);

% t1 = hcurve(n).yearsp1;
% t2 = hcurve(n).yearsp2;
% trt1 = tcurve(n).tryearsp1;
% trt2 = tcurve(n).tryearsp2;
% stt1 = tcurve(n).styearsp1;
% stt2 = tcurve(n).styearsp2;
% tr1 = hcurve(n).trdatp1;
% tr2 = hcurve(n).trdatp2;
% st1 = hcurve(n).stdatp1*100;
% st2 = hcurve(n).stdatp2*100;

t1 = hc(k).styrs;
t2 = hc(k).tryrs;
tr = hc(k).tr;
st = hc(k).st*100;
tmin = max([min(hc(k).styrs), min(hc(k).tryrs)]);
tmax = min([max(hc(k).styrs), max(hc(k).tryrs)]);

%create plot (divergence test)

    %axis([tcurve(i).min, tcurve(i).max,
    %min(tcurve(i).temp)-5,max(tcurve(i).temp)+5]);
plot(t1, st,'r', t2, tr-trmean,'b');
%     hold on
%     line([1950 1950],[-250 250]);
%     hold off
    
            % tmin = max(1000,S(k).Initial_year);
%             xlim([tmin tmax]);
%             ylabel(char(S(k).Type),'FontName','Times','FontSize',[12])
%             xlabel('Time','FontName','Times','FontSize',[12])
%             % title
% 	   % title([Type(k),unit(k),'in',Site(k),'(',Lat(k),'^\circ',Lon(k),'\circ',')']);
%             rec=[num2str(k) ')'];
%             obj=char(S(k).Object);
%             Type=char(S(k).Type);
%             Lon = S(k).Lon;
%             if Lon<180
%                 Lon_str=[num2str(round(S(k).Lon)),'{}^{\circ}E'];
%             else
%                 Lon_str=[num2str(360-round(S(k).Lon)),'{}^{\circ}W'];
%             end
%             Lat=S(k).Lat;
%             if (Lat>=0)
%                 Lat_str=[num2str(round(Lat)),'{}^{\circ}N'];
%             else
%                 Lat_str=[num2str(-round(Lat)),'{}^{\circ}S'];
%             end
%             str=[rec,' ', obj,' ',Type,' in ', deblank(Sites(k,:)),' (',Lon_str,',',Lat_str,')'];
%             title(str,'FontName','Times','FontSize',[14])
        end
    end
    % filen=['./figures_db/coral_Database_', num2str(n),'_to_',num2str(n+5),'.pdf'];
    %pause
    % export_fig(filen,'-cmyk','-r300');
end
