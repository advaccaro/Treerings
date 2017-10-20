%
close all;  clear all;
% addpath(genpath('/Users/jianghaw/Paleo/CFR_work/matlib'))

load('global_proxy_database_pre1800.mat');
load JEG_graphics
% load proxy files
% ==========================
% Step 1 
% Make the color figures for coral Data anomaly plots, 4 to a page
nr=length(S);
Sites=char({S.Site});
%
for n=1:6:nr
    if n+5 < nr
        fig(['Proxy Records ' num2str(n:n+5)]), clf
    else
        fig(['Proxy Records ' num2str(n:nr)]), clf
    end
    
    %orient landscape
    %k=p;
    for k=n:n+5
        if k<=nr
            % S(k).Site
            subplot(3,2,k-n+1)
            trmean = nmean(hcurve(n).trtemp);

            t1 = hcurve(k).yearsp1;
            t2 = hcurve(k).yearsp2;
            % trt1 = tcurve(n).tryearsp1;
            % trt2 = tcurve(n).tryearsp2;
            % stt1 = tcurve(n).styearsp1;
            % stt2 = tcurve(n).styearsp2;
            tr1 = hcurve(k).trdatp1;
            tr2 = hcurve(k).trdatp2;
            st1 = hcurve(k).stdatp1*100;
            st2 = hcurve(k).stdatp2*100;


            figure;
            %axis([tcurve(i).min, tcurve(i).max,
            %min(tcurve(i).temp)-5,max(tcurve(i).temp)+5]);
            plot(t1,tr1-trmean,'g',t1,st1,'r',t2,tr2-trmean,'g--',t2,st2,'r--');
            hold on
            line([1950 1950],[-250 250]);
            hold off
    

            % tmin = max(1000,S(k).Initial_year);
            xlim([hcurve(k).tmin hcurve(k).tmax]);
            %ylabel(char(S(k).Type),'FontName','Times','FontSize',[12])
            %xlabel('Time','FontName','Times','FontSize',[12])
            % title
	   % title([Type(k),unit(k),'in',Site(k),'(',Lat(k),'^\circ',Lon(k),'\circ',')']);
            %rec=[num2str(k) ')'];
            %obj=char(S(k).Object);
            %Type=char(S(k).Type);
            %Lon = S(k).Lon;
            %if Lon<180
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

% ================================
% Step 2: Create Color Map

% form coral matrix
year_min=[S.Initial_year];
year_max=[S.Final_year];

%create colormap
% load('greengiant.mat')
c = jet;
hla = zeros(nr,1);
Nc  = size(c,1);
tm  = linspace(0,max(year_min),Nc);

tr = [0:max(year_max)];
nt = length(tr); % number of observations
proxy_mat = NaN(nt,nr);

% CalcuLate the age of each proxy
for k=1:nr
    ind = find(tr <= year_max(k) & tr >= year_min(k));
    proxy_mat(ind,k) = 1;
end

% assign colors to proxy types
prox_col(1,:) = bright_green;
prox_col(2,:) = skyblue;
prox_col(3,:) = hotpink;
prox_col(4,:) = maroon;
prox_col(5,:) = mid_blue;
% Define proxy symbols for each Type
icon{1}= '^'; % 'Tree Ring';
icon{2}= '*'; % 'Ice core';
icon{5}= 'o'; %icon{4,2} = 'Coral';
icon{3}= 'd'; % 'Speleothem';
icon{4}= 's'; % 'Sediment';

%==========================================================
%  Now by age
%==========================================================
c_Lon = [S(:).Lon]';
c_Lat = [S(:).Lat]';
Object = {S(:).Object}';

%
fig('Database by Age'),clf
subplot(3,1,1:2)
% prepare canvas
orient landscape
m_proj('Robinson','clong',200);
m_coast('patch',ligr);

m_grid('box','off','xtick',6,'ytick',9,'xlabeldir','end','xticklabels',[], 'fontsize',5,'fontname','Times');


hl=zeros(nr,1);
%
for k=1:nr
    if strcmp(Object(k),'TREE')
        col = prox_col(1,:); mark = icon{1}; pcode(k) = 1; kt = k; % markerc=firebrick;
    elseif strcmp(Object(k),'ICE');
        col = prox_col(2,:); mark = icon{2}; pcode(k) = 2; ki = k;% markerc='k';
    elseif strcmp(Object(k),'SPELEO')
        col = prox_col(3,:); mark = icon{3}; pcode(k) = 3; kl = k;% markerc=deep_sky;
    elseif strcmp(Object(k),'SEDI')
        col = prox_col(4,:); mark = icon{4}; pcode(k) = 4; ks = k;
    elseif strcmp(Object(k), 'CORAL')
        col = prox_col(5,:); mark = icon{5}; pcode(k) = 5; kc = k; 
    end
    
    Lonp=S(k).Lon;
    Latp=S(k).Lat;
    hl(k)=m_line(Lonp,Latp,'marker',mark,'color',col,'MarkerFaceColor',col,'linewidth',[1],'MarkerSize',[5],'linestyle','none');
end

hr=[hl(kt) hl(ki) hl(kl) hl(ks) hl(kc)];

title(['Recent contributed Proxies: ',num2str(nr), ' records'],style_t{:});

lab{1}='Tree Ring';
lab{2}='Ice Core';
lab{3}='Speleothem';
lab{4}='Sediment core';
lab{5}='Coral';

[LEGH,OBJH,OUTH,OUTM]=legend(hr,lab{:}); pause% 4 is for Lower right-hand corner
set(LEGH,style_l{:});
%  print
legend boxoff

subplot(3,1,3)
for j = 1:5
    nproxy(:,j) = nsum(proxy_mat(:,pcode==j),2);
end
hb=bar(tr,sq(nproxy),'stacked');
axis([0 2000 0 nr])
fancyplot_deco('Proxy availability over time','Time','# proxies');

%Set colors
for k = 1:5
    set(hb(k),'EdgeColor',prox_col(k,:),'FaceColor',prox_col(k,:));
end
[LEGH,~,~,OUTM]=legend(hb,lab{:},'location','NorthWest'); % 4 is for Lower right-hand corner
set(LEGH,style_l{:},'box','off');

%orient landscape
% \hepta_figprint('./figures_db/coral_db_map_by_age');
%hepta_figprint('./figs/Proxy_Database_map_by_age');

%fig('Database in color'),clf

orient landscape
% hepta_figprint('./figs/global_proxy_Database_map');
export_fig('./figs/global_proxy_Database_map.pdf','-cmyk','-r200')
