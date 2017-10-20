close all;  clear all;
load('/Users/jeg/Documents/Science/Research/LM_synthesis/data/sensor/coral_db/coral_db_tropical_year.mat'); 
load JEG_graphics
% Select only published records
tier = [S.tier]';
S0 = S;
S = S0(find(tier == 1)); % lost 9 records

res = [S.res]';
x = [1:12];
fig('Resolution Histogram')
hist(res,x);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',ligr)
fancyplot_deco('Temporal resolution in Coral Database','\Delta t (months)','Count');
export_fig('./figs/coral_db_resolution.pdf','-r200');

% Make the color figures for coral data anomaly plots, 6 to a page
nr=length(S);
sites=char({S.site});
%
for n=1:6:nr
    if n+5 < nr
        fig(['Coral Records' num2str(n:n+5)]), clf
    else
        fig(['Coral Records' num2str(n:nr)]), clf
    end
    
    %orient landscape
    %k=p;
    for k=n:n+5
        if k<=nr
            S(k).site
            subplot(3,2,k-n+1)
            tr=S(k).chron_reg;
            Xa=S(k).anom;
            dt = res(k);
            plot(tr,Xa,'r-'), grid; axis tight;
            tmin = max(1000,S(k).year_i);
            xlim([tmin S(k).year_f]);
            ylabel(char(S(k).type),'FontName','Times','FontSize',[10])
            xlabel('Time','FontName','Times','FontSize',[10])
            % title
            rec=[num2str(k) ')'];
            obj=char(S(k).class)
            type=char(S(k).type);
            lon_str=[num2str(round(S(k).lon)),'{}^{\circ}E'];
            lat=S(k).lat;
            if (lat>=0)
                lat_str=[num2str(round(lat)),'{}^{\circ}N']
            else
                lat_str=[num2str(-round(lat)),'{}^{\circ}S']
            end
            str=[rec,': ', obj,' ',type,' in ', deblank(sites(k,:)),' (',lon_str,',',lat_str,'), \Deltat = ',int2str(dt)];
            title(str,'FontName','Times','FontSize',[10])
        end
    end
filen=['./figs/coral_database_', num2str(n),'_to_',num2str(n+5),'.pdf'];
%pause
export_fig(filen,'-cmyk','-r300');
end

% form coral matrix
year_min=[S.year_i];
year_max=[S.year_f];

%create colormap
load('greengiant.mat')
c=greengiant;
hla=zeros(nr,1);
Nc = size(greengiant,1);
tm = linspace(min(year_min),max(year_min),Nc);

tr = [min(year_min):max(year_max)];
nt = length(tr); % number of observations
coral_mat = NaN(nt,nr);

for k=1:nr 
    ind = find(tr <= year_max(k) & tr >= year_min(k));
    coral_mat(ind,k) = 1; 
end

% assign colors to proxy types
prox_col(1,:) = ligr;
prox_col(2,:) = 1.2*dkgr;
prox_col(3,:) = bck;

%==========================================================
%  Now by age	
%==========================================================
c_lon = [S(:).lon]';
c_lat = [S(:).lat]';
type = {S(:).type}';

% fix longitudes:
c_lon(c_lon <0) = 360 + c_lon(c_lon <0);
%
fig('Database by Age'),clf
subplot(3,1,1:2)
% prepare canvas
orient landscape
m_proj('Robinson','clong',170,'lat',[-40 40]); 
m_coast('patch',ligr);
m_grid('box','off','xtick',6,'ytick',9,'xlabeldir','end','xticklabels',[], 'fontsize',5,'fontname','Times');

for k=1:nr
    if strcmp(type(k),'\delta^{18}O')
        mark='o'; pcode(k) = 1; markerc=firebrick;
    elseif strcmp(type(k),'Sr/Ca');
        mark='s'; pcode(k) = 2; markerc='k';
    else
        mark = '^'; pcode(k) = 3;markerc=firebrick;
    end
	% find closest age point
	[~,imin]=min(abs(tm(:)-year_min(k)));
	hla(k)=m_line(c_lon(k),c_lat(k),'marker',mark,'MarkerEdgeColor',markerc,'MarkerFaceColor',c(imin,:),'Color',c(imin,:),'linewidth',[1],'MarkerSize',[12],'linestyle','none');
    %
    if strcmp(type(k),'\delta^{18}O')
        h1 = hla(k); kd = k;
    elseif strcmp(type(k),'Sr/Ca');
        h2 = hla(k); ks = k;
    else
        h3 = hla(k); kf = k;
    end
end
% colorbar
caxis([min(tm) , max(tm)]), 
colormap(greengiant), caxis([min(tm),  max(tm)]);
colorbar2('horiz','Most ancient age resolved');

% legend
h = [h1 h2 h3]; lab = [type(kd) type(ks) 'other'];
[LEGH,OBJH,OUTH,OUTM] = legend(h ,lab{:});

set(LEGH,style_l{:},'location','Southeast'), pause
legend boxoff
title(['a) Proxy representation by age  (',num2str(nr), ' records)'],style_t{:});
%
subplot(3,1,3)
for j = 1:3
   ncoral(:,j) = nsum(coral_mat(:,pcode==j),2);
end

hb=bar(tr,sq(ncoral),'stacked')
axis([1500 2000 0 nr])
fancyplot_deco('b) Proxy availability over time','Time','# proxies');
%Set colors
for k = 1:3
    set(hb(k),'EdgeColor',prox_col(k,:),'FaceColor',prox_col(k,:));
end

[LEGH,OBJH,OUTH,OUTM]=legend(hb,lab{:},'location','NorthWest'); % 4 is for Lower right-hand corner
set(LEGH,style_l{:},'box','off'); 

hepta_figprint('./figs/coral_db_map_by_age');

% ============================================
%  MATRIX PLOT
% ============================================
year_i = [S(:).year_i]; % update age
year_f = [S(:).year_f];
tc = [min(year_i):max(year_f)];

tcmin = tc(find(ncoral >= 5,1,'first'));
tm = [tcmin:1995]; % too few records beforehand


%  Select DJF season 
for r = 1:nr
   d = S(r).anom;
   t = S(r).chron_reg;
   if res(r) < 12
      if res(r) <= 6 & res(r) > 3  % if subannual, obtain cold season averages
         ms = 10; me = 3;
      elseif res(r) <= 3
         ms = 12; me = 2;
      end
      [da,ta,ts]=intra_annual_avg2(d,t,ms,me,1);
      %plot(t,d,'k-',ta,da,'rx');
      %fancyplot_deco(S(r).site, 'Time', S(r).type);
      S(r).data_annual = da;
      S(r).chron_annual = ta;
   elseif res(r) == 12  
      S(r).data_annual = d;
      S(r).chron_annual = round(t); %+1 ?
   end
end

% put all records in one matrix 
nm = length(tm); % number of years
coral = nan(nm,nr);
sitesf = {S(:).site};

for r = 1:nr
	t = S(r).chron_annual;
	d = S(r).data_annual;	
   [tt,dc,ind] = consolidator(t,d);
   ti = intersect(tt,tm);
	coral(ismember(tm,ti),r)= dc(ismember(tt,ti));
	clear d t 	
   if strcmpi(S(r).type,'\delta^{18}O') | strcmpi(S(r).type,'km^3');
      plab{r} = ['$', char(S(r).type), '$'];
   else
      plab{r} = char(S(r).type);
   end
   llon{r} = [sprintf('%d',round(c_lon(r))),'^{\circ} E']; 
end
% sort by longitude
%[~,ilon] = sort(c_lon);
ilon = [1:nr];
v = [1:nr]; ind = find(tm >= 1500);
fig('Proxy matrix'),clf
pcolor(tm(ind), v, standardize(coral(ind,ilon))'), shading flat
set(gca,'Ytick',v,'YtickLabel',plab(ilon),'FontSize',13); plotTickLatex2D('xlabeldy',0.00);
hTitle=title('Proxy Matrix View');
hYLabel=ylabel([]); hXLabel=xlabel('Time');
FontSize = 16;
set(gca,'FontName','Times','FontSize', round(FontSize*0.71));
set([hTitle, hXLabel, hYLabel], 'FontName' , 'Palatino');
set([hXLabel, hYLabel],'FontSize', round(FontSize*0.86));
set(hTitle, 'FontSize', FontSize, 'FontWeight' , 'bold');
% plot site name
lim = get(gca,'Xlim'); offset = 0.02;
text(repmat(lim(2)+offset*(lim(2)-lim(1)),length(v),1),v, ...
  sitesf,'VerticalAlignment','middle',...
  'HorizontalAlignment','left','rotation',45);
% text(repmat(lim(2)+offset*(lim(2)-lim(1)),length(v),1),v, ...
%    llon(ilon),'VerticalAlignment','middle',...
%    'HorizontalAlignment','left','rotation',45);
colorbar2('horiz')
% Print 
hepta_figprint('./figs/coral_matrix');

%
sites=char({S.site});
for r = 1:nr
  if (~isempty(strfind(ddeblank(sites(r,:)),'Mahe')))
     imahe = r;
  end
end
