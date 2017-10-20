% The purpose of this script is to see how the realistic Signal-to-noise ratio is globally distributed. % One may think that we have pretty good SNR in general, however, in reality SNR is ~0.25 globally.
% addpath(genpath('/Users/jianghaw/Paleo/CFR_work/matlib/'));
% addpath('/Users/jianghaw/Paleo/CFR_work/data/')

load JEG_graphics

%preallocate for nproxy
nproxy1 = nan(1,length(group(1).nproxy));
nproxy2 = nan(1,length(group(2).nproxy));
nproxy3 = nan(1,length(group(3).nproxy));
nproxy = vertcat(nproxy1, nproxy2, nproxy3);
%nproxy = vertcat(nproxy1,nproxy2,nproxy3);

tm = zeros(3,group(1).ny);
hl = zeros(3,1); 

proxycolor = zeros(3,3);

icon{1,1}= [1 0 0];   icon{1,2}= '^'; icon{1,3} = 'Raw Data';
icon{2,1}= [1 1 0];   icon{2,2}= 'v'; icon{2,3} = 'Quality Controlled Data';
icon{3,1}= [0 1 0];   icon{3,2}= '*'; icon{3,3} = 'Screened Data';


fig('Tree-Ring MXD Availability'),clf
ax1 = subplot(3,1,1:2);
hold on;
m_proj('Robinson','clong',180);
m_grid('xtick',[0:60:360],'tickdir','out','ytick',[-90:30:90], 'color',dkgr, 'fontsize',8,'fontname','Times New Roman');

m_coast('color','k');

for j = 1:3 %change to 3
    np = length(group(j).nproxy);
    for i = 1:np
    nproxy(j,i) = group(j).nproxy(i);
    %ind(j)= find(pcode == j, 1, 'last' );  
    end
    
     group(j).lon(group(j).lon<0) = group(j).lon(group(j).lon<0) + 360;
    
    tm(j,:) = group(1).tm;
end

for j = 1:3
    %[~,imin]=min(abs(tm(:)-tmin(j)));
    
    proxycolor(j,:) = icon{j,1};
    hl(j)=m_line(group(j).lon,group(j).lat,'color',icon{j,1},'marker',icon{j,2},'MarkerFaceColor',icon{j,1},'MarkerSize',7,'LineStyle','none');
end
hold on;
caxis([group(1).year_i , group(1).year_f]), 
colormap(proxycolor), caxis([group(1).year_i,  group(1).year_f]);
hold on;
h = colorbar2('horiz','Most ancient age resolved');
set(h,'position',[0.13 0.27 0.775 0.04075])
title('Tree-Ring MXD Data','FontWeight','bold','FontSize',14,'FontName','Times');

[LEGH,OBJH,OUTH,OUTM]=legend(hl(:),icon{:,3});pause; % 4 is for Lower right-hand corner
set(LEGH,'FontName','Times','FontSize',10);
set(OUTH,'Markersize',9);
legend boxoff


%export_fig('Proxy_by_age.pdf','-cmyk','-r300')

subplot(3,1,3);

hb = bar(tm(1,:),nproxy','grouped');
axis([300 2000 0 1200])
fancyplot_deco('Global Tree Ring MXD availability (512 records)','Time','# proxies');

[LEGH,OBJH,OUTH,OUTM]=legend(hb,icon{:,3},'Location','NorthWest'); % 4 is for Lower right-hand corner
set(LEGH,'FontName','Times','FontSize',9);
legend boxoff
% colormap(proxycolor)
orient landscape


%export_fig('temporal_proxy_avail.pdf','-r300','-cmyk')
%hepta_figprint('M08_proxy_availability')



