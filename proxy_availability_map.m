% The purpose of this script is to see how the realistic Signal-to-noise ratio is globally distributed. % One may think that we have pretty good SNR in general, however, in reality SNR is ~0.25 globally.
% addpath(genpath('/Users/jianghaw/Paleo/CFR_work/matlib/'));
% addpath('/Users/jianghaw/Paleo/CFR_work/data/')
% b = load('pproxy_mann08_sparse_realisticSNR.mat');
% ptype = b.ptype; np = size(ptype,2);
% proxy = b.pproxy{2}(:,:,1);
% pcode = zeros(1,np);
% lon = b.plon; lat = b.plat; % coordinates of the real proxies
% time = b.ptime; nt = length(time);
% 
load JEG_graphics
nproxy = zeros(nt,8);
pcode(round(ptype/1000) == 9) = 1;
pcode(round(ptype/1000) == 8) = 3;
pcode(round(ptype/1000) == 7) = 4;
pcode(round(ptype/100) == 75) = 2;
pcode(round(ptype/1000) == 6) = 5;
pcode(round(ptype/1000) == 5) = 6;
pcode(round(ptype/1000) == 4) = 7;
pcode(round(ptype/1000) == 3 | round(ptype/1000) == 2) = 8;
proxycolor = zeros(8,3);

icon{1,1}=bright_green; icon{1,2}= '^'; icon{1,3} = 'Tree Ring Width';
icon{2,1}=Lawn_green;   icon{2,2}= 'v'; icon{2,3} = 'Tree Ring MXD';
icon{3,1}=skyblue;      icon{3,2}= '*'; icon{3,3} = 'Ice core';
icon{4,1}=ornj;         icon{4,2}= 'o'; icon{4,3} = 'Coral';
icon{5,1}=hotpink;      icon{5,2}= 'd'; icon{5,3} = 'Speleothem';
icon{6,1}=dkgr;         icon{6,2}= 'p'; icon{6,3} = 'Documentary';
icon{7,1}=maroon;       icon{7,2}= 's'; icon{7,3} = 'Sediment';
icon{8,1}=blue;         icon{8,2}= 'h'; icon{8,3} = 'Composite';

c=greengiant;
hla=zeros(n,1);
Nc = size(greengiant,1);
n = length(pcode);
for k = 1:n
    ind = find(~isnan(proxy(:,k)),1,'first');
    tmin(k) = time(ind);
end
tm = linspace(min(tmin),max(tmin),Nc);

fig('Proxy Availability'),clf
 subplot(3,1,1:2)
m_proj('Robinson','clong',180);
m_grid('xtick',[0:60:360],'tickdir','out','ytick',[-90:30:90], 'color',dkgr, 'fontsize',8,'fontname','Times New Roman');

m_coast('color','k');

for k = 1:8
    nproxy(:,k) = sum(~isnan(proxy(:,pcode == k)),2);
    ind(k)= find(pcode == k, 1, 'last' );  

end

lon(lon<0) = lon(lon<0) + 360;
for k = 1:n
    [~,imin]=min(abs(tm(:)-tmin(k)));
    
    hl(k)=m_line(lon(k),lat(k),'color',c(imin,:),'marker',icon{pcode(k),2},'MarkerFaceColor',c(imin,:),'MarkerSize',7,'linestyle','none');
end
caxis([min(tm) , max(tm)]), 
colormap(greengiant), caxis([min(tm),  max(tm)]);
h = colorbar2('horiz','Most ancient age resolved');
set(h,'position',[0.13 0.27 0.775 0.04075])
title('Proxy in the M08 Network','FontWeight','bold','FontSize',14,'FontName','Times');

[LEGH,OBJH,OUTH,OUTM]=legend(hl(ind),icon{:,3});pause; % 4 is for Lower right-hand corner
set(LEGH,'FontName','Times','FontSize',9);
set(OUTH,'Color','k','MarkerFaceColor','k','Markersize',7);
legend boxoff

export_fig('Proxy_by_age.pdf','-cmyk','-r300')

 subplot(3,1,3)
for i = 1:8
    proxycolor(i,:) = icon{i,1};
end
hb = bar(time,nproxy,'stacked');
axis([850 2000 0 1200])
fancyplot_deco('Proxy availability in M08 Network (1138 records)','Time','# proxies');

[LEGH,OBJH,OUTH,OUTM]=legend(hb,icon{:,3},'Location','NorthWest'); % 4 is for Lower right-hand corner
set(LEGH,'FontName','Times','FontSize',9);
legend boxoff
colormap(proxycolor)
orient landscape

export_fig('temporal_proxy_avail.pdf','-r300','-cmyk')
%hepta_figprint('M08_proxy_availability')



