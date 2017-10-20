% The purpose of this script is to see how the realistic Signal-to-noise ratio is globally distributed. % One may think that we have pretty good SNR in general, however, in reality SNR is ~0.25 globally.
% addpath(genpath('/Users/jianghaw/Paleo/CFR_work/matlib/'));
% addpath('/Users/jianghaw/Paleo/CFR_work/data/')

load JEG_graphics

%preallocate for nproxy
nproxy1 = nan(1,length(hgroup(1).nproxy));
nproxy2 = nan(1,length(hgroup(2).nproxy));
nproxy3 = nan(1,length(hgroup(3).nproxy));
nproxy = vertcat(nproxy1, nproxy2, nproxy3);


tm = zeros(3,hgroup(1).ny);
hl = zeros(3,1); 

proxycolor = zeros(3,3);

icon{1,1}= dark_green;   icon{1,2}= '^';
icon{1,3} = sprintf('Raw MXD, n = %u', ng1); 
icon{1,4} = 'Raw MXD';
icon{2,1}= ornj;   icon{2,2}= 'v';
icon{2,3} = sprintf('MXD with station pair, n = %u', ng2);
icon{2,4} = 'MXD with station pair';
icon{3,1}= Lawn_green;   icon{3,2}= 'o';
icon{3,3} = sprintf('Divergence-free MXD, n = %u', ng3);
icon{3,4} = 'Divergence-free MXD';


fig('Tree-Ring MXD Availability'),clf
subplot(3,1,1:2)
m_proj('Robinson','clong',0);
m_grid('xtick',[-180:60:180],'tickdir','out','ytick',[-90:30:90], 'color',dkgr, 'fontsize',10,'fontname','Times New Roman');

m_coast('color','k');


for j = 1:3
    np = length(hgroup(j).nproxy);
    for i = 1:np
        nproxy(j,i) = hgroup(j).nproxy(i);
    end
    %group(j).lon(group(j).lon>180) = group(j).lon(group(j).lon>180) - 360;
    tm(j,:) = hgroup(1).tm;
end

for j = 1:3
    proxycolor(j,:) = icon{j,1};
    hl(j)=m_line(hgroup(j).lon,hgroup(j).lat,'color',icon{j,1},'marker',icon{j,2},'MarkerFaceColor',icon{j,1},'MarkerSize',8.5,'LineStyle','none');
end



[LEGH,OBJH]=legend(hl(:),icon{:,3});%pause; 
set(LEGH,'FontName','Times','FontSize',14,'Location','SouthEast');
set( findobj(OBJH,'Type','line'), 'Markersize', 12)
legend('boxon')

caxis([hgroup(1).year_i , hgroup(1).year_f]), 
colormap(proxycolor), caxis([hgroup(1).year_i,  hgroup(1).year_f]);
h = colorbar2('horiz','Most ancient age resolved');
set(h,'position',[0.13 0.27 0.775 0.04075])
title('Spatial Tree-Ring MXD Availability','FontWeight','bold','FontSize',14,'FontName','Times');



% export_fig('Proxy_by_age.pdf','-cmyk','-r300')

subplot(3,1,3)


hb = bar(tm(1,:),nproxy','grouped');
axis([1200 1980 0 600])
fancyplot_deco('Tree-Ring MXD availability over time','Time','# proxies');

[LEGH,OBJH,OUTH,OUTM]=legend(hb,icon{:,4},'Location','NorthWest'); % 4 is for Lower right-hand corner
set(LEGH,'FontName','Times','FontSize',14);
legend boxoff
%xlim([1200 2000]);
% colormap(proxycolor)
%orient landscape

pause;legend('boxoff');

orient landscape

%print -painters -dpdf -cmyk -r1000 tr_availmap1.pdf 
%export_fig('tr_availmap1.pdf','-r300','-cmyk')
%hepta_figprint('M08_proxy_availability')



