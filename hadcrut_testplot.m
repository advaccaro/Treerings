%hadcrut test plotting
fig('Tree-Ring MXD Availability'),clf
subplot(3,1,1:2)
m_proj('Robinson','clong',180);
m_grid('xtick',[0:60:360],'tickdir','out','ytick',[-90:30:90], 'color',dkgr, 'fontsize',10,'fontname','Times New Roman');

m_coast('color','k');

H.lon(H.lon<0) = H.lon(H.lon<0) + 360;

m_line(H.lon,H.lat,'color','b','marker','*','MarkerFaceColor','b','MarkerSize',8.5,'LineStyle','none');
