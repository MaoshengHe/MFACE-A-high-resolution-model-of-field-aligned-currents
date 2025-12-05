clear 
load  exampleData.mat
doy =0*Coor_2dGrid.MLT+(T-datenum(year(T),1,0));
imfBx = 0*Coor_2dGrid.MLT+INDEX(1);
imfBy = 0*Coor_2dGrid.MLT+INDEX(2);
imfBz = 0*Coor_2dGrid.MLT+INDEX(3);
Vsw = 0*Coor_2dGrid.MLT+INDEX(4);
AE= 0*Coor_2dGrid.MLT+INDEX(5);
%--% run MFACE
[J, ~, EOF1,~,~] = MFACE_v1(Coor_2dGrid.MLT,Coor_2dGrid.MLAT,'doy',doy,'imfBy',imfBy,...
    'imfBz',imfBz,'imfBx',imfBx,'Vsw',Vsw,'AE',AE); %,'imfBx',imfBx
%--% FAC map
Handles.J=pcolor(Coor_2dGrid.xgrid,Coor_2dGrid.ygrid,J);
set(Handles.J,'LineStyle','none');hold on;
%--% Coor_2dGrid.latitude grid
Handles.plotGrid=plot(grid_show.x,grid_show.y,'k-');
%--% coast map
Handles.coast=plot(Maps.x,Maps.y,'k');
%--% set axis
set(gca,'xlim',[-1 1]*max(Coor_2dGrid.xgrid(:)))%prctile(xxx(:),[0.01,99.99])
set(gca,'ylim',[-1 1]*max(Coor_2dGrid.ygrid(:)))%prctile(yyy(:),[0.01,99.99])
set(gca,'clim',[-1 1]*prctile(abs(J(:)),[99.95]))%
set(gca,'visible','off');    set(gca,'color','none'),
%--% missing drivers
if any(isnan(INDEX)),
    set(gca,'color','y'),
    Handles.missing_input=Text1234(gca,2,[sprintf('\n'),sprintf('\n'),'Missing Driver'],14);
    set(Handles.missing_input,'color','r')
end
axis equal;    %axis tight;
%--% Label Local Time
for ilt=1:length(LT_label.timeNum)
    Handles.LT{ilt}=text(LT_label.x(ilt),LT_label.y(ilt),0,datestr(LT_label.timeNum(ilt),'HH:MM'),...
        'Color','m',   'FontSize',11,'VerticalAlignment','middle','HorizontalAlignment','center');
end
 colormap(COLORMAP);
%--% colorbar
Handles.colorbar=colorbar(gca,'south');
set(Handles.colorbar,'Position', [0.15 0.15 0.1879 0.0156],'xcolor','m','tickdir','out')
xlabel(Handles.colorbar,'J_{Downward} (\muAm^-^2)', 'FontSize',11,...
    'Units','normalized','Position',[1.1 1.2 0], 'HorizontalAlignment','left','Color','m');
%--% print time
Handles.time= text('Parent',gca,'Units','normalized',...
    'VerticalAlignment','top',    'HorizontalAlignment','left',    'FontSize',12,...
    'String',['FAC predicted by MFACE at ',sprintf('\n'),datestr(T,31)],...
    'Position',[0 1 0]);
set(Handles.time,'color','m')
set(gcf);%[232   20   560   488]11.1041*.5
set(gcf,'Units','centimeters','Position',[0 0 13 13])
set(gcf,'PaperPositionMode','manual','PaperUnits',...
    get(gcf,'Units'),'PaperPosition',get(gcf,'Position'),'Resize','on','InvertHardcopy','off','color','w','clipping','off'); 
print('-dpng','-r900','MFACE_Example.png');