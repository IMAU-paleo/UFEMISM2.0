clc
clear all
close all

filename = '../../results_20230605_001/test_remapping_mesh2grid_output.nc';

wa = 200;
ha = 200;
margins_hor = [25, 25, 25];
margins_ver = [25, 25];
H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: 1
  for j = 1: 2
    set( H.Ax{ i,j},'xtick',[],'ytick',[],'clim',[-1,1],'xlim',[0,1],'ylim',[0,1]);
  end
end

d_grid_ex = ncread( filename,'d_grid_ex');
d_grid    = ncread( filename,'d_grid');

image('parent',H.Ax{ 1,1},'xdata',[0,1],'ydata',[0,1],'cdata',d_grid_ex','cdatamapping','scaled');
image('parent',H.Ax{ 1,2},'xdata',[0,1],'ydata',[0,1],'cdata',d_grid'  ,'cdatamapping','scaled');