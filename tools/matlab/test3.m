clc
clear all
close all


H.Fig = figure('position',[100,100,820,800],'color','w'); % Position in pixels: [x_links,y_onder,breedte,hoogte]
H.Ax( 1,1) = axes('parent',H.Fig,'units','pixels','position',[ 80,450,320,325]);
H.Ax( 1,2) = axes('parent',H.Fig,'units','pixels','position',[480,450,320,325]);
H.Ax( 2,1) = axes('parent',H.Fig,'units','pixels','position',[ 80, 50,320,325]);
H.Ax( 2,2) = axes('parent',H.Fig,'units','pixels','position',[480, 50,320,325]);

for i = 1: 2
  for j = 1: 2
    set( H.Ax( i,j),'fontsize',24,'xgrid','on','ygrid','on')
    if (i == 2)
      xlabel( H.Ax( i,j),'x-coordinate (m)')
    end
    if (j == 1)
      ylabel( H.Ax( i,j),'y-coordinate (m)')
    end
  end
end