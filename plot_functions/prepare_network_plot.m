function [ x,y,u,v ] = prepare_network_plot( net_fluxes, scaling, ps, mode, fontsize, font_offset, rotation_mode, no_numbers )
    
    if isnan(scaling)
        scaling = max(max(net_fluxes));
    end
    
    if ~exist('rotation_mode','var')
        warning('Unknown rotation_mode: setting default to vertical')
        rotation_mode = 'vertical';
    end
    
    rw = 1;
    rh = 3*rw;
    dw = 1.3*rh;
    dh = dw;
    rx1 = -1.5*dw-2*rw;
    rx2 = -dw/2-rw;
    rx3 = dw/2;
    rx4 = 1.5*dw+rw;
    ry1 = dh+rh/2;
    ry2 = -rh/2;
    ry3 = -dh-rh-rh/2;
        
    if strcmp(mode,'new')  % if new plot, draw the configurations and their numbers once
        
        ylim([-9 9])
        axis equal
        box on

        if ~any(any(isnan(ps))) && size(ps,1) == 1
            % Draw probabilities inside rectangles
            fill_color = [1 0.5 0.3];
            rectangle('Position',[rx1  ry2  rw rh*ps(1)],'FaceColor',fill_color,'EdgeColor',fill_color); hold on;
            rectangle('Position',[rx2  ry1  rw rh*ps(4)],'FaceColor',fill_color,'EdgeColor',fill_color)
            rectangle('Position',[rx2  ry2  rw rh*ps(3)],'FaceColor',fill_color,'EdgeColor',fill_color)
            rectangle('Position',[rx2  ry3  rw rh*ps(2)],'FaceColor',fill_color,'EdgeColor',fill_color)
            rectangle('Position',[rx3  ry1  rw rh*ps(5)],'FaceColor',fill_color,'EdgeColor',fill_color)
            rectangle('Position',[rx3  ry2  rw rh*ps(6)],'FaceColor',fill_color,'EdgeColor',fill_color)
            rectangle('Position',[rx3  ry3  rw rh*ps(7)],'FaceColor',fill_color,'EdgeColor',fill_color)
            rectangle('Position',[rx4  ry2  rw rh*ps(8)],'FaceColor',fill_color,'EdgeColor',fill_color)
        elseif ~any(any(isnan(ps))) && size(ps,1) == 2
            % Draw probabilities inside rectangles for act and rep
            % act
            fill_color = [0.3 0.8 0.3];
            rectangle('Position',[rx1  ry2  rw/2 rh*ps(1,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01); hold on;
            rectangle('Position',[rx2  ry1  rw/2 rh*ps(4,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx2  ry2  rw/2 rh*ps(3,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx2  ry3  rw/2 rh*ps(2,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx3  ry1  rw/2 rh*ps(5,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx3  ry2  rw/2 rh*ps(6,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx3  ry3  rw/2 rh*ps(7,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx4  ry2  rw/2 rh*ps(8,1)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            % rep
            fill_color = [1 0.5 0.3];
            rectangle('Position',[rx1+rw/2  ry2  rw/2 rh*ps(1,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx2+rw/2  ry1  rw/2 rh*ps(4,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx2+rw/2  ry2  rw/2 rh*ps(3,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx2+rw/2  ry3  rw/2 rh*ps(2,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx3+rw/2  ry1  rw/2 rh*ps(5,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx3+rw/2  ry2  rw/2 rh*ps(6,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx3+rw/2  ry3  rw/2 rh*ps(7,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
            rectangle('Position',[rx4+rw/2  ry2  rw/2 rh*ps(8,2)],'FaceColor',fill_color,'EdgeColor',fill_color,'LineWidth',0.01)
        end
        
        % Draw rectangles
        rectangle('Position',[rx1  ry2  rw rh]); hold on;
        rectangle('Position',[rx2  ry1  rw rh])
        rectangle('Position',[rx2  ry2  rw rh])
        rectangle('Position',[rx2  ry3  rw rh])
        rectangle('Position',[rx3  ry1  rw rh])
        rectangle('Position',[rx3  ry2  rw rh])
        rectangle('Position',[rx3  ry3  rw rh])
        rectangle('Position',[rx4  ry2  rw rh])

        % Draw disks inside rectangles
        rcd = 0.3;  % distance between rectangle and circle
        rcr = rw - 2*rcd;
        rectangle('Position',[rx1+rcd  ry2+rcd rcr rcr],'Curvature',[1 1],'FaceColor','k')
        rectangle('Position',[rx1+rcd  ry2+rcd+rw rcr rcr],'Curvature',[1 1],'FaceColor','k')
        rectangle('Position',[rx1+rcd  ry2+rcd+2*rw rcr rcr],'Curvature',[1 1],'FaceColor','k')

        rectangle('Position',[rx2+rcd  ry1+rcd+rw rcr rcr],'Curvature',[1 1],'FaceColor','k')
        rectangle('Position',[rx2+rcd  ry1+rcd+2*rw rcr rcr],'Curvature',[1 1],'FaceColor','k')

        rectangle('Position',[rx2+rcd  ry2+rcd rcr rcr],'Curvature',[1 1],'FaceColor','k')
        rectangle('Position',[rx2+rcd  ry2+rcd+2*rw rcr rcr],'Curvature',[1 1],'FaceColor','k')

        rectangle('Position',[rx2+rcd  ry3+rcd rcr rcr],'Curvature',[1 1],'FaceColor','k')
        rectangle('Position',[rx2+rcd  ry3+rcd+rw rcr rcr],'Curvature',[1 1],'FaceColor','k')

        rectangle('Position',[rx3+rcd  ry3+rcd rcr rcr],'Curvature',[1 1],'FaceColor','k')

        rectangle('Position',[rx3+rcd  ry2+rcd+rw rcr rcr],'Curvature',[1 1],'FaceColor','k')

        rectangle('Position',[rx3+rcd  ry1+rcd+2*rw rcr rcr],'Curvature',[1 1],'FaceColor','k')
    
        if ~exist('no_numbers','var') || no_numbers == 0
            if strcmp('rotation_mode','vertical')
                sx = 0.85+font_offset;  % needs to be adapted with fontsize to move the numbers in front of the boxes
                sx2 = 0.2;
                sy = 0.0;
                text(rx1-sx,ry2+rh/2-sy,'1','FontSize',fontsize);
                text(rx2-sx,ry1+rh/2-sy,'4','FontSize',fontsize);
                text(rx2+rw+sx2,ry2+rh/2-sy,'3','FontSize',fontsize);
                text(rx2-sx,ry3+rh/2-sy,'2','FontSize',fontsize);
                text(rx3+rw+sx2,ry3+rh/2-sy,'7','FontSize',fontsize);
                text(rx3-sx,ry2+rh/2-sy,'6','FontSize',fontsize);
                text(rx3+rw+sx2,ry1+rh/2-sy,'5','FontSize',fontsize);
                text(rx4+rw+sx2,ry2+rh/2-sy,'8','FontSize',fontsize);
            else  % x and y are effectivly swapped in the later camroll(-90) command
                sx = 0.6;  % needs to be adapted with fontsize to move the numbers in front of the boxes
                sx2 = 0.5+font_offset;
                sy = 0.2;
                text(rx1-sx,ry2+rh/2-sy,'1','FontSize',fontsize);
                text(rx2-sx,ry1+rh/2-sy,'4','FontSize',fontsize);
                text(rx2+rw+sx2,ry2+rh/2-sy,'3','FontSize',fontsize);
                text(rx2-sx,ry3+rh/2-sy,'2','FontSize',fontsize);
                text(rx3+rw+sx2,ry3+rh/2-sy,'7','FontSize',fontsize);
                text(rx3-sx,ry2+rh/2-sy,'6','FontSize',fontsize);
                text(rx3+rw+sx2,ry1+rh/2-sy,'5','FontSize',fontsize);
                text(rx4+rw+sx2,ry2+rh/2-sy,'8','FontSize',fontsize);
            end
        end
        if strcmp(rotation_mode,'horizontal')
            camroll(-90)
            xlim([-9 9])
        end
    elseif strcmp(mode,'add')
    else
        error('Unknown mode')
    end
    
    %% Calculate arrows
    d = rw/2;
    l = dh-2*d;
    dg = 1/sqrt(2)*l;  % diagonal
    ds = (dh-l/sqrt(2))/2;
    
    Astates = [1,4;      1,3;     1,2;      4,5;     3,5;      3,7;      4,6;      2,6;      2,7;     5,8;      6,8;     7,8];  % the corresponding D reaction is 1->2, ...
    Dx =      [rx1+rw+ds rx1+rw+d rx1+rw+ds rx2+rw+d rx2+rw+ds rx2+rw+ds rx2+rw+ds rx2+rw+ds rx2+rw+d rx3+rw+ds rx3+rw+d rx3+rw+ds];
    Dy =      [ry2+rh+ds ry2+rh/2 ry2-ds    ry1+rh/2 ry2+rh+ds ry2-ds    ry1-ds    ry3+rh+ds ry3+rh/2 ry1-ds    ry2+rh/2 ry3+rh+ds];
    Du =      [dg        l        dg        l        dg        dg        dg        dg        l        dg        l        dg];
    Dv =      [dg        0       -dg        0        dg       -dg       -dg        dg        0       -dg        0        dg];
    %quiver(Dx,Dy,Du,Dv,0,'LineWidth',2)
    
    Dstates = [Astates(:,2) Astates(:,1)];
    Ax = Dx + Du;
    Ay = Dy + Dv;
    Au = -Du;
    Av = -Dv;
    %quiver(Ax,Ay,Au,Av,0,'LineWidth',2)
    
    Sstates = [ 3,4;     2,3;     6,5;     7,6; ];
    Sx =      [ rx2+rw/2 rx2+rw/2 rx3+rw/2 rx3+rw/2];
    Sy =      [ ry1-d    ry2-d    ry1-d    ry2-d];
    Su =      [ 0        0        0        0];
    Sv =      [-l       -l       -l       -l];
    %quiver(Sx,Sy,Su,Sv,0,'LineWidth',2)
   
    Sstates = [Sstates; Sstates(:,2) Sstates(:,1)];
    Sx2 = Sx + Su;
    Sy2 = Sy + Sv;
    Su2 = -Su;
    Sv2 = -Sv;
    Sx = [Sx Sx2];
    Sy = [Sy Sy2];
    Su = [Su Su2];
    Sv = [Sv Sv2];
    
    %quiver(Sx,Sy,Su,Sv,0,'LineWidth',2)
    
    Rx = NaN*ones(8);
    Ry = NaN*ones(8);
    Rv = NaN*ones(8);
    Ru = NaN*ones(8);  
    for i=1:size(Astates,1)
        Rx(Astates(i,1),Astates(i,2)) = Ax(i);
        Ry(Astates(i,1),Astates(i,2)) = Ay(i);
        Ru(Astates(i,1),Astates(i,2)) = Au(i);
        Rv(Astates(i,1),Astates(i,2)) = Av(i);
        Rx(Dstates(i,1),Dstates(i,2)) = Dx(i);
        Ry(Dstates(i,1),Dstates(i,2)) = Dy(i);
        Ru(Dstates(i,1),Dstates(i,2)) = Du(i);
        Rv(Dstates(i,1),Dstates(i,2)) = Dv(i);
    end
    for i=1:size(Sstates,1)
        Rx(Sstates(i,1),Sstates(i,2)) = Sx(i);
        Ry(Sstates(i,1),Sstates(i,2)) = Sy(i);
        Ru(Sstates(i,1),Sstates(i,2)) = Su(i);
        Rv(Sstates(i,1),Sstates(i,2)) = Sv(i);
    end
    x = Rx(net_fluxes>0);
    y = Ry(net_fluxes>0);
    u = Ru(net_fluxes>0).*net_fluxes(net_fluxes>0)./scaling;
    v = Rv(net_fluxes>0).*net_fluxes(net_fluxes>0)./scaling;
     

    
end
