function draw_arrow(x_0, y_0, alpha, d, mycolor)

    if ~exist('mycolor','var')
        mycolor = [0 0 0];
    end

    for i=1:7
        xArrow = [ 0     0.8 0.8     1 0.8     0.8 0];                   % X-coords of arrow edge points
        yArrow = [-d/2 -d/2 -d/2-0.1 0 d/2+0.1 d/2 d/2];   % Y-coords of arrow edge points
        xArrow_rot(i) = cos(alpha)*xArrow(i) - sin(alpha)*yArrow(i) + x_0;
        yArrow_rot(i) = sin(alpha)*xArrow(i) + cos(alpha)*yArrow(i) + y_0;
    end
    %hArrow = fill(xArrow_rot,yArrow_rot, [0 0 0]);    % Plot a red arrow (leads to wierd separation lines in pdf print)
    patch(xArrow_rot([1 2 6]),yArrow_rot([1 2 6]),mycolor,'EdgeColor',mycolor,'LineWidth',0.01)
    patch(xArrow_rot([1 6 7]),yArrow_rot([1 6 7]),mycolor,'EdgeColor',mycolor,'LineWidth',0.01)
    patch(xArrow_rot([3 4 5]),yArrow_rot([3 4 5]),mycolor,'EdgeColor',mycolor,'LineWidth',0.01)
    
end