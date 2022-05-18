plot3(1,1,1,'o');
hold on;
for i = 1:max(size(grid))
        coord = grid(i).coord;
        diff = coord(2,1)-coord(1,1);
        vertices = [coord(1,:);coord(1,1)+diff coord(1,2) coord(1,3);...
                    coord(1,1) diff+coord(1,2) coord(1,3);coord(1,1) coord(1,2) diff+coord(1,3);...
                    coord(2,:);coord(2,1)-diff coord(2,2) coord(2,3);...
                    coord(2,1) coord(2,2)-diff coord(2,3);coord(2,1) coord(2,2) coord(2,3)-diff;];
        plot3(vertices(:,1),vertices(:,2),vertices(:,3),'*r');
        plot3([vertices(1,1),vertices(2,1)],[vertices(1,2),vertices(2,2)],[vertices(1,3),vertices(2,3)],'-b');
        plot3([vertices(1,1),vertices(3,1)],[vertices(1,2),vertices(3,2)],[vertices(1,3),vertices(3,3)],'-b');
        plot3([vertices(1,1),vertices(4,1)],[vertices(1,2),vertices(4,2)],[vertices(1,3),vertices(4,3)],'-b');
        plot3([vertices(5,1),vertices(6,1)],[vertices(5,2),vertices(6,2)],[vertices(5,3),vertices(6,3)],'-b');
        plot3([vertices(5,1),vertices(7,1)],[vertices(5,2),vertices(7,2)],[vertices(5,3),vertices(7,3)],'-b');
        plot3([vertices(5,1),vertices(8,1)],[vertices(5,2),vertices(8,2)],[vertices(5,3),vertices(8,3)],'-b');    
        plot3([vertices(4,1),vertices(6,1)],[vertices(4,2),vertices(6,2)],[vertices(4,3),vertices(6,3)],'-b');
        plot3([vertices(2,1),vertices(7,1)],[vertices(2,2),vertices(7,2)],[vertices(2,3),vertices(7,3)],'-b');
        plot3([vertices(2,1),vertices(8,1)],[vertices(2,2),vertices(8,2)],[vertices(2,3),vertices(8,3)],'-b');
        plot3([vertices(3,1),vertices(8,1)],[vertices(3,2),vertices(8,2)],[vertices(3,3),vertices(8,3)],'-b');
        plot3([vertices(3,1),vertices(6,1)],[vertices(3,2),vertices(6,2)],[vertices(3,3),vertices(6,3)],'-b');
        plot3([vertices(4,1),vertices(7,1)],[vertices(4,2),vertices(7,2)],[vertices(4,3),vertices(7,3)],'-b');
    end