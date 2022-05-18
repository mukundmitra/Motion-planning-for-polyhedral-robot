function [new_coord] = decompose_grid(coord)
%DECOMPOSE_GRID Summary of this function goes here
%   Detailed explanation goes here

d = coord(2,1)-coord(1,1);
v = [coord(1,:);coord(1,1)+d coord(1,2) coord(1,3);...
            coord(1,1) d+coord(1,2) coord(1,3);coord(1,1) coord(1,2) d+coord(1,3);...
            coord(2,:);coord(2,1)-d coord(2,2) coord(2,3);...
            coord(2,1) coord(2,2)-d coord(2,3);coord(2,1) coord(2,2) coord(2,3)-d;];
%Here
cell1 = [v(1,1) v(1,2) v(1,3)+(d/2); v(1,1)+(d/2) v(1,2)+(d/2) v(1,3)+d];
cell2 = [v(1,1)+(d/2) v(1,2) v(1,3)+(d/2);v(1,1)+d v(1,2)+(d/2) v(1,3)+d];
cell3 = [v(1,1)+(d/2) v(1,2)+(d/2) v(1,3)+(d/2);v(5,:)];
cell4 = [v(1,1) v(1,2)+(d/2) v(1,3)+(d/2);v(1,1)+(d/2) v(1,2)+d v(1,3)+d];
cell5 = [v(1,:);v(1,1)+(d/2) v(1,2)+(d/2) v(1,3)+(d/2)];
cell6 = [v(1,1)+(d/2) v(1,2) v(1,3);v(1,1)+d v(1,2)+(d/2) v(1,3)+(d/2)];
cell7 = [v(1,1)+(d/2) v(1,2)+(d/2) v(1,3);v(1,1)+d v(1,2)+d v(1,3)+(d/2)];
cell8 = [v(1,1) v(1,2)+(d/2) v(1,3);v(1,1)+(d/2) v(1,2)+d v(1,3)+(d/2)];
new_coord(1).coord = cell1;
new_coord(2).coord = cell2;
new_coord(3).coord = cell3;
new_coord(4).coord = cell4;
new_coord(5).coord = cell5;
new_coord(6).coord = cell6;
new_coord(7).coord = cell7;
new_coord(8).coord = cell8;

end

