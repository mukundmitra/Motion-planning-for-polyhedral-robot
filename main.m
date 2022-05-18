clear all;
clc;

tStart = tic;  
map3D = occupancyMap3D;
[xGround,yGround,zGround] = meshgrid(0:30,0:30,0);
xyzGround = [xGround(:) yGround(:) zGround(:)];
occval = 0;
setOccupancy(map3D,xyzGround,occval);

inflate(map3D,0.0001);
[xBuilding1,yBuilding1,zBuilding1] = meshgrid(3:8,13:20,0:20);
[xBuilding2,yBuilding2,zBuilding2] = meshgrid(15:17,3:8,0:11);
[xBuilding3,yBuilding3,zBuilding3] = meshgrid(12:17,13:25,0:22);
[xBuilding4,yBuilding4,zBuilding4] = meshgrid(20:23,9:14,0:17);

xyzBuildings = [xBuilding1(:) yBuilding1(:) zBuilding1(:);...
                xBuilding2(:) yBuilding2(:) zBuilding2(:);...
                xBuilding3(:) yBuilding3(:) zBuilding3(:);...
                xBuilding4(:) yBuilding4(:) zBuilding4(:)];

obs = 0.65;
updateOccupancy(map3D,xyzBuildings,obs)
xlabel('X')
ylabel('Y')
zlabel('Z')
zlim([-5 35])
show(map3D)
hold on;

grid(1).coord = [0,0,0;30,30,30];
grid(1).complexity = comvox(grid(1).coord,map3D);
delpos = 0;

% for i = 1:max(size(grid))
%     coord = grid(i).coord;
%     diff = coord(2,1)-coord(1,1);
%     vertices = [coord(1,:);coord(1,1)+diff coord(1,2) coord(1,3);...
%                 coord(1,1) diff+coord(1,2) coord(1,3);coord(1,1) coord(1,2) diff+coord(1,3);...
%                 coord(2,:);coord(2,1)-diff coord(2,2) coord(2,3);...
%                 coord(2,1) coord(2,2)-diff coord(2,3);coord(2,1) coord(2,2) coord(2,3)-diff;];
%     plot3(vertices(:,1),vertices(:,2),vertices(:,3),'*r');
%     plot3([vertices(1,1),vertices(2,1)],[vertices(1,2),vertices(2,2)],[vertices(1,3),vertices(2,3)],'-b');
%     plot3([vertices(1,1),vertices(3,1)],[vertices(1,2),vertices(3,2)],[vertices(1,3),vertices(3,3)],'-b');
%     plot3([vertices(1,1),vertices(4,1)],[vertices(1,2),vertices(4,2)],[vertices(1,3),vertices(4,3)],'-b');
%     plot3([vertices(5,1),vertices(6,1)],[vertices(5,2),vertices(6,2)],[vertices(5,3),vertices(6,3)],'-b');
%     plot3([vertices(5,1),vertices(7,1)],[vertices(5,2),vertices(7,2)],[vertices(5,3),vertices(7,3)],'-b');
%     plot3([vertices(5,1),vertices(8,1)],[vertices(5,2),vertices(8,2)],[vertices(5,3),vertices(8,3)],'-b');    
%     plot3([vertices(4,1),vertices(6,1)],[vertices(4,2),vertices(6,2)],[vertices(4,3),vertices(6,3)],'-b');
%     plot3([vertices(2,1),vertices(7,1)],[vertices(2,2),vertices(7,2)],[vertices(2,3),vertices(7,3)],'-b');
%     plot3([vertices(2,1),vertices(8,1)],[vertices(2,2),vertices(8,2)],[vertices(2,3),vertices(8,3)],'-b');
%     plot3([vertices(3,1),vertices(8,1)],[vertices(3,2),vertices(8,2)],[vertices(3,3),vertices(8,3)],'-b');
%     plot3([vertices(3,1),vertices(6,1)],[vertices(3,2),vertices(6,2)],[vertices(3,3),vertices(6,3)],'-b');
%     plot3([vertices(4,1),vertices(7,1)],[vertices(4,2),vertices(7,2)],[vertices(4,3),vertices(7,3)],'-b');
% end

%% complex cell test
tic

[check,pos] = check_com(grid);
stop = 1;
while check
    if pos~=max(size(grid))
        current = 1;
        for i = pos+1:max(size(grid))
            temp(current).coord = grid(i).coord;
            temp(current).complexity = grid(i).complexity;
            current = current+1;
        end
        new_coord = decompose_grid(grid(pos).coord);
        for j = 0:7
            grid(pos+j).coord = new_coord(j+1).coord;
            grid(pos+j).complexity = 0;
        end
        for i = 1:current-1
            grid(pos+7+i).coord = temp(i).coord;
            grid(pos+7+i).complexity = temp(i).complexity;
        end
    else
        new_coord = decompose_grid(grid(pos).coord);
        for i = 0:7
            grid(pos+i).coord = new_coord(i+1).coord;
            grid(pos+i).complexity = 0;
        end
    end
    
    for i = pos:max(size(grid))
        res_ = comvox(grid(i).coord,map3D);
        if res_ == 1
            grid(i).complexity = 1;
            break
        end
    end
    stop = stop+1;
    [check,pos] = check_com(grid);
end

toc
%% valid nodes generation

tic
valid_nodes = [];
unvalid_nodes = [];
for i = 1:max(size(grid))
    for j = 1:2
        if(grid(i).coord(j,3)==0)
            grid(i).coord(j,3)= 0.1;
        end
        if getOccupancy(map3D,ceil(grid(i).coord(j,:)))<0.65
            valid_nodes = [valid_nodes;grid(i).coord(j,:)];
        else
            unvalid_nodes = [unvalid_nodes;grid(i).coord(j,:)];
        end
    end
end

toc
%% plot nodes

save('valid_nodes.mat');

% plot3(1,1,1,'o');
% hold on;
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
hold on;
% plot3(valid_nodes(:,1),valid_nodes(:,2),valid_nodes(:,3),'o');

%% graph generation

tic

nodes = [];
file = load('valid_nodes.mat');

for i= 1:length(valid_nodes)
%     newNode = nodeGen(file.valid_nodes(i,1),file.valid_nodes(i,2),file.valid_nodes(i,3),i);
    newNode = nodeGen(valid_nodes(i,1),valid_nodes(i,2),valid_nodes(i,3),i);
    nodes = [nodes newNode];
end

s = [];
t = [];
w = [];

xdiff = [];
ydiff = [];
zdiff = [];

xmin = 2;
ymin = 5;
zmin = 11;

nodesCopy = nodes;
T = struct2table(nodesCopy);

sortedZ = sortrows(T, [1 2 3]); 
zdiff = table2struct(sortedZ);

sortedY = sortrows(T, [1 3 2]); 
ydiff = table2struct(sortedY);

sortedX = sortrows(T, [2 3 1]); 
xdiff = table2struct(sortedX);

for i = 1:length(zdiff)
    if(i<length(zdiff))
        if((zdiff(i).x == zdiff(i+1).x) & (zdiff(i).y == zdiff(i+1).y))
            flag = false;
            if(zdiff(i+1).z - zdiff(i).z < zmin)
                flag = true;
            elseif(zdiff(i+1).z - zdiff(i).z >= zmin)
                indicator = 1;
                if(obsCheck(zdiff(i),zdiff(i+1),indicator,map3D))
                    flag = false;
                else 
                    flag = true;
                end
            end
            if(flag)
                s = [s zdiff(i).id];
                t = [t zdiff(i+1).id];
                weight = zdiff(i+1).z - zdiff(i).z;
                w = [w weight];
            end
        end
    end
end

for i = 1:length(ydiff)
    if(i<length(ydiff))
        if((ydiff(i).x == ydiff(i+1).x) & (ydiff(i).z == ydiff(i+1).z))
            flag = false;
            if(ydiff(i+1).y - ydiff(i).y < ymin)
                flag = true;
            elseif(ydiff(i+1).y - ydiff(i).y >= ymin)
                indicator = 2;
                if(obsCheck(ydiff(i),ydiff(i+1),indicator,map3D))
                    flag = false;
                else 
                    flag = true;
                end
            end
            if(flag)
                s = [s ydiff(i).id];
                t = [t ydiff(i+1).id];
                weight = ydiff(i+1).y - ydiff(i).y;
                w = [w weight];
            end
        end
    end
end

for i = 1:length(xdiff)
    if(i<length(xdiff))
        if(xdiff(i).y == xdiff(i+1).y) & (xdiff(i).z == xdiff(i+1).z)
            flag = false;
            if(xdiff(i+1).x - xdiff(i).x < xmin)
                flag = true;
            elseif(xdiff(i+1).x - xdiff(i).x >= xmin)
                indicator = 3;
                if(obsCheck(xdiff(i),xdiff(i+1),indicator,map3D))
                    flag = false;
                else 
                    flag = true;
                end
            end
            if(flag)
                s = [s xdiff(i).id];
                t = [t xdiff(i+1).id];
                weight = xdiff(i+1).x - xdiff(i).x;
                w = [w weight];
            end
        end
    end
end

G = graph(s,t,w);

toc 
%% dijkstra
tic

A = adjacency(G);
[cost path] = dijkstra(A,2,300);

cost;
path;

toc

%% plot path and connectivity graph

waypoints = nodes(path(1,:));

start = [waypoints(end).x,waypoints(end).y,waypoints(end).z]
goal = [waypoints(1).x,waypoints(1).y,waypoints(1).z]

for i = 1:max(size(waypoints))
    waypoints_coord(i,:) = [waypoints(i).x waypoints(i).y waypoints(i).z];
end

hold on
plot3(waypoints_coord(:,1),waypoints_coord(:,2),waypoints_coord(:,3),'-g','linewidth',5)
hold off

figure
plot(G)

tEnd = toc(tStart)
%%
function [res,pos] = check_com(grid)
    res = false;
    pos = 0;
    for i = 1:max(size(grid))
        if grid(i).complexity == 1
            res = true;
            pos = i;
            return;
        end
    end
    disp('Res is ');disp(res);disp('Pos is ');disp(pos);
end

%%
function flag = obsCheck(node1, node2, indicator,map3D)
    flag = false;
    %check edge for z coordinates
    if(indicator == 1)        
      for i = node1.z:2:node2.z
          newNode = nodeGen(node1.x,node1.y,i,-1);
          if(getOccupancy(map3D,ceil([newNode.x,newNode.y,newNode.z])) >=0.65)
              flag = true;
              return;
          else
              flag = false;
          end
      end
    %check edge for y coordinates    
    elseif(indicator == 2)
      for i = node1.y:2:node2.y
          newNode = nodeGen(node1.x, i, node1.z,-1);
          if(getOccupancy(map3D,ceil([newNode.x,newNode.y,newNode.z])) >=0.65)
              flag = true;
              return;
          else
              flag = false;
          end
      end
    %check edge for x coordinates
    elseif(indicator == 3)
      for i = node1.x:2:node2.x
          newNode = nodeGen(i,node1.y,node1.z,-1);
          if(getOccupancy(map3D,ceil([newNode.x,newNode.y,newNode.z])) >=0.65)
              flag = true;
              return;
          else
              flag = false;
          end
      end
    end
end
%%
function node = nodeGen(x,y,z,id)
    node.x = x;
    node.y = y;
    node.z = z;
    node.id = id;
end