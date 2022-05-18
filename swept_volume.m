%% New imp with prev Graph generation
clear all;
clc;
 

building1 = [3 13 0;3 20 0;8 13 0;8 20 0;3 13 20;3 20 20;8 13 20;8 20 20];
building2 = [15 3 0;15 8 0;17 3 0;17 8 0;15 3 11;15 8 11;17 3 11;17 8 11];
building3 = [12 13 0;12 25 0;17 13 0;17 25 0;12 13 22;12 25 22;17 13 22;17 25 22];
building4 = [20 9 0;20 14 0;23 9 0;23 14 0;20 9 17;20 14 17;23 9 17;23 14 17];

% th = input("Enter rotation of robot w.r.t. x-axis (CCW, degrees)"); %input theta

th1 = 0; th2 = 180;
th = linspace(th1,th2,ceil(th2-th1));
% rot = [cosd(th) -sin(th) 0;sind(th) cosd(th) 0;0 0 1];
Ro=[0 0 0; 0 1 0; 1 1 0; 1 0 0;0 0 2; 0 1 2; 1 1 2; 1 0 2];
s = size(Ro);

tStart = tic; 
R = [0 0 0;0 0 2];
for i = 1:size(th,2)
    rot = [cosd(th(i)) -sind(th(i)) 0;sind(th(i)) cosd(th(i)) 0;0 0 1];
    for k = 1:s(1)
        R_new(k,:) = rot*Ro(k,:)';
    end
    R_new([1,5],:) = [];
    R = [R;R_new];
end


A = -R;


[S1,D1]=minksum(building1,A);
[S2,D2]=minksum(building2,A);
[S3,D3]=minksum(building3,A);
[S4,D4]=minksum(building4,A);

b{1} = alphaShape(S1(:,1),S1(:,2),S1(:,3),Inf);
b{2} = alphaShape(S2(:,1),S2(:,2),S2(:,3),Inf);
b{3} = alphaShape(S3(:,1),S3(:,2),S3(:,3),Inf);
b{4} = alphaShape(S4(:,1),S4(:,2),S4(:,3),Inf);

figure;
plot3(0,0,0,'o');
hold on;

for i = 1:max(size(b))
    plot(b{i},'FaceColor','red','FaceAlpha',0)
end

grid(1).coord = [0,0,0;30,30,30];
grid(1).complexity = comvox_(grid(1).coord,b);
%% Complex Cell test



[check1,pos] = check_com(grid);
stop = 1;
while check1
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
        res_ = comvox_(grid(i).coord,b);
        if res_ == 1
            grid(i).complexity = 1;
            break
        end
    end
    stop = stop+1;
    [check1,pos] = check_com(grid);
end

tEnd = toc(tStart)

%% Valid Nodes generation

tic
valid_nodes = [];
for i = 1:max(size(grid))
    pts = grid(i).coord;
    check1_ = check(b,pts);
    np = 0; p = [];
    [np,p] = out_edge_check(pts);
    if check1_(1) == 0
        valid_nodes = [valid_nodes;pts(1,:)];
    end
    if check1_(2) == 0
        valid_nodes = [valid_nodes;pts(2,:)];
    end
    if np~=0
        for j = 1:size(p,1)
            if isequal(p(j,:),pts(1,:)) || isequal(p(j,:),pts(2,:))
                p(j,:) = [NaN NaN NaN];
                np = np - 1;
            end
        end
        p(any(isnan(p),2),:) = []; 
        if np==0
            continue
        end
        check_o = check(b,p);
        for j = 1:max(size(check_o))
            if check_o(j) == 0
                valid_nodes = [valid_nodes;p(j,:)];
            end
        end
    end
end

toc
%% Graph Generation
tic

nodes = [];

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
                if(obsCheck(zdiff(i),zdiff(i+1),indicator,b))
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
                if(obsCheck(ydiff(i),ydiff(i+1),indicator,b))
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
                if(obsCheck(xdiff(i),xdiff(i+1),indicator,b))
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

x_graph = T.x;
y_graph = T.y;
z_graph = T.z;

plot(G,'XData',x_graph,'YData',y_graph,'ZData',z_graph)

%% Mapping arbitrary point to grid vertex

start = [10 10 20];
goal = [1 17 13];
% goal = [25 10 25];
if check(b,start)==1 || check(b,goal)==1
    disp("Path doesn't exist");
    return
end
start_ = [];goal_=[];

for i = 1:max(size(grid))
    curr = grid(i).coord;
    xlim = [curr(1,1) curr(2,1)];
    ylim = [curr(1,2) curr(2,2)];
    zlim = [curr(1,3) curr(2,3)];
    if start(1)<=xlim(2) & start(1)>=xlim(1)
        if start(2)<=ylim(2) & start(2)>=ylim(1)
            if start(3)<=zlim(2) & start(3)>=zlim(1)
                start_ = [start_;i];
            end
        end
    end
    if goal(1)<=xlim(2) & goal(1)>=xlim(1)
        if goal(2)<=ylim(2) & goal(2)>=ylim(1)
            if goal(3)<=zlim(2) & goal(3)>=zlim(1)
                goal_ = [goal_;i];
            end
        end
    end
end

s = start_(1);
g = goal_(1);
vs = find_vertices(grid(s).coord);
vg = find_vertices(grid(g).coord);

ocvs = check(b,vs);
ocvg = check(b,vg);
start_n = [];
goal_n = [];
for i=1:max(size(ocvs))
    if ocvs(i)==0
        start_n = vs(i,:);
    end
end

for i=1:max(size(ocvg))
    if ocvg(i)==0
        goal_n = vg(i,:);
    end
end
sn = [];
gn = [];
for i = 1:max(size(nodes))
    if nodes(i).x == start_n(1)
        if nodes(i).y == start_n(2)
            if nodes(i).z == start_n(3)
                sn = i;
            end
        end
    end
    
    if nodes(i).x == goal_n(1)
        if nodes(i).y == goal_n(2)
            if nodes(i).z == goal_n(3)
                gn = i;
            end
        end
    end
end

%% Dijkstra

A = adjacency(G);
[cost path] = dijkstra(A,sn,gn);

cost;
path;
waypoints = nodes(path(1,:));

for i = 1:max(size(waypoints))
    waypoints_coord(i,:) = [waypoints(i).x waypoints(i).y waypoints(i).z];
end

plot3([start(1) start_n(1)],[start(2) start_n(2)],[start(3) start_n(3)],'-b','linewidth',5);
plot3(waypoints_coord(:,1),waypoints_coord(:,2),waypoints_coord(:,3),'-g','linewidth',5)
plot3([goal(1) goal_n(1)],[goal(2) goal_n(2)],[goal(3) goal_n(3)],'-b','linewidth',5);

%% Functions
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

function node = nodeGen(x,y,z,id)
    node.x = x;
    node.y = y;
    node.z = z;
    node.id = id;
end

function flag = obsCheck(node1, node2, indicator,b)
    flag = false;
    %check edge for z coordinates
    if(indicator == 1)        
      for i = node1.z:2:node2.z
          newNode = nodeGen(node1.x,node1.y,i,-1);
          if(check(b,[newNode.x,newNode.y,newNode.z]) ==1)
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
          if(check(b,[newNode.x,newNode.y,newNode.z]) ==1)
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
          if(check(b,[newNode.x,newNode.y,newNode.z]) ==1)
              flag = true;
              return;
          else
              flag = false;
          end
      end
    end
end

function vertices = find_vertices(coord)
    diff = coord(2,1)-coord(1,1);
    vertices = [coord(1,:);coord(1,1)+diff coord(1,2) coord(1,3);...
                coord(1,1) diff+coord(1,2) coord(1,3);coord(1,1) coord(1,2) diff+coord(1,3);...
                coord(2,:);coord(2,1)-diff coord(2,2) coord(2,3);...
                coord(2,1) coord(2,2)-diff coord(2,3);coord(2,1) coord(2,2) coord(2,3)-diff;];
end

function [np,p] = out_edge_check(pts)
    v = find_vertices(pts);
    np = 0;
    p = [];
    for i = 1:8
        tmp = v(i,:);
        if tmp(1) == 30
           if tmp(2) == 0
               np = np+1;
               p = [p;v(i,:)];
               continue;
           end
           if tmp(2) == 30
               np = np+1;
               p = [p;v(i,:)];
               continue
           end
           if tmp(3) == 0
               np = np+1;
               p = [p;v(i,:)];
               continue
           end
        end
        
        if tmp(2) == 30
            if tmp(3) == 30 || tmp(3) == 0
                np = np+1;
                p = [p;v(i,:)];
                continue
            end
        end
        
        if tmp(1) == 0 & tmp(3) == 30
            np = np+1;
            p = [p;v(i,:)];
            continue
        end
    end
end