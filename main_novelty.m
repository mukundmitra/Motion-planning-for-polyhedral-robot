clear all;
clc;

building1 = [3 13 0;3 20 0;8 13 0;8 20 0;3 13 20;3 20 20;8 13 20;8 20 20];
building2 = [15 3 0;15 8 0;17 3 0;17 8 0;15 3 11;15 8 11;17 3 11;17 8 11];
building3 = [12 13 0;12 25 0;17 13 0;17 25 0;12 13 22;12 25 22;17 13 22;17 25 22];
building4 = [20 9 0;20 14 0;23 9 0;23 14 0;20 9 17;20 14 17;23 9 17;23 14 17];


th = 30;
rot = [cosd(th) -sin(th) 0;sind(th) cosd(th) 0;0 0 1];
R=[0 0 0; 0 1 0; 1 1 0; 1 0 0;0 0 2; 0 1 2; 1 1 2; 1 0 2];
s = size(R);

for k = 1:s(1)
    R_new(k,:) = rot*R(k,:)';
end

A = -R_new;

[S1,D1]=minksum(building1,A);
[S2,D2]=minksum(building2,A);
[S3,D3]=minksum(building3,A);
[S4,D4]=minksum(building4,A);

[k1,vol1] = convhulln(S1);
[k2,vol2] = convhulln(S2);
[k3,vol3] = convhulln(S3);
[k4,vol4] = convhulln(S4);

figure
trisurf(k1,S1(:,1),S1(:,2),S1(:,3),'FaceColor','cyan')
hold on
trisurf(k2,S2(:,1),S2(:,2),S2(:,3),'FaceColor','cyan')
hold on
trisurf(k3,S3(:,1),S3(:,2),S3(:,3),'FaceColor','cyan')
hold on
trisurf(k4,S4(:,1),S4(:,2),S4(:,3),'FaceColor','cyan')

v11 = S1(k1(:,1),:);
v12 = S1(k1(:,2),:);
v13 = S1(k1(:,3),:);

v21 = S2(k2(:,1),:);
v22 = S2(k2(:,2),:);
v23 = S2(k2(:,3),:);

v31 = S3(k3(:,1),:);
v32 = S3(k3(:,2),:);
v33 = S3(k3(:,3),:);

v41 = S4(k4(:,1),:);
v42 = S4(k4(:,2),:);
v43 = S4(k4(:,3),:);

f1{1} = v11; f1{2} = v12; f1{3} = v13;
f2{1} = v21; f2{2} = v22; f2{3} = v23;
f3{1} = v31; f3{2} = v32; f3{3} = v33;
f4{1} = v41; f4{2} = v42; f4{3} = v43;

boundary(1).surface = k1;
boundary(2).surface = k2;
boundary(3).surface = k3;
boundary(4).surface = k4;

boundary(1).vertices = S1;
boundary(2).vertices = S2;
boundary(3).vertices = S3;
boundary(4).vertices = S4;



%% Complex Cell test
tStart = tic;
grid(1).coord = [0,0,0;30,30,30];
grid(1).complexity = comvox__(grid(1).coord,boundary);



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
        res_ = comvox__(grid(i).coord,boundary);
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
    check1_ = check_(boundary,pts);
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
        check_o = check_(boundary,p);
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
                if(obsCheck(zdiff(i),zdiff(i+1),indicator,boundary))
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
                if(obsCheck(ydiff(i),ydiff(i+1),indicator,boundary))
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
                if(obsCheck(xdiff(i),xdiff(i+1),indicator,boundary))
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
goal = [25 10 25];
if check_(boundary,start)==1 || check_(boundary,goal)==1
    disp("Path doesn't exist");
    quit
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

ocvs = check_(boundary,vs);
ocvg = check_(boundary,vg);
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
%% Novel rayIntersection Approach
function [path] = novelty(start,goal,wp,map3D)
    disp("Inside novelty")
    tot = size(wp);
    loop = tot(1);
    sp = [start(1) start(2) start(3) 1 0 0 0];
    check = loop;
    %fp = flip(wp);
    next = 1;
    for j=1:loop
        az = atan2(wp(j,2)-goal(2),wp(j,1)-goal(1));
        el = atan2(wp(j,3)-goal(3),sqrt(pdist([wp(j,:);goal],'euclidean')^2 - (wp(j,3)-goal(3))^2));
        maxrange = pdist([wp(j,:);start],'euclidean');
        directions = [az el];
        [~,is] = rayIntersection(map3D,sp,directions,maxrange);
        if is==-1
            next = j;
            break
        end
        next = j;
    end
    
    for i = 1:loop
        z = loop-i+1;
        az = atan2(wp(z,2)-start(2),wp(z,1)-start(1));
        el = atan2(wp(z,3)-start(3),sqrt(pdist([wp(z,:);start],'euclidean')^2 - (wp(z,3)-start(3))^2));
        maxrange = pdist([wp(z,:);start],'euclidean');
        directions = [az el];
        [~,is] = rayIntersection(map3D,sp,directions,maxrange);
        if is==-1
            check = z;
            break
        end
        check = z;
    end
    path = [start;wp(check:loop,:)];
end

%% Other functions
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
          if(check_(b,[newNode.x,newNode.y,newNode.z]) ==1)
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
          if(check_(b,[newNode.x,newNode.y,newNode.z]) ==1)
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
          if(check_(b,[newNode.x,newNode.y,newNode.z]) ==1)
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