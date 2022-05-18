%% New imp with new Graph generation
clear all;
clc;
 

building1 = [3 13 0;3 20 0;8 13 0;8 20 0;3 13 20;3 20 20;8 13 20;8 20 20];
building2 = [15 3 0;15 8 0;17 3 0;17 8 0;15 3 11;15 8 11;17 3 11;17 8 11];
building3 = [12 13 0;12 25 0;17 13 0;17 25 0;12 13 22;12 25 22;17 13 22;17 25 22];
building4 = [20 9 0;20 14 0;23 9 0;23 14 0;20 9 17;20 14 17;23 9 17;23 14 17];

% th = input("Enter rotation of robot w.r.t. x-axis (CCW, degrees)"); %input theta

th = 30;
rot = [cosd(th) -sin(th) 0;sind(th) cosd(th) 0;0 0 1];
R=[0 0 0; 0 1 0; 1 1 0; 1 0 0;0 0 2; 0 1 2; 1 1 2; 1 0 2];
s = size(R);

tStart = tic; 

for k = 1:s(1)
    R_new(k,:) = rot*R(k,:)';
end


A = -R_new;


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
    plot(b{i})
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
%% Valid Node and Graph Generation
tic
valid_nodes = [];
nodes = [];
hash = hashtable;
count = 0;
for i = 1:max(size(grid))
    pts = grid(i).coord;
    check1_ = check(b,pts);
    np = 0; p = [];
    [np,p] = out_edge_check(pts);
    if check1_(1) == 0
        count = count+1;
        valid_nodes = [valid_nodes;pts(1,:)];
        newNode = nodeGen(valid_nodes(count,1),valid_nodes(count,2),valid_nodes(count,3),count);
        nodes = [nodes newNode];
        hash = put(hash,num2str(pts(1,:)),count);
    end
    if check1_(2) == 0
        count = count+1;
        valid_nodes = [valid_nodes;pts(2,:)];
        newNode = nodeGen(valid_nodes(count,1),valid_nodes(count,2),valid_nodes(count,3),count);
        nodes = [nodes newNode];
        hash = put(hash,num2str(pts(2,:)),count);
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
                count = count+1;
                valid_nodes = [valid_nodes;p(j,:)];
                newNode = nodeGen(valid_nodes(count,1),valid_nodes(count,2),valid_nodes(count,3),count);
                nodes = [nodes newNode];
                hash = put(hash,num2str(p(j,:)),count);
            end
        end
    end
end

% Creating graph from nodes
s = []; t = []; w = [];
for i = 1:max(size(grid))
    vs = find_vertices(grid(i).coord);
    [ns,vns] = find_valid(hash,vs);
    [s,t,w,hash] = AddEdges(hash,s,t,w,grid(i).coord,ns,vns);
end
G = graph(s,t,w);
toc

T = struct2table(nodes);
x_graph = T.x;
y_graph = T.y;
z_graph = T.z;

plot(G,'XData',x_graph,'YData',y_graph,'ZData',z_graph)

%% Functions
function [so,to,wo,hash] = AddEdges(hash,s,t,w,coord,ns,vns)
    so = s;
    to = t;
    wo = w;
    diff = coord(2,1) - coord(1,1);
    if ismember(1,vns)
        if ismember(3,vns)
            if ~iskey(hash,num2str(ns(find(vns==3)))+" to "+num2str(ns(find(vns==1))))
                so = [so ns(find(vns==1))]; to = [to ns(find(vns==3))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==3)))+" to "+num2str(ns(find(vns==1))),[1]);
                hash  = put(hash,num2str(ns(find(vns==1)))+" to "+num2str(ns(find(vns==3))),[1]);
            end
        end
        if ismember(2,vns)
            if ~iskey(hash,num2str(ns(find(vns==2)))+" to "+num2str(ns(find(vns==1))))
                so = [so ns(find(vns==1))]; to = [to ns(find(vns==2))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==2)))+" to "+num2str(ns(find(vns==1))),[1]);
                hash  = put(hash,num2str(ns(find(vns==1)))+" to "+num2str(ns(find(vns==2))),[1]);
            end
        end
        if ismember(4,vns)
            if ~iskey(hash,num2str(ns(find(vns==4)))+" to "+num2str(ns(find(vns==1))))
                so = [so ns(find(vns==1))]; to = [to ns(find(vns==4))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==4)))+" to "+num2str(ns(find(vns==1))),[1]);
                hash  = put(hash,num2str(ns(find(vns==1)))+" to "+num2str(ns(find(vns==4))),[1]);
            end
        end
    end
    
    if ismember(6,vns)
        if ismember(3,vns)
            if ~iskey(hash,num2str(ns(find(vns==3)))+" to "+num2str(ns(find(vns==6))))
                so = [so ns(find(vns==6))]; to = [to ns(find(vns==3))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==3)))+" to "+num2str(ns(find(vns==6))),[1]);
                hash  = put(hash,num2str(ns(find(vns==6)))+" to "+num2str(ns(find(vns==3))),[1]);
            end
        end
        if ismember(5,vns)
            if ~iskey(hash,num2str(ns(find(vns==5)))+" to "+num2str(ns(find(vns==6))))
                so = [so ns(find(vns==6))]; to = [to ns(find(vns==5))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==5)))+" to "+num2str(ns(find(vns==6))),[1]);
                hash  = put(hash,num2str(ns(find(vns==6)))+" to "+num2str(ns(find(vns==5))),[1]);
            end
        end
        if ismember(4,vns)
            if ~iskey(hash,num2str(ns(find(vns==4)))+" to "+num2str(ns(find(vns==6))))
                so = [so ns(find(vns==6))]; to = [to ns(find(vns==4))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==4)))+" to "+num2str(ns(find(vns==6))),[1]);
                hash  = put(hash,num2str(ns(find(vns==6)))+" to "+num2str(ns(find(vns==4))),[1]);
            end
        end
    end
    
    if ismember(7,vns)
        if ismember(2,vns)
            if ~iskey(hash,num2str(ns(find(vns==2)))+" to "+num2str(ns(find(vns==7))))
                so = [so ns(find(vns==7))]; to = [to ns(find(vns==2))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==2)))+" to "+num2str(ns(find(vns==7))),[1]);
                hash  = put(hash,num2str(ns(find(vns==7)))+" to "+num2str(ns(find(vns==2))),[1]);
            end
        end
        if ismember(5,vns)
            if ~iskey(hash,num2str(ns(find(vns==5)))+" to "+num2str(ns(find(vns==7))))
                so = [so ns(find(vns==7))]; to = [to ns(find(vns==5))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==5)))+" to "+num2str(ns(find(vns==7))),[1]);
                hash  = put(hash,num2str(ns(find(vns==7)))+" to "+num2str(ns(find(vns==5))),[1]);
            end
        end
        if ismember(4,vns)
            if ~iskey(hash,num2str(ns(find(vns==4)))+" to "+num2str(ns(find(vns==7))))
                so = [so ns(find(vns==7))]; to = [to ns(find(vns==4))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==4)))+" to "+num2str(ns(find(vns==7))),[1]);
                hash  = put(hash,num2str(ns(find(vns==7)))+" to "+num2str(ns(find(vns==4))),[1]);
            end
        end
    end
    
    if ismember(8,vns)
        if ismember(2,vns)
            if ~iskey(hash,num2str(ns(find(vns==2)))+" to "+num2str(ns(find(vns==8))))
                so = [so ns(find(vns==8))]; to = [to ns(find(vns==2))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==2)))+" to "+num2str(ns(find(vns==8))),[1]);
                hash  = put(hash,num2str(ns(find(vns==8)))+" to "+num2str(ns(find(vns==2))),[1]);
            end
        end
        if ismember(5,vns)
            if ~iskey(hash,num2str(ns(find(vns==5)))+" to "+num2str(ns(find(vns==8))))
                so = [so ns(find(vns==8))]; to = [to ns(find(vns==5))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==5)))+" to "+num2str(ns(find(vns==8))),[1]);
                hash  = put(hash,num2str(ns(find(vns==8)))+" to "+num2str(ns(find(vns==5))),[1]);
            end
        end
        if ismember(3,vns)
            if ~iskey(hash,num2str(ns(find(vns==3)))+" to "+num2str(ns(find(vns==8))))
                so = [so ns(find(vns==8))]; to = [to ns(find(vns==3))]; wo = [wo diff];
                hash  = put(hash,num2str(ns(find(vns==3)))+" to "+num2str(ns(find(vns==8))),[1]);
                hash  = put(hash,num2str(ns(find(vns==8)))+" to "+num2str(ns(find(vns==3))),[1]);
            end
        end
    end
        
end

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

function [nodes,vertices] = find_valid(hash,vs)
    nodes = [];
    vertices = [];
    for i = 1:size(vs,1)
        if iskey(hash,num2str(vs(i,:)))
            nodes = [nodes;get(hash,num2str(vs(i,:)))];
            vertices = [vertices;i];
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
