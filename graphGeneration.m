
nodes = [];
file = load('valid_nodes.mat');

for i= 1:length(file.valid_nodes)
    newNode = nodeGen(file.valid_nodes(i,1),file.valid_nodes(i,2),file.valid_nodes(i,3),i);
    nodes = [nodes newNode];
end

s = [];
t = [];
weights = [];
% names = 

for i = 1:length(nodes)
% for i = 1:100
%     if(i<5507)
        neighbors = [];
%         if(nodes(i).x == nodes(i+1).x) & (nodes(i).y == nodes(i+1).y)
%             s = [s nodes(i).id];
%             t = [t nodes(i+1).id]; 
%             weights = [weights 1];
%         end
        
        for j = 1:length(nodes)
%         for j = 1:100
            if(nodes(i).x - nodes(j).x == 2) & (nodes(i).y == nodes(j).y) & (nodes(i).z == nodes(j).z)
                neighbors = [neighbors nodes(j)];
            elseif (nodes(i).y - nodes(j).y == 2) & (nodes(i).x == nodes(j).x) & (nodes(i).z == nodes(j).z)
                neighbors = [neighbors nodes(j)];
            elseif (nodes(i).z - nodes(j).z == 2) & (nodes(i).x == nodes(j).x) & (nodes(i).y == nodes(j).y)
                neighbors = [neighbors nodes(j)];
            end
        end

        for k = 1:length(neighbors)
            s = [s nodes(i).id];
            t = [t neighbors(k).id];
            weights = [weights 1];
        end 
end

G = graph(s,t,weights);
plot(G);

function node = nodeGen(x,y,z,id)
    node.x = x;
    node.y = y;
    node.z = z;
    node.id = id;
end

% s = [1 1 1 2 2 3 3 4 5 5 6 7];
% t = [2 4 8 3 7 4 6 5 6 8 7 8];
% weights = [10 10 1 10 1 10 1 1 12 12 12 12];
% names = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
% G = graph(s,t,weights,names);
% plot(G,'EdgeLabel',G.Edges.Weight);