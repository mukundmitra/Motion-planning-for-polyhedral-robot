function [sampled_voxel] = freeVox(vertices)
sampled_voxel = [];
for i = vertices(1,1):2:vertices(2,1)
    for j = vertices(1,2):2:vertices(3,2)
        for k = vertices(1,3):2:vertices(4,3)
            sampled_voxel = [sampled_voxel;i,j,k];
        end
    end
end
end