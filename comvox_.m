function [res] = comvox_(coord,obs)
%COMVOX Summary of this function goes here
%   Detailed explanation goes here
diff = coord(2,1)-coord(1,1);
if diff<5
    res = 0;
    return;
end
vertices = [coord(1,:);coord(1,1)+diff coord(1,2) coord(1,3);...
            coord(1,1) diff+coord(1,2) coord(1,3);coord(1,1) coord(1,2) diff+coord(1,3);...
            coord(2,:);coord(2,1)-diff coord(2,2) coord(2,3);...
            coord(2,1) coord(2,2)-diff coord(2,3);coord(2,1) coord(2,2) coord(2,3)-diff;];
occ = check(obs,vertices);
res = 0; % Assuming cell is not complex

%Checking for complex Voxel
if max(size(find(occ==0)))==8
    sampled_voxel = freeVox(vertices);
    s = size(sampled_voxel);
    for i = 1:s(1)
        if check(obs,sampled_voxel(i,:))== 1
            res = 1;
            return
        end
    end
elseif max(size(find(occ==1)))==8
    sampled_voxel = freeVox(vertices);
    s = size(sampled_voxel);
    for i = 1:s(1)
        if check(obs,sampled_voxel(i,:))==0
            res = 1;
            return
        end
    end
else
    res = 0;
end
face1xy = [coord(1,:);vertices(2,:);vertices(8,:);vertices(3,:)];
occ1xy = check(obs,face1xy);
face2xy = [vertices(4,:);vertices(7,:);vertices(5,:);vertices(6,:)];
occ2xy = check(obs,face2xy);
face1yz = [vertices(1,:);vertices(3,:);vertices(6,:);vertices(4,:)];
occ1yz = check(obs,face1yz);
face2yz = [vertices(2,:);vertices(8,:);vertices(5,:);vertices(7,:)];
occ2yz = check(obs,face2yz);
face1zx = [vertices(1,:);vertices(4,:);vertices(7,:);vertices(2,:)];
occ1zx = check(obs,face1zx);
face2zx = [vertices(3,:);vertices(6,:);vertices(5,:);vertices(8,:)];
occ2zx = check(obs,face2zx);

%Checking for face ambiguity
res = face_ambig(occ1xy);
if res == 1
    return
end

res = face_ambig(occ2xy);
if res == 1
    return
end

res = face_ambig(occ1yz);
if res == 1
    return
end

res = face_ambig(occ2yz);
if res == 1
    return
end

res = face_ambig(occ1zx);
if res == 1
    return
end

res = face_ambig(occ2zx);
if res == 1
    return
end

%Checking for complexe faces
if max(size(find(check(obs,face1xy)==0)))==4 || max(size(find(check(obs,face1xy)==1)))==4
    if occ1xy(1) == 0
        sampled_face = freeFace(face1xy);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==1
                res = 1;
                return
            end
        end
    else
        sampled_face = freeFace(face1xy);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==0
                res = 1;
                return
            end
        end
    end
else
    res = 0;
end

if max(size(find(check(obs,face2xy)==0)))==4 || max(size(find(check(obs,face2xy)==1)))==4
    if occ2xy(1) == 0
        sampled_face = freeFace(face2xy);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==1
                res = 1;
                return
            end
        end
    else
        sampled_face = freeFace(face2xy);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==0
                res = 1;
                return
            end
        end
    end
else
    res = 0;
end

if max(size(find(check(obs,face1yz)==0)))==4 || max(size(find(check(obs,face1yz)==1)))==4
    if occ1yz(1) == 0
        sampled_face = freeFace(face1yz);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==1
                res = 1;
                return
            end
        end
    else
        sampled_face = freeFace(face1yz);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==0
                res = 1;
                return
            end
        end
    end
else
    res = 0;
end

if max(size(find(check(obs,face2yz)==0)))==4 || max(size(find(check(obs,face2yz)==1)))==4
    if occ2yz(1) == 0
        sampled_face = freeFace(face2yz);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==1
                res = 1;
                return
            end
        end
    else
        sampled_face = freeFace(face2yz);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==0
                res = 1;
                return
            end
        end
    end
else
    res = 0;
end

if max(size(find(check(obs,face1zx)==0)))==4 || max(size(find(check(obs,face1zx)==1)))==4
    if occ1zx(1) == 0
        sampled_face = freeFace(face1zx);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==1
                res = 1;
                return
            end
        end
    else
        sampled_face = freeFace(face1zx);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==0
                res = 1;
                return
            end
        end
    end
else
    res = 0;
end

if max(size(find(check(obs,face2zx)==0)))==4 || max(size(find(check(obs,face2zx)==1)))==4
    if occ2zx(1) == 0
        sampled_face = freeFace(face2zx);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==1
                res = 1;
                return
            end
        end
    else
        sampled_face = freeFace(face2zx);
        s = size(sampled_face);
        for i = 1:s(1)
            if check(obs,sampled_face(i,:))==0
                res = 1;
                return
            end
        end
    end
else
    res = 0;
end

%Checking for ambiguity Voxel
res = voxel_ambig(occ);

end

