function [res] = face_ambig(occ)
%FACE_AMBIG Summary of this function goes here
%   Detailed explanation goes here

res = 0;
ray13 = 0;
ray24 = 0;

if occ(1)==0
    if occ(3)==0
        ray13 = -1;
    else
        ray13 = 0;
        return;
    end
else
    if occ(3)==1
        ray13 = 1;
    else
        ray13 = 0;
        return;
    end
end

if occ(2)==0
    if occ(4)==0
        ray24 = -1;
    else
        ray24 = 0;
        return;
    end
else
    if occ(4)==1
        ray24 = 1;
    else
        ray24 = 0;
        return;
    end
end

if ray13==ray24
    res = 0;
else
    res = 1;
end

end

