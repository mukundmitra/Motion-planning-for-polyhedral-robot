function [res] = voxel_ambig(occ)
%VOXEL_AMBIG Summary of this function goes here
%   Detailed explanation goes here

res = 0;
ray15 = 0;
ray48 = 0;
ray37 = 0;
ray26 = 0;

check = 0;

if occ(5)==0
    if occ(1)==0
        ray15 = -1;
    else
        ray15 = 0;
    end
else
    if occ(1)==1
        ray15 = 1;
    else
        ray15 = 0;
    end
end

if ray15==0
    check = 0;
else
    check = 1;
end

if occ(4)==0
    if occ(8)==0
        ray48 = -1;
    else
        ray48 = 0;
    end
else
    if occ(8)==1
        ray48 = 1;
    else
        ray48 = 0;
    end
end

if check == 1
    if ray15 == ray48
        res = 0;
    else
        res = 1;
        return
    end
else
    if ray48==0
        check = 0;
    else
        check = 1;
    end
end

if occ(3)==0
    if occ(7)==0
        ray37 = -1;
    else
        ray37 = 0;
    end
else
    if occ(7)==1
        ray37 = 1;
    else
        ray37 = 0;
    end
end

if check == 1
    if ray37 == ray48
        res = 0;
    else
        res = 1;
        return
    end
else
    if ray37==0
        check = 0;
    else
        check = 1;
    end
end

if occ(2)==0
    if occ(6)==0
        ray26 = -1;
    else
        ray26 = 0;
    end
else
    if occ(6)==1
        ray26 = 1;
    else
        ray26 = 0;
    end
end

if check == 1
    if ray37 == ray26
        res = 0;
    else
        res = 1;
        return
    end
end

end

