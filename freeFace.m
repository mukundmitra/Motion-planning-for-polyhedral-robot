function [sampled_face] = freeFace(face)
sampled_face = [];
if face(1,1)==face(2,1) && face(1,1)==face(3,1) && face(1,1)==face(4,1)
    for i = face(1,2):2:face(2,2)
        for j = face(1,3):2:face(4,3)
            sampled_face = [sampled_face;face(1,1),i,j];
        end
    end
elseif face(1,2)==face(2,2) && face(1,2)==face(3,2) && face(1,2)==face(4,2)
    for i = face(1,1):2:face(4,1)
        for j = face(1,3):2:face(2,3)
            sampled_face = [sampled_face;i,face(1,2),j];
        end
    end
elseif face(1,3)==face(2,3) && face(1,3)==face(3,3) && face(1,3)==face(4,3)
    for i = face(1,1):2:face(2,1)
        for j = face(1,2):2:face(4,2)
            sampled_face = [sampled_face;i,j,face(1,3)];
        end
    end
end
end