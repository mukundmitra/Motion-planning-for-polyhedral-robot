function [res] = check(obs,pts)
%CHECK Summary of this function goes here
%   Detailed explanation goes here

ns = size(pts);
np = ns(1);
res = false(np,1);
for i = 1:np
    for j = 1:max(size(obs))
        if inShape(obs{j},pts(i,:)) == 1
            res(i) = 1;
        end
    end
end

end