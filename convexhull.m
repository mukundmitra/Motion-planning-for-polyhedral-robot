% plot3(build_1(:,1),build_1(:,2),build_1(:,3),'*')
% I1_lower = find(build_1(:,3)== min(build_1(:,3)));
% build_1_lower = build_1(I1_lower,:);
% ind_1_lower = convhull(build_1_lower(:,[1,2]));
% conv_1_lower = build_1_lower(ind_1_lower,:);
% 
% 
% I1_upper = find(build_1(:,3)== max(build_1(:,3)));
% build_1_upper = build_1(I1_upper,:);
% ind_1_upper = convhull(build_1_upper(:,[1,2]));
% conv_1_upper = build_1_upper(ind_1_upper,:);
% 
% hold on
% plot3(conv_1_lower(:,1),conv_1_lower(:,2),conv_1_lower(:,3))
% plot3(conv_1_upper(:,1),conv_1_upper(:,2),conv_1_upper(:,3))
% plot3(S1(:,1),S1(:,2),S1(:,3),'o')
% axis([-3 31 -3 31 -3 31])


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


%%

%build1 = ones(3,3,max(size(k1)));
build1 = [];
is_vis1 = false(max(size(S1)),1);
for i = 1:max(size(k1))
    if is_vis1(k1(i,1))==0
        build1 = [build1;S1(k1(i,1),:)];
        is_vis1(k1(i,1))=1;
    else
        continue
    end
    if is_vis1(k1(i,2))==0
        build1 = [build1;S1(k1(i,2),:)];
        is_vis1(k1(i,2))=1;
    else
        continue
    end
    if is_vis1(k1(i,3))==0
        build1 = [build1;S1(k1(i,3),:)];
        is_vis1(k1(i,3))=1;
    else
        continue
    end
    %build1(:,:,i) = [S1(k1(i,1),:);S1(k1(i,2),:);S1(k1(i,3),:)]; 
end


%%
build2 = ones(3,3,max(size(k2)));

for i = 1:max(size(k2))
    build2(:,:,i) = [S2(k2(i,1),:);S2(k2(i,2),:);S2(k2(i,3),:)]; 
end

build3 = ones(3,3,max(size(k3)));

for i = 1:max(size(k3))
    build3(:,:,i) = [S3(k3(i,1),:);S3(k3(i,2),:);S3(k3(i,3),:)]; 
end

build4 = ones(3,3,max(size(k4)));

for i = 1:max(size(k4))
    build4(:,:,i) = [S4(k4(i,1),:);S4(k4(i,2),:);S4(k4(i,3),:)]; 
end


