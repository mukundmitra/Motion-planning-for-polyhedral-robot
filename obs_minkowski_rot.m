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

% plot3(A(:,1),A(:,2),A(:,3),'*',R_new(:,1),R_new(:,2),R_new(:,3),'s',building1(:,1),building1(:,2),building1(:,3),'o',S1(:,1),S1(:,2),S1(:,3),'d')

build_1(1,:) = [min(S1(:,1)),min(S1(:,2)),min(S1(:,3))];
build_1(2,:) = [max(S1(:,1)),min(S1(:,2)),min(S1(:,3))];
build_1(3,:) = [min(S1(:,1)),max(S1(:,2)),min(S1(:,3))];
build_1(4,:) = [max(S1(:,1)),max(S1(:,2)),min(S1(:,3))];
build_1(5,:) = [min(S1(:,1)),min(S1(:,2)),max(S1(:,3))];
build_1(6,:) = [max(S1(:,1)),min(S1(:,2)),max(S1(:,3))];
build_1(7,:) = [min(S1(:,1)),max(S1(:,2)),max(S1(:,3))];
build_1(8,:) = [max(S1(:,1)),max(S1(:,2)),max(S1(:,3))];

build_2(1,:) = [min(S2(:,1)),min(S2(:,2)),min(S2(:,3))];
build_2(2,:) = [max(S2(:,1)),min(S2(:,2)),min(S2(:,3))];
build_2(3,:) = [min(S2(:,1)),max(S2(:,2)),min(S2(:,3))];
build_2(4,:) = [max(S2(:,1)),max(S2(:,2)),min(S2(:,3))];
build_2(5,:) = [min(S2(:,1)),min(S2(:,2)),max(S2(:,3))];
build_2(6,:) = [max(S2(:,1)),min(S2(:,2)),max(S2(:,3))];
build_2(7,:) = [min(S2(:,1)),max(S2(:,2)),max(S2(:,3))];
build_2(8,:) = [max(S2(:,1)),max(S2(:,2)),max(S2(:,3))];

build_3(1,:) = [min(S3(:,1)),min(S3(:,2)),min(S3(:,3))];
build_3(2,:) = [max(S3(:,1)),min(S3(:,2)),min(S3(:,3))];
build_3(3,:) = [min(S3(:,1)),max(S3(:,2)),min(S3(:,3))];
build_3(4,:) = [max(S3(:,1)),max(S3(:,2)),min(S3(:,3))];
build_3(5,:) = [min(S3(:,1)),min(S3(:,2)),max(S3(:,3))];
build_3(6,:) = [max(S3(:,1)),min(S3(:,2)),max(S3(:,3))];
build_3(7,:) = [min(S3(:,1)),max(S3(:,2)),max(S3(:,3))];
build_3(8,:) = [max(S3(:,1)),max(S3(:,2)),max(S3(:,3))];

build_4(1,:) = [min(S4(:,1)),min(S4(:,2)),min(S4(:,3))];
build_4(2,:) = [max(S4(:,1)),min(S4(:,2)),min(S4(:,3))];
build_4(3,:) = [min(S4(:,1)),max(S4(:,2)),min(S4(:,3))];
build_4(4,:) = [max(S4(:,1)),max(S4(:,2)),min(S4(:,3))];
build_4(5,:) = [min(S4(:,1)),min(S4(:,2)),max(S4(:,3))];
build_4(6,:) = [max(S4(:,1)),min(S4(:,2)),max(S4(:,3))];
build_4(7,:) = [min(S4(:,1)),max(S4(:,2)),max(S4(:,3))];
build_4(8,:) = [max(S4(:,1)),max(S4(:,2)),max(S4(:,3))];