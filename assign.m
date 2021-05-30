clear all; close all;
%% calibration part for estimating 0deg-Galaxy camera parameter
calibimgset =imageSet(fullfile('C:\Users\anstn\Desktop\대학원 2학년\대학원강의\김동환교수님 고급로봇공학\과제2\','calib'));
imageFileNames = calibimgset.ImageLocation;
[imagePoints, boardSize] = detectCheckerboardPoints(imageFileNames);
squareSizeInMM = 16;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
I = readimage(calibimgset,1);
imageSize = [size(I, 1),size(I, 2)];
params = estimateCameraParameters(imagePoints,worldPoints, ...
    'ImageSize',imageSize);
%showReprojectionErrors(params);
%figure;
%showExtrinsics(params);
intmatG=params.IntrinsicMatrix;
%% calibration part for estimating 45deg-Iphone camera parameter
calibimgset =imageSet(fullfile('C:\Users\anstn\Desktop\대학원 2학년\대학원강의\김동환교수님 고급로봇공학\과제2\','calib2'));
imageFileNames = calibimgset.ImageLocation;
[imagePoints, boardSize] = detectCheckerboardPoints(imageFileNames);
squareSizeInMM = 16;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
I = readimage(calibimgset,1);
imageSize = [size(I, 1),size(I, 2)];
params = estimateCameraParameters(imagePoints,worldPoints, ...
    'ImageSize',imageSize);
%showReprojectionErrors(params);
%figure;
%showExtrinsics(params);
intmatA=params.IntrinsicMatrix;
%% epipolar line searching process
deg00_img=imread('C:\Users\anstn\Desktop\대학원 2학년\대학원강의\김동환교수님 고급로봇공학\과제2\G000.jpg');
deg45_img=imread('C:\Users\anstn\Desktop\대학원 2학년\대학원강의\김동환교수님 고급로봇공학\과제2\A045.jpg');

fx1 = intmatG(1);
fy1 = intmatG(5);
u01 = size(deg00_img,2)/2;
v01 = size(deg00_img,1)/2;
c1Va=[1/fx1, 0, -u01/fx1; 0, 1/fy1, -v01/fy1; 0, 0, 1];
n=[1;0;0];
o=[1;cosd(0);sind(0)];
a=[0;-sind(0);cosd(0)];
p=[40;200*sind(0);200*cosd(0)]; % 카메라 평균 높이 40mm, 타겟간의 거리 약 200mm, 방향 0도에서 촬영 (원래는 35로줌 둘다)
rMc1=[eye(3),[sum(-p.*n);sum(-p.*o);sum(-p.*a)];[0,0,0],[1]];

fx2 = intmatA(1);
fy2 = intmatA(5);
u02 = size(deg45_img,2)/2;
v02 = size(deg45_img,1)/2;
c2Vb=[1/fx2, 0, -u02/fx2; 0, 1/fy2, -v02/fy2; 0, 0, 1];
n=[1;0;0];
o=[1;cosd(45);sind(45)];
a=[0;-sind(45);cosd(45)];
p=[60;200*sind(45);200*cosd(45)]; % 카메라 평균 높이 60mm, 타겟간의 거리 약 200mm, 방향 45도에서 촬영
rMc2=[[1,0,0;0,cos(45),sin(45);0,-sin(45),cos(45)],[sum(-p.*n);sum(-p.*o);sum(-p.*a)];[0,0,0,1]];

c1Mc2=inv(rMc1)*rMc2;
c1Rc2=c1Mc2(1:3,1:3);
c1Tc2=c1Mc2(1:3,4);
sT=[0,-c1Tc2(1),c1Tc2(2);c1Tc2(3),0,-c1Tc2(3);-c1Tc2(2),c1Tc2(1),0];
F=transpose(c1Va)*transpose(sT)*c1Rc2*c2Vb;

%% Select a point of which we want to find its correpondence.
figure(1);
imshow(deg00_img);
[u1,v1]=ginput(1);
u1 = ceil(u1);
v1 = ceil(v1);
Equ=[u1,v1,1]*F;

figure(2);
imshow(deg45_img);
x=[1:size(deg45_img,2)];
hold on;
plot(x,round(Equ(1)*x/Equ(2)+Equ(3)/Equ(2))); % 이미지 좌표계의 y의 방향이 x의 방향과 반대이므로 더해줘야함
hold off;
%% True correspondence searching process
searchspace=[];
% y픽셀 300~2324 내에서만 검색(하얀부분이 계속 검색되었음)
y_up_bound = 400;
y_down_bound = 800;
for x = 1:size(deg45_img,2)
    if round(Equ(1)*x/Equ(2)+Equ(3)/Equ(2)) >= 1*y_up_bound && round(Equ(1)*x/Equ(2)+Equ(3)/Equ(2)) <= size(deg45_img,1)-1*y_down_bound
        searchspace = [searchspace;[x,round(Equ(1)*x/Equ(2)+Equ(3)/Equ(2))]];
    end
end
visit = zeros(size(deg45_img,1),size(deg45_img,2));
searchrange=100;
calculrange=250;
tgtxy = [0,0];
state = 3; %1 = L1, 2 = L2, 3 = Pearson Correlation Coefficient
if state == 3
    tgt = -2;
else
    tgt=1000000;
end

for i = 1:length(searchspace)
    disp(i)
    for x=-searchrange:searchrange
        for y=-searchrange:searchrange
            % 중심의 좌표 (x,y)가 searchspace(i,1)+x, searchspace(i,2)+y인 점에서
            % calculrange 안에 들어오는 모든 점의 L2 Norm을 계산하는 코드.
            if searchspace(i,1)+x >= 1 && searchspace(i,1)+x <= size(deg45_img,2) && searchspace(i,2) + y >= 1 && searchspace(i,2) + y <= size(deg45_img,1) && visit(searchspace(i,2) + y,searchspace(i,1) + x) == 0
                up_limit = max(1,searchspace(i,2)+y-calculrange);
                down_limit = min(size(deg45_img,1),searchspace(i,2)+y+calculrange);
                left_limit = max(1,searchspace(i,1)+x-calculrange);
                right_limit = min(size(deg45_img,2),searchspace(i,1)+x+calculrange);
                tgt45 = deg45_img(up_limit:down_limit,left_limit:right_limit,:);
                tgt00 = deg00_img(v1-(searchspace(i,2)+y-up_limit):v1-(searchspace(i,2)+y-down_limit),u1-(searchspace(i,1)+x-left_limit):u1-(searchspace(i,1)+x-right_limit),:);
                if state == 3
                    cnt = Patch_PCC(tgt00,tgt45);
                    if cnt > tgt
                        tgt = cnt;
                        tgtxy = [searchspace(i,1)+x,searchspace(i,2)+y];
                        disp(tgtxy)
                        disp(tgt)
                    end
                else
                    if state == 1
                        cnt = Patch_L1_Norm(tgt00,tgt45);
                    elseif state == 2
                        cnt = Patch_L2_Norm(tgt00,tgt45);
                    end
                    if cnt < tgt
                        tgt = cnt;
                        tgtxy = [searchspace(i,1)+x,searchspace(i,2)+y];
                        disp(tgtxy)
                        disp(tgt)
                    end
                end
                visit(searchspace(i,2)+y,searchspace(i,1)+x) = 1;
            end
        end
    end
end

%% Plotting result
figure(3);
imshow([deg00_img,deg45_img]);
hold on;
s1=scatter(u1,v1);
s1.LineWidth = 2;
s1.MarkerEdgeColor = 'g';
s2=scatter(size(deg00_img,2)+tgtxy(1),tgtxy(2));
s2.LineWidth = 2;
s2.MarkerEdgeColor = 'r';
hold off;

%% Cost Function Set
%L1 Norm
function cost = Patch_L1_Norm(cmp,tgt)
cost = sum(abs(cmp-tgt),'all');
ele = size(tgt,1)*size(tgt,2);
cost = cost / ele;
end

%L2 Norm
function cost = Patch_L2_Norm(cmp,tgt)
cost = sum((cmp-tgt).^2,'all');
ele = size(tgt,1)*size(tgt,2);
cost = cost / ele;
end

%Pearson Correlation Coefficient
function cost = Patch_PCC(cmp,tgt)
cmp=im2double(cmp);
tgt=im2doulbe(tgt);
cost = [];
for i=1:size(cmp,3)
    %enum=sum((tgt(:,:,i)-mean(tgt(:,:,i),'all').*(cmp(:,:,i)-mean(cmp(:,:,i),'all'))),'all')
    %denum=
    ecost = sum((tgt(:,:,i)-mean(tgt(:,:,i),'all').*(cmp(:,:,i)-mean(cmp(:,:,i),'all'))),'all')/(sum((tgt(:,:,i)-mean(tgt(:,:,i),'all').^2),'all')*sum((cmp(:,:,i)-mean(cmp(:,:,i),'all').^2),'all'))^0.5;
    cost = [cost,ecost];
end
cost = norm(cost);
end