clear;
clc;
close all;

%%
%轨迹、地图、墙壁信息整合部分


%读取运动轨迹并平滑

load('trajBefore.mat');
smooth = 200;
x=trajBefore;
lenX = floor(length(x)/smooth);

for i = 1:lenX
    x_mean(:,i) = mean(x(:,(i-1)*smooth+1:i*smooth),2);
end

x0 = x_mean(2,8:end - 8);
y0 = x_mean(1,8:end - 8);


%轨迹与旋转
%按指定距离距离平移到最佳位置


RMtheta = 30*pi/180;   %本案例中初始旋转角度为30
x = x0.*cos(RMtheta) - y0.*sin(RMtheta);
y = x0.*sin(RMtheta) + y0.*cos(RMtheta);



%给定PDR初始点位置
%按指定距离距离平移到最佳位置
%本案例为向右平移8.09，向下平移-27.22

x = x + 8.09;   %本案例中初始点坐标为（8.09，-27.22），给定坐标系为
y = y - 27.22;


X(1,:) = x;
X(2,:) = y;

%PDR有关数据
T = length(X);
[distance,theta] = GeneratePDRData(X(1,:),X(2,:),T);

%转换为轨迹向量矩阵

M1=zeros(length(X)-1,4);
for i=1:length(X)-1
    M1(i,1)=X(1,i);
    M1(i,2)=X(2,i);
    M1(i,3)=X(1,i+1);
    M1(i,4)=X(2,i+1);    
end

%读取地图
 
% map = 'B7.bmp';    %输入此图片，有一定概率会发生粒子数目骤减为0的情况，因此需要重新滤波
map = 'F8配色1.jpg';   %输入此图片，效果好，粒子数目不会骤减为0，且能检测到门的存在，允许粒子进入房间

rgb=imread(map);
[t1,t2,t3]=size(rgb);
wall = extractline(map);%全幅显示地图

%调用LSD程序lsd,读取地图，获取墙壁等特征向量
%LSD程序输入：地图图片
%LSD程序输出：N*5矩阵，N代表图片中直线线段的个数，（i,1：2）第i个直线线段的端点1的X、Y坐标，（i,3：4）直线线段i的端点2的X、Y坐标，（i,5）直线线段i的线宽

wall_vectors = lsd(map);
M2=zeros(length(wall_vectors),4);
wall_vectors=wall_vectors';

%转换为墙壁向量矩阵

M2(:,1)=wall_vectors(:,1);
M2(:,2)=wall_vectors(:,3);
M2(:,3)=wall_vectors(:,2);
M2(:,4)=wall_vectors(:,4);

%计算真实轨迹与地图之间的比例尺，并将地图与墙壁向量进行比例转换
% 本案例中进行了若干手动调整的比例参数（如4、6、20），建议根据实际地图手动调整，待改进

MaxX=max(X');
MinY=min(X');

scaleX=(t1/round(MaxX(1))+6);
scaleY=-1*(t2/(floor(abs(MinY(2)))+20));
wall_new(2,:)=-wall(1,:)/(t2/(floor(abs(MinY(2)))+20));
wall_new(1,:)=wall(2,:)/(t1/round(MaxX(1))+6)+4;

Xmin=min(wall_new(1,:));  %地图界限
Xmax=max(wall_new(1,:));   
Ymin=min(wall_new(2,:));
Ymax=max(wall_new(2,:));

M2(:,1)=M2(:,1)/scaleX+4;
M2(:,3)=M2(:,3)/scaleX+4;
M2(:,2)=M2(:,2)/scaleY;
M2(:,4)=M2(:,4)/scaleY;

lines=wall_vectors';
lines(1,:)=lines(1,:)/scaleX+4;
lines(2,:)=lines(2,:)/scaleX+4;
lines(3,:)=lines(3,:)/scaleY;
lines(4,:)=lines(4,:)/scaleY;

%计算轨迹向量与墙壁向量的交叉
%调用LSI（直线交叉）程序lineSegmentIntersec,计算轨迹向量组与墙壁向量组的交点
%LSI程序输入：N1*4矩阵，N2*4矩阵，矩阵的每行代表直线线段，矩阵的4列分别表示线段的端点坐标[x1 y1 x2 y2]
%LSI程序输出：结构体，包含交叉检验矩阵intAdjacencyMatrix（N1*N2）、交叉点坐标矩阵intMatrixXY、交叉点距离intNormalizedDistance1To2、平行线段parAdjacencyMatrix、重合线段coincAdjacencyMatrix等

out=lineSegmentIntersect(M1,M2);
intersectPoint=out.intMatrixXY;

%%
%画图
figure(1)

%原地图按比例缩放后画出，MarkerSize（黑点）尽量小，画出的地图效果梗好
plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;

hold on;
h1=plot(X(1,:),X(2,:),'b');axis equal;
hold on;
h2=scatter(intersectPoint(1,:),intersectPoint(2,:),'filled','d','MarkerEdgeColor','r','MarkerFaceColor','r');
title('室内行走轨迹及穿墙点检测');
legend([h1,h2],'行走轨迹','穿墙点');
xlabel('X/m');ylabel('Y/m');

%%
%while型粒子滤波（一旦粒子数目为0.将一切重置，重新粒子滤波）

N = 1000;   %粒子总数
P_r = 0.1;     %%%粒子半径
Q = 0.1;      %过程噪声
w = zeros(N,1).*1/N;       %每个粒子的权重,此处全设为0，满足进入while型循环条件

while(sum(w)==0)
    
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, T);  %所有粒子的中心位置
PCenter(:,1) = X(:,1);
w = ones(N,1).*1/N;       %每个粒子的权重

    %粒子群初始化
    for i  =  1:N
        P(:,i) = X(:,1) + wgn(2, 1, 10*log10(P_r));
    end
    
    figure(2)
    clf;
    plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;
    hold on;
    h1=plot(P(1,:),P(2,:),'r.');
    title('粒子滤波');
    legend([h1],'粒子');
    xlabel('X/m');ylabel('Y/m');
    WallPerc(1) = 1;  %粒子重采样指标
    
    %粒子群随步数更新
    for k = 1 : T-1
        
        %粒子更新
        for i = 1 : N
            Xp1 = P(:,i);
            P(:, i) = P(:,i) + distance(k) * [(cos(theta(k))); sin(theta(k))] + wgn(2, 1, 10*log10(Q));
            Xp2 = P(:,i);
            M=zeros(1,4);
            M(1,1:2)=Xp1';
            M(1,3:4)=Xp2';
            checkM=lineSegmentIntersect(M,M2);
            if Xp2(1) <= Xmin||Xp2(1)>=Xmax||Xp2(2)<=Ymin||Xp2(2)>=Ymax    %设置粒子不能越出地图界限
                w(i) = 0;
            else if checkM.intAdjacencyMatrix == 0     %检测是否穿墙
                    w(i) = 1;
                else
                    w(i) = 0;
                end
            end
        end
       
        WallPerc(k+1) = sum(w)/N;
               
        %地图信息修正后的粒子云
        Pf1(1,:) = P(1,:).*w';
        Pf1(2,:) = P(2,:).*w';
              
        %实时将新粒子群画出
%         clf;
%         plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;
        hold on;
        h1=plot(Pf1(1,:),Pf1(2,:),'r.');
        title('粒子滤波');
        legend([h1],'粒子');
        xlabel('X/m');ylabel('Y/m');
        pause(0.05);
        
        P(1,:) = P(1,:).*w';
        P(2,:) = P(2,:).*w';
        
        w_average = w./sum(w);
        Pf2(1,:) = P(1,:).*w_average';
        Pf2(2,:) = P(2,:).*w_average';
        PCenter(:, k+1) = sum(Pf2, 2);      %定位结果更新
        
         %检测到粒子数目为0，跳出粒子更新迭代，重新滤波
        if sum(w)==0
            break;
        end
        
        %重采样方法1，随机化方法
%         if WallPerc(k+1) <= 0.5
%             for h = 1 : N
%                 wmax = 2 *max(w)*rand;  
%                 index = randi(N, 1);
%                 while(wmax > w(index))
%                     wmax = wmax - w(index);
%                     index = index + 1;
%                     if index > N
%                         index = 1;
%                     end
%                 end
%                 w(h)=w(index);
%                 P(:, h) = P(:, index);     %得到新粒子
%             end
%             WallPerc(k+1) = sum(w)/N;
%         end
        
        %重采样方法2，将粒子按权重排序，将权重为1的粒子取代权重为0的粒子
        
        if WallPerc(k+1) <= 0.5
                [w_order,order]=sort(w);
                [w_order2,order2]=sort(w,'descend');
                number=N-sum(w);  %权重为0的粒子数目
                w2=w;
                flag=round(number/2);   %将权重为1的粒子取代半数权重为0的粒子，保证重采样过后有N/2个粒子
                w2(order(1:flag),1)=1;
                if(sum(w)>=flag)
                    P(:,order(1:flag))=P(:,order2(1:flag));
                else
                    if rem(flag,sum(w))~=0
                        times=ceil(flag/sum(w));
                        for h=1:times-1
                            P(:,order(1+(h-1)*sum(w):h*sum(w)))=P(:,order2(1:sum(w)));
                        end
                        P(:,order(1+h*sum(w):flag))=P(:,order2(1:rem(flag,sum(w))));
                    else
                        times=flag/sum(w);
                        for h=1:times
                            P(:,order(1+(h-1)*sum(w):h*sum(w)))=P(:,order2(1:sum(w)));
                        end
                    end
                end
                w=w2;
                WallPerc(k+1) = sum(w)/N;
        end
    end   
end

%%
%粒子滤波后的结果，设置第一个粒子中心为起始点
PCenter(:,1) = X(:,1);

%%
%画图

figure(3);
h1=plot(X(1,:),X(2,:),'b');
hold on;
plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;
hold on;

%将墙壁向量在图中画出
% for i = 1:size(lines, 2)
%      plot(lines(1:2, i), lines(3:4, i), 'LineWidth', lines(5, i) / 2, 'Color', [1, 0, 0]);
%      hold on;
%      axis equal;
% end

% scatter(intersectPoint(1,:),intersectPoint(2,:),'filled','d','MarkerEdgeColor','b','MarkerFaceColor','b');
h2=plot(PCenter(1,:),PCenter(2,:),'g');
title('粒子滤波前后定位轨迹对比');
legend([h1,h2],'粒子滤波前','粒子滤波后');
xlabel('X/m');ylabel('Y/m');


% 
% %按地图原像素画图
% figure(4);
% nCordX=(X(1,:)-4)*scaleX;
% nCordY=X(2,:)*scaleY;
% 
% nPCenter(1,:)=(PCenter(1,:)-4)*scaleX;
% nPCenter(2,:)=(PCenter(2,:))*scaleY;
% 
% imshow(rgb);
% hold on;
% h2=plot(nPCenter(1,:),nPCenter(2,:),'g')
% hold on;
% h1=plot(nCordX,nCordY,'b');
% title('粒子滤波前后定位轨迹对比');
% legend([h1,h2],'粒子滤波前','粒子滤波后');
% xlabel('X/m');ylabel('Y/m');





