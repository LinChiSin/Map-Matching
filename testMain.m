clear;
clc;
close all;

%%
%�켣����ͼ��ǽ����Ϣ���ϲ���


%��ȡ�˶��켣��ƽ��

load('trajBefore.mat');
smooth = 200;
x=trajBefore;
lenX = floor(length(x)/smooth);

for i = 1:lenX
    x_mean(:,i) = mean(x(:,(i-1)*smooth+1:i*smooth),2);
end

x0 = x_mean(2,8:end - 8);
y0 = x_mean(1,8:end - 8);


%�켣����ת
%��ָ���������ƽ�Ƶ����λ��


RMtheta = 30*pi/180;   %�������г�ʼ��ת�Ƕ�Ϊ30
x = x0.*cos(RMtheta) - y0.*sin(RMtheta);
y = x0.*sin(RMtheta) + y0.*cos(RMtheta);



%����PDR��ʼ��λ��
%��ָ���������ƽ�Ƶ����λ��
%������Ϊ����ƽ��8.09������ƽ��-27.22

x = x + 8.09;   %�������г�ʼ������Ϊ��8.09��-27.22������������ϵΪ
y = y - 27.22;


X(1,:) = x;
X(2,:) = y;

%PDR�й�����
T = length(X);
[distance,theta] = GeneratePDRData(X(1,:),X(2,:),T);

%ת��Ϊ�켣��������

M1=zeros(length(X)-1,4);
for i=1:length(X)-1
    M1(i,1)=X(1,i);
    M1(i,2)=X(2,i);
    M1(i,3)=X(1,i+1);
    M1(i,4)=X(2,i+1);    
end

%��ȡ��ͼ
 
% map = 'B7.bmp';    %�����ͼƬ����һ�����ʻᷢ��������Ŀ���Ϊ0������������Ҫ�����˲�
map = 'F8��ɫ1.jpg';   %�����ͼƬ��Ч���ã�������Ŀ�������Ϊ0�����ܼ�⵽�ŵĴ��ڣ��������ӽ��뷿��

rgb=imread(map);
[t1,t2,t3]=size(rgb);
wall = extractline(map);%ȫ����ʾ��ͼ

%����LSD����lsd,��ȡ��ͼ����ȡǽ�ڵ���������
%LSD�������룺��ͼͼƬ
%LSD���������N*5����N����ͼƬ��ֱ���߶εĸ�������i,1��2����i��ֱ���߶εĶ˵�1��X��Y���꣬��i,3��4��ֱ���߶�i�Ķ˵�2��X��Y���꣬��i,5��ֱ���߶�i���߿�

wall_vectors = lsd(map);
M2=zeros(length(wall_vectors),4);
wall_vectors=wall_vectors';

%ת��Ϊǽ����������

M2(:,1)=wall_vectors(:,1);
M2(:,2)=wall_vectors(:,3);
M2(:,3)=wall_vectors(:,2);
M2(:,4)=wall_vectors(:,4);

%������ʵ�켣���ͼ֮��ı����ߣ�������ͼ��ǽ���������б���ת��
% �������н����������ֶ������ı�����������4��6��20�����������ʵ�ʵ�ͼ�ֶ����������Ľ�

MaxX=max(X');
MinY=min(X');

scaleX=(t1/round(MaxX(1))+6);
scaleY=-1*(t2/(floor(abs(MinY(2)))+20));
wall_new(2,:)=-wall(1,:)/(t2/(floor(abs(MinY(2)))+20));
wall_new(1,:)=wall(2,:)/(t1/round(MaxX(1))+6)+4;

Xmin=min(wall_new(1,:));  %��ͼ����
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

%����켣������ǽ�������Ľ���
%����LSI��ֱ�߽��棩����lineSegmentIntersec,����켣��������ǽ��������Ľ���
%LSI�������룺N1*4����N2*4���󣬾����ÿ�д���ֱ���߶Σ������4�зֱ��ʾ�߶εĶ˵�����[x1 y1 x2 y2]
%LSI����������ṹ�壬��������������intAdjacencyMatrix��N1*N2����������������intMatrixXY����������intNormalizedDistance1To2��ƽ���߶�parAdjacencyMatrix���غ��߶�coincAdjacencyMatrix��

out=lineSegmentIntersect(M1,M2);
intersectPoint=out.intMatrixXY;

%%
%��ͼ
figure(1)

%ԭ��ͼ���������ź󻭳���MarkerSize���ڵ㣩����С�������ĵ�ͼЧ������
plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;

hold on;
h1=plot(X(1,:),X(2,:),'b');axis equal;
hold on;
h2=scatter(intersectPoint(1,:),intersectPoint(2,:),'filled','d','MarkerEdgeColor','r','MarkerFaceColor','r');
title('�������߹켣����ǽ����');
legend([h1,h2],'���߹켣','��ǽ��');
xlabel('X/m');ylabel('Y/m');

%%
%while�������˲���һ��������ĿΪ0.��һ�����ã����������˲���

N = 1000;   %��������
P_r = 0.1;     %%%���Ӱ뾶
Q = 0.1;      %��������
w = zeros(N,1).*1/N;       %ÿ�����ӵ�Ȩ��,�˴�ȫ��Ϊ0���������while��ѭ������

while(sum(w)==0)
    
P = zeros(2, N);    %��������Ⱥ
PCenter = zeros(2, T);  %�������ӵ�����λ��
PCenter(:,1) = X(:,1);
w = ones(N,1).*1/N;       %ÿ�����ӵ�Ȩ��

    %����Ⱥ��ʼ��
    for i  =  1:N
        P(:,i) = X(:,1) + wgn(2, 1, 10*log10(P_r));
    end
    
    figure(2)
    clf;
    plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;
    hold on;
    h1=plot(P(1,:),P(2,:),'r.');
    title('�����˲�');
    legend([h1],'����');
    xlabel('X/m');ylabel('Y/m');
    WallPerc(1) = 1;  %�����ز���ָ��
    
    %����Ⱥ�沽������
    for k = 1 : T-1
        
        %���Ӹ���
        for i = 1 : N
            Xp1 = P(:,i);
            P(:, i) = P(:,i) + distance(k) * [(cos(theta(k))); sin(theta(k))] + wgn(2, 1, 10*log10(Q));
            Xp2 = P(:,i);
            M=zeros(1,4);
            M(1,1:2)=Xp1';
            M(1,3:4)=Xp2';
            checkM=lineSegmentIntersect(M,M2);
            if Xp2(1) <= Xmin||Xp2(1)>=Xmax||Xp2(2)<=Ymin||Xp2(2)>=Ymax    %�������Ӳ���Խ����ͼ����
                w(i) = 0;
            else if checkM.intAdjacencyMatrix == 0     %����Ƿ�ǽ
                    w(i) = 1;
                else
                    w(i) = 0;
                end
            end
        end
       
        WallPerc(k+1) = sum(w)/N;
               
        %��ͼ��Ϣ�������������
        Pf1(1,:) = P(1,:).*w';
        Pf1(2,:) = P(2,:).*w';
              
        %ʵʱ��������Ⱥ����
%         clf;
%         plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;
        hold on;
        h1=plot(Pf1(1,:),Pf1(2,:),'r.');
        title('�����˲�');
        legend([h1],'����');
        xlabel('X/m');ylabel('Y/m');
        pause(0.05);
        
        P(1,:) = P(1,:).*w';
        P(2,:) = P(2,:).*w';
        
        w_average = w./sum(w);
        Pf2(1,:) = P(1,:).*w_average';
        Pf2(2,:) = P(2,:).*w_average';
        PCenter(:, k+1) = sum(Pf2, 2);      %��λ�������
        
         %��⵽������ĿΪ0���������Ӹ��µ����������˲�
        if sum(w)==0
            break;
        end
        
        %�ز�������1�����������
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
%                 P(:, h) = P(:, index);     %�õ�������
%             end
%             WallPerc(k+1) = sum(w)/N;
%         end
        
        %�ز�������2�������Ӱ�Ȩ�����򣬽�Ȩ��Ϊ1������ȡ��Ȩ��Ϊ0������
        
        if WallPerc(k+1) <= 0.5
                [w_order,order]=sort(w);
                [w_order2,order2]=sort(w,'descend');
                number=N-sum(w);  %Ȩ��Ϊ0��������Ŀ
                w2=w;
                flag=round(number/2);   %��Ȩ��Ϊ1������ȡ������Ȩ��Ϊ0�����ӣ���֤�ز���������N/2������
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
%�����˲���Ľ�������õ�һ����������Ϊ��ʼ��
PCenter(:,1) = X(:,1);

%%
%��ͼ

figure(3);
h1=plot(X(1,:),X(2,:),'b');
hold on;
plot(wall_new(1,:),wall_new(2,:),'k.','MarkerSize',0.01);axis equal;
hold on;

%��ǽ��������ͼ�л���
% for i = 1:size(lines, 2)
%      plot(lines(1:2, i), lines(3:4, i), 'LineWidth', lines(5, i) / 2, 'Color', [1, 0, 0]);
%      hold on;
%      axis equal;
% end

% scatter(intersectPoint(1,:),intersectPoint(2,:),'filled','d','MarkerEdgeColor','b','MarkerFaceColor','b');
h2=plot(PCenter(1,:),PCenter(2,:),'g');
title('�����˲�ǰ��λ�켣�Ա�');
legend([h1,h2],'�����˲�ǰ','�����˲���');
xlabel('X/m');ylabel('Y/m');


% 
% %����ͼԭ���ػ�ͼ
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
% title('�����˲�ǰ��λ�켣�Ա�');
% legend([h1,h2],'�����˲�ǰ','�����˲���');
% xlabel('X/m');ylabel('Y/m');





