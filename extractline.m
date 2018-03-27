function [linedata] = extractline(cad_in)

a = imread(cad_in);
[h,w,r] = size(a);
linedata = [];
n = 1;

for i = 1:h
    for j = 1:w
        if a(i,j,1) <= 230 || a(i,j,2) <= 230 || a(i,j,3) <= 230  %�޸��˼��ǰ�ɫ�����ص�ļ����ֵ
            linedata(1,n) = i;
            linedata(2,n) = j;
            n = n + 1;
        end
    end
end
 
end