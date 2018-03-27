function [distance, theta] = GeneratePDRData(x,y,T)

for i=1:T-1
    if (x(i+1)-x(i)) > 0 && (y(i+1)-y(i)) > 0
        theta(i) = atan((y(i+1)-y(i))/(x(i+1)-x(i)));
    else if (x(i+1)-x(i)) < 0 && (y(i+1)-y(i)) > 0
        theta(i) = atan((y(i+1)-y(i))/(x(i+1)-x(i))) + pi;
    else if (x(i+1)-x(i)) < 0 && (y(i+1)-y(i)) < 0
        theta(i) = atan((y(i+1)-y(i))/(x(i+1)-x(i))) + pi;
    else
        theta(i) = atan((y(i+1)-y(i))/(x(i+1)-x(i))) + 2*pi;
        end
        end
    end
end
distance = sqrt(diff(x).^2+diff(y).^2);