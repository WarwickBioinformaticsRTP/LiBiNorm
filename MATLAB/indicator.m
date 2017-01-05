function [ output ] = indicator( x, c1, c2 )
%output = zeros(length(x));
for i = 1:length(x)
if ((x(i)>c1) && (x(i)<c2))
    output(i) = 1;
else output(i) = 0;
end
end
end

