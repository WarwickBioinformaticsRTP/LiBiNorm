function [LogL] = FLL_Deng(param, data)

d = 10^(param(1));
h = 10^(param(2));
%t1 = 10^(param(3));
t2 = 10^(param(3));
%sig = 10^(param(5));
x = data(1, :);
l = data(2, :);
freq_l = data(3, :);


    f_frag = (x> h).*(x < l-h).*exp(-t2*(x+h)) + ...
    exp(-t2*(x))/d;

if (isempty(find(f_frag == 0, 1)))

    norm =  (2*h<l).*(exp(-2*h*t2) - exp(-l*( t2)))/t2 + ...
        (1-exp(-l*( t2)))/t2/d;

    LogL = -2*sum(log(f_frag./norm)./freq_l);

else
LogL = realmax;

end
end

