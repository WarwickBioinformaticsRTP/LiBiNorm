function [LogL] = FLL_Tang(param, data)

d = 10^(param(1));
h = 10^(param(2));
t1 = 10^(param(3));
t2 = 10^(param(4));
%sig = 10^(param(5));
x = data(1, :);
l = data(2, :);
freq_l = data(3, :);

f_frag = (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d;

if (isempty(find(f_frag == 0, 1)))

    norm =  (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
        (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;

    LogL = -2*sum(log(f_frag./norm)./freq_l);
else
    LogL = realmax;

end

end

