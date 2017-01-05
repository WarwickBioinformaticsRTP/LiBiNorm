function [LogL] = FLL_ModelA(param, data)

d = 10^(param(1));
h = 10^(param(2));

%sig = 10^(param(5));
x = data(1, :);
l = data(2, :);
freq_l = data(3, :);

f_frag = (x> h).*(x < l-h) + ...
    1/d;

norm =  (2*h<l).*(l-2*h) + ...
        l/d;

LogL = -2*sum(log(f_frag./norm)./freq_l);


end

