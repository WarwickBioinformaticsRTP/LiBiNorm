function plus = Tang_logL_plus(x, l, d, h, t1, t2)

%if (x == 0) 
%    f_53 = t1.*exp(-t1*(l-x)) + exp(-t1*l);
%else
   % f_53 = t1.*exp(-t1*(l-x));
%end
f_frag = indicator(x, h, l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
    1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d;
%f_frag = 0;

%LL = @(xx, ll) t1.*exp(-t1*(ll-xx)) ...
%    + indicator(xx, h, ll-h)./(t1+ t2) .* (t1.*exp(-t1*(ll-xx) - 2*t2*h- t1*h) + t2.*exp(-t1*l-t2*(xx+h)))/d;

%norm = integral(@(xx)LL(xx, l), 0, l);
if (2*h<l)
   % norm = 1 -exp(-l*t1) + (exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/d/(t1 + t2);
     norm =  (exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;
else
   % norm = 1 -exp(-l*t1);
   norm =  (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;
end

plus = log(( f_frag)/norm);

end



