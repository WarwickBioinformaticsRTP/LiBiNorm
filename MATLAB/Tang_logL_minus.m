function minus = Tang_logL_minus(x, l, d, h, t1, t2)

%if (x == l) 
%    f_35 = t2/(t1+ t2) .* (exp(-t1*(l-x))+ exp(-t1*l-t2*x)) + 1/(t1 + t2).*(t1*exp(-t1*(l-x)) + t2*exp(-t1*l -t2*l+(l-x)));
%else
  %  f_35 = t2/(t1+ t2) .* (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*x));
%end

f_frag = indicator(x, h, l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
    1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d;
%f_frag = 0;

%LL = @(xx, ll) t2/(t1+ t2) .* (t1*exp(-t1*(ll-xx)) + t2*exp(-t1*ll-t2*xx)) ...
%    + indicator(xx, h, ll-h)./(t1+ t2) .* (t1.*exp(-t1*(ll-xx) - 2*t2*h- t1*h) + t2.*exp(-t1*ll-t2*(xx+h)))/d;

%norm = integral(@(xx)LL(xx, l), 0, l);
if (2*h<l)
 %   norm = (t2 - t2*exp(-l*(t1 + t2)))/(t1 + t2) + (exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/d/(t1 + t2);
 norm =  (exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;
else
  %  norm = (t2 - t2*exp(-l*(t1 + t2)))/(t1 + t2);
  norm =   (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;
end

minus = log((  f_frag)/norm);

end

