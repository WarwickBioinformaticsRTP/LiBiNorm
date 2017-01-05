function plus = Dan_logL_plus(x, l, d, h, t1, t2)

%if (x == 0) 
%    f_53 = t1.*exp(-t1*(l-x)) + exp(-t1*l);
%else
 %   f_53 = t1.*exp(-t1*(l-x));
%end
f_frag = indicator(x, h, l-h)./t1/(t1+ t2).*exp(-l*t1 - (2*h + x)*(t1+ t2)).* (exp(t1*l) - exp(t1*(x+h))).*(exp(x*(t1 +t2)) - exp(h*(t1 + t2))) + ...
    1./t1/(t1+ t2) .*exp(-l*t1 - x*(t1+ t2)).* (exp(t1*l) - exp(t1*x)).*(exp(x*(t1 +t2)) - 1)/d;
%f_frag = 0;

%LL = @(xx, ll) t1.*exp(-t1*(ll-xx)) ...
%    + indicator(xx, h, ll-h)./(t1+ t2) .* (t1.*exp(-t1*(ll-xx) - 2*t2*h- t1*h) + t2.*exp(-t1*l-t2*(xx+h)))/d;

%norm = integral(@(xx)LL(xx, l), 0, l);
if (2*h<l)
   % norm = 1 -exp(-l*t1) + (exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/d/(t1 + t2);
     norm =  (exp(-l*t1 - 2*h*t2)*(t1 + t2)^2 - exp(-l*(t1 + t2))*t1^2 + t1*t2*exp(-2*h*(t1 + t2))*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1 + t2)^2/t1^2/t2 + ...
         (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;
else
   % norm = 1 -exp(-l*t1);
   norm =  (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;
end

plus = log(( f_frag)/norm);

end



