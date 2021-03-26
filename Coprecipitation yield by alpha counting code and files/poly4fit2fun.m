function [fun] = poly4fit2fun(R)
%fits E versus depth

A(:,1)=R(:,1)./10000; %micrometer
A(:,2)=R(:,2).*1E-6;  %MeV

fittt=fit(A(:,2),A(:,1),'poly4');
coeffs1 = coeffvalues(fittt);
A = coeffs1(1,1);
B = coeffs1(1,2);
C = coeffs1(1,3);
D = coeffs1(1,4);
E = coeffs1(1,5);

%%%%%%%%%%%%%%%%fitting function above to pA vs Elec power
fun=[A B C D E];

%syms fun 
%fun = @(x) A*x.^4 + B*x.^3 + C*x.^2 + D*x+ E;

end
