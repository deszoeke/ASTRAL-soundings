function thw = theta_w(p,Temp,qv)
% theta_w(p[Pa],Temp[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008 and 
% eqn. 43 of Bolton, 1980, MWR, 108, 1046-1053.

Rd = 287.04;
Cp = 1005.7;

% Davies-Jones rational function coefficients
a0=7.101574;
a1=-20.68208;
a2= 16.11182;
a3=  2.574631;
a4= -5.205688;
b1= -3.552497;
b2=  3.781782;
b3= -0.6899655;
b4= -0.5929340;

% compute Theta_e from Bolton as approx of thw
Pi = (p/100000).^((Rd/Cp)*(1-0.28*qv));
th = Temp./Pi;
Tl = Tlcl(p,Temp,qv);
thw = th .* exp((3376./Tl - 2.54).*qv.*(1 + 0.81*qv));

% evaluates Davies-Jones 2008 rational function fit for thw [here in K]
ii=thw>=173.15;
X=thw(ii)/273.15;
X2=X.*X;
X3=X2.*X;
X4=X2.*X2;
thw(ii)=thw(ii)-exp((a0+a1*X+a2*X2+a3*X3+a4*X4)./(1+b1*X+b2*X2+b3*X3+b4*X4));

% % constants
% global Rd Cp 
% A = 2675;  %K
% C = 273.15; %K
% lam = Cp/Rd;
% ep=0.622
% es0=6.112 % hPa
% a=17.67;
b=243.5;

% % Davies-Jones 2008 method resembles Bolton
% if the<=257
%     dlnesdlnthe=a*b/((Temp-C+b).*(Temp-C+b));
%     thw=the-C-A*rs/(1+A*rs.*dlnesdlnthe));
% elseif the
%     
% elseif 377<=the & the<674
% 
% end

% doesn't work at all
% % initialize thw
% % thw=max(min(-35,Twet(th-273.15,qv,1e5,5)),85)+273.15;
% thw=th;
% 
% for iter=1:niter
%     % es(T) in hPa from es.m Buck
%     esat=6.1121*(1.0007 + 3.46e-8*1000).*exp((17.502*(thw-273.15))./(240.97 + thw-273.15));
%     %hPa                          %hPa
%     rs=622.*esat./(1000-esat); % saturation mixing ratio in hPa, Bolton eqn. 41
%     the=thw.*exp( (3.376./thw-0.00254).*rs.*(1+0.81e-3*rs) );
%     rath=thetae./the; % ratio of Bolton's thetae to guess thetae
%     thw=thw.*rath.^0.287;
% end
    