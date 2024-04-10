%% Q 1 \
clear all 
clc
n0 = 1.545;
theta_1 = 47*pi/180;
alpha = 48*pi/180;
lambda = 1450e-9;
%a
% its need to be with the bruster angle the polarization should be Parallel polarization
theta_b = atan(n0)*(180/pi); 
%% 1 b
% Linear polarizer with 'a' angle and than hitting the mirror with
% theta1
I_0 = 5;
I_1 = 5/2;
E_1 = sqrt(I_1);
E_1_p = E_1*cos(alpha);
E_1_s = E_1*sin(alpha);

theta_i = theta_1;
theta_t = asin(sin(theta_1)/n0);
t_p = (2*sin(theta_t)*cos(theta_i))/(sin(theta_i+theta_t)*cos(theta_t-theta_i));
t_s = (2*sin(theta_t)*cos(theta_i))/(sin(theta_i+theta_t));
T_p = n0*cos(theta_t)/cos(theta_i)*t_p^2;
T_s = n0*cos(theta_t)/cos(theta_i)*t_s^2;

E_2_p = t_p*E_1_p;
E_2_s = t_s*E_1_s;
I_2 = abs(E_1_s)^2*T_s+abs(E_1_p)^2*T_p; 

%% 1 c:
syms d t 
c =3*1e8;
dlambda = 0.5e-9;
ti = asin(n0*sin(t));
rs =  sin(t-ti)/sin(t+ti);
f1 = c/lambda;
f2 = c/(lambda +dlambda);
df = abs(f1-f2);
eq1 = (1-rs^2)^2 == 0;
t = double([0,1,0]*solve(eq1))
eq2 = df == c/(2*n0*d*cos(t));
d = double((solve(eq2)))
%% d
a = [1e25,1e27,1e28];
% l_vec = (linspace(lambda-5*dlambda,lambda+5*dlambda,1e5)).';
lambda_vec_um = (lambda-5*dlambda:dlambda*1e-4:lambda+5*dlambda).';
v = c./lambda_vec_um;
n = @(v) n0 -a./v.^2;
% n = @(v) n0;
ti = real(asin(n0*sin(t)));
tt = asin(1./n(v)*sin(ti));
rs = @(tt,ti) sin(tt-ti)./sin(tt+ti);
delta = 4*pi*n(v).*d.*cos(tt)./lambda_vec_um;
T = (1 -rs(tt,ti).^2).^2./((1- rs(tt,ti).^2).^2 +4*rs(tt,ti).^2.*sin(delta/2).^2);
figure
semilogy(lambda_vec_um*1e9,T+1e-33)
legend(['a = 10^{26}';'a = 10^{27}';'a = 10^{28}'],'fontsize',12)
xlabel("$\lambda$ [nm]",'interpreter','latex','fontsize',16)
ylabel("$T$ [a.u]",'interpreter','latex','fontsize',16)
set(gca,'fontsize',12)


