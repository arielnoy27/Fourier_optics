
%% Q 2
% a
% we will use Lens formula
clc

n_2 = 1.545;
n_1 = 1;
n_3 = 1;
R_1 = 5*10^(-3);
R_2 = 6*10^(-3);
f = ((n_2 - n_1)/R_1 + (n_3 - n_2)/R_2)^(-1); % f = 55.04 mm
% we can get a real image with this lens as long as u > f the object is
% placed further from the focal lens distance

%% b
lambda_arra = 0.4:0.01:1;
lambda_array = lambda_arra*1e-6;

n_BK7 = 1 + (1.03961212 * lambda_array.^2) ./ (lambda_array.^2 - 0.00600069867) + ...
    0.231792344 .* lambda_array ./ (lambda_array.^2 - 0.0200179144) + ...
    (1.01046945 .* lambda_array) ./ (lambda_array.^2 - 103.560653);

f = zeros(size(lambda_array)); 

for i = 1:length(lambda_array)
    n_BK7_i = 1 + (1.03961212 * lambda_array(i)^2) / (lambda_array(i)^2 - 0.00600069867) + ...
        0.231792344 * lambda_array(i) / (lambda_array(i)^2 - 0.0200179144) + ...
        (1.01046945 * lambda_array(i)) / (lambda_array(i)^2 - 103.560653);
    
    f(i) = ((n_BK7_i - n_1) / R_1 + (n_3 - n_BK7_i) / R_2)^(-1);
end

% Plot n_BK7
subplot(2,1,1);
plot(lambda_array, n_BK7);
xlabel('Wavelength (\mum)');
ylabel('n_{BK7}');
title('Index of Refraction (n_{BK7}) vs. Wavelength');

% Plot f
subplot(2,1,2);
plot(lambda_array, f);
xlabel('Wavelength (\mum)');
ylabel('Focal Length (m)');
title('Focal Length vs. Wavelength');

%% c 1+2
R = 25;
d1 = 70;
syms d2;
lambda_1 = 595*10^(-3);
lambda_2 = 487*10^(-3);
lambda_array = lambda_1;

n_BK7 = sqrt(1 + (1.03961212 * lambda_array^2)/ (lambda_array^2 - 0.00600069867) + ...
    (0.231792344 .* lambda_array)/ (lambda_array^2 - 0.0200179144) + ...
    (1.01046945 * lambda_array) / (lambda_array^2 - 103.560653));




M = calculateM(d1,d2,n_BK7,R);

% Define the equation
B = M(1,2);
equation = B == 0;

% Solve the equation
solution = solve(equation, d2);

% Display the solution
disp(solution);

% 2 c 2
r = 1;
theta = -pi/5:pi/25:pi/5;
x_ = -1:0.1:1;
% Initialize arrays to store beam positions after each matrix multiplication
positions_after_A = [];
positions_after_B = [];
positions_after_C = [];
positions_after_D = [];
lambda_array =lambda_1;
for i = 1:length(theta)
    bean_in = [r;theta(i)];
    bean_A = matrixA(70)*bean_in;
    positions_after_A = [positions_after_A, bean_A];
    bean_B = matrixB(R,n_BK7)*bean_A;
    positions_after_B = [positions_after_B, bean_B];
    bean_C = matrixC(R,n_BK7)*bean_B;
    positions_after_C = [positions_after_C, bean_C];
    bean_D = matrixD(32)*bean_C;
    positions_after_D = [positions_after_D, bean_D];
     % Define the x and y coordinates of the four points
     x = [0, 70, 71, 102];
     y = [bean_in(1), bean_A(1), bean_B(1), bean_D(1)];
     
     % Plot the points
     figure(1)
     plot(x, y, '--', 'MarkerSize', 100); 
     xlabel('X axis mm');
     ylabel('Y axis mm');
     title('Graph of beam path');
     hold on
end
hold off

for j = 1:2
    if j == 1
        lambda_array = lambda_1;
        n_BK7 = sqrt(1 + (1.03961212 * lambda_array^2)/ (lambda_array^2 - 0.00600069867) + ...
    (0.231792344 .* lambda_array)/ (lambda_array^2 - 0.0200179144) + ...
    (1.01046945 * lambda_array) / (lambda_array^2 - 103.560653));
    end
    if j == 2
        lambda_array = lambda_2;
        n_BK7 = sqrt(1 + (1.03961212 * lambda_array^2)/ (lambda_array^2 - 0.00600069867) + ...
    (0.231792344 .* lambda_array)/ (lambda_array^2 - 0.0200179144) + ...
    (1.01046945 * lambda_array) / (lambda_array^2 - 103.560653));
    end

    for i = 1:length(x_)
        bean_in = [x_(i);0];
        bean_A = matrixA(70)*bean_in;
        positions_after_A = [positions_after_A, bean_A];
        bean_B = matrixB(R,n_BK7)*bean_A;
        positions_after_B = [positions_after_B, bean_B];
        bean_C = matrixC(R,n_BK7)*bean_B;
        positions_after_C = [positions_after_C, bean_C];
        bean_D = matrixD(32)*bean_C;
        positions_after_D = [positions_after_D, bean_D];
         % Define the x and y coordinates of the four points
         x = [0, 70, 71, 102];
         y = [bean_in(1), bean_A(1), bean_B(1), bean_D(1)];
         
         if j == 1
         % Plot the points
         figure(2)
         plot(x, y, '--', 'MarkerSize', 100,'Color','b'); 
         hold on
         xlabel('X axis');
         ylabel('Y axis');
         title('Graph of beam path');
         
         else
             % Plot the points
             figure(2)
             plot(x, y, '--', 'MarkerSize', 100,'Color','r'); 
             legend('lambda_1','Location','south')
             legend('\color{blue}lambda_1','\color{red}lambda_2')
             xlabel('X axis');
             ylabel('Y axis');
             title('Graph of beam path');
         end
         
       
    end
end
hold off
%% 2 d 1

n_BK7 = sqrt(1 + (1.03961212 * lambda_array^2)/ (lambda_array^2 - 0.00600069867) + ...
    (0.231792344 .* lambda_array^2)/ (lambda_array^2 - 0.0200179144) + ...
    (1.01046945 * lambda_array^2) / (lambda_array^2 - 103.560653));

n_F2 = sqrt(1 + (1.34533359* lambda_array^2)/( lambda_array^2 - 0.00997743871) +...
    (0.209073176*lambda_array^2)/(lambda_array^2 - 0.0470450767) + ...
    (0.937357162*lambda_array^2)/(lambda_array^2 - 111.886764));


d1 = 0.7*10^(-3);
d2 = 1.15*10^(-3);
lambda1 = 445*10^(-6);
lambda2 = 739*10^(-6);
R = 35*10^(-3);
syms R3; % R3 = 21.15

for i = 1:2
    if i == 1
        lambda_array = lambda1;
        n_BK7 = sqrt(1 + (1.03961212 * lambda_array^2)/ (lambda_array^2 - 0.00600069867) + ...
        (0.231792344 .* lambda_array^2)/ (lambda_array^2 - 0.0200179144) + ...
        (1.01046945 * lambda_array^2) / (lambda_array^2 - 103.560653));

        n_F2 = sqrt(1 + (1.34533359* lambda_array^2)/( lambda_array^2 - 0.00997743871) +...
        (0.209073176*lambda_array^2)/(lambda_array^2 - 0.0470450767) + ...
        (0.937357162*lambda_array^2)/(lambda_array^2 - 111.886764));
        n_BK7 = 1;
        n_F2 = 0.99;
        M1 = [1,0;(1-n_BK7)/(n_BK7*R),1/n_BK7];
        M2 = [1,d1;0,1];
        M3 = [1,0;(n_BK7-n_F2)/(-R*n_F2),n_BK7/n_F2];
        M4 = [1,d2;0,1];
        M5 = [1,0;(n_F2-1)/(-R3),n_F2];
        Mtot1 = M5*M4*M3*M2*M1;
        C1 = Mtot1(2,1);
    end
    if i == 2
        lambda_array = lambda2;
        n_BK7 = sqrt(1 + (1.03961212 * lambda_array^2)/ (lambda_array^2 - 0.00600069867) + ...
        (0.231792344 .* lambda_array^2)/ (lambda_array^2 - 0.0200179144) + ...
        (1.01046945 * lambda_array^2) / (lambda_array^2 - 103.560653));

        n_F2 = sqrt(1 + (1.34533359* lambda_array^2)/( lambda_array^2 - 0.00997743871) +...
        (0.209073176*lambda_array^2)/(lambda_array^2 - 0.0470450767) + ...
        (0.937357162*lambda_array^2)/(lambda_array^2 - 111.886764));

        M1 = [1,0;(1-n_BK7)/(n_BK7*R),1/n_BK7];
        M2 = [1,d1;0,1];
        M3 = [1,0;(n_BK7-n_F2)/(-R*n_F2),n_BK7/n_F2];
        M4 = [1,d2;0,1];
        M5 = [1,0;(n_F2-1)/(-R3),n_F2];
        Mtot2 = M5*M4*M3*M2*M1;
        C2 = Mtot2(2,1);

    end
    
end

equation = C1 - C2 == 0;

% Solve the equation
R3_solution = solve(equation, R3);

% Display the solution
disp(solution);
%R3 = 21.15;

%% 2 d 2
N1_f2 = 1.34533359;
N2_f2 = 0.209073176;
N3_f2 = 0.937357162;
D1_f2 = -0.00997743871;
D2_f2 = -0.0470450767;
D3_f2 = -111.886764;
% 
% 
% n = sqrt(1 + N1.*lambda.^2./(lambda.^2- D1) + N2.*lambda.^2./(lambda.^2 -D2) + N3.*lambda.^2./(lambda.^2- D3));
nf2 = @(l) sqrt(1 + N1_f2.*l.^2./(l.^2- D1_f2) + N2_f2.*l.^2./(l.^2 -D2_f2) + N3_f2.*l.^2./(l.^2- D3_f2));


N1 = 1.03961212;
N2 = 0.231792344;
N3 = 1.01046945;
D1 = 0.00600069867;
D2 = 0.0200179144;
D3 = 103.560653;
nbk7 = @(l) sqrt(1 + N1*l.^2./(l.^2 - D1) + N2*l.^2./(l.^2 - D2) + N3*l.^2./(l.^2 - D3));
% WaveLength range
wavelengths_nm = 400:1000;
focal_lengths = zeros(1, length(wavelengths_nm));
% Constants 

lambda1 = 445*10^(-6);
lambda2 = 739*10^(-6);

R = 24*1e-3;
d1 = 0.62*1e-3;
d2 = 1.65*1e-3;
R3 = -0.037065717171796109168813960092515;
f_avg = 0;
n_bk7_avg = 0;
for i = 1:length(wavelengths_nm)
    lambda_f = wavelengths_nm(i) / 1000; % Convertion to appropriate units 

    n1 = nbk7(lambda_f);
    n2 = nf2(lambda_f);
    M1=[1 0;(1-n1)/(n1*R) 1/n1];
    M2 = [1 d1; 0 1];
    M3 = [1 0; (n1-n2)/n2*-R n1/n2];
    M4= [1 d2; 0 1];
    M5=[1 0;(n2-1)/R3 n2];
    M_tot = M5*M4*M3*M2*M1;
    focal_lengths(i) = -1 / M_tot(2, 1);
    f_avg = focal_lengths(i) +f_avg;
    n_bk7_avg = n_bk7_avg + n1;
end

figure(4)
plot(wavelengths_nm, focal_lengths, '--', 'MarkerSize', 100,'Color','b'); 
hold on;
xlabel('lambda nm');
ylabel('focal length mm');
title('focal length as function of lambda ');
xlabel('Wave-Length (nm)', 'Interpreter', 'latex');
ylabel('Focal length (m)', 'Interpreter', 'latex');
title('Focal Length  bk7 and f2)', 'Interpreter', 'latex');
grid on;



%% 2 d 2 המשך
syms d_before d_after
tot_mat = [1 d_before; 0 1]*M_tot*[1 d_after; 0 1];
eq1 = tot_mat(1, 2) == 0;
eq2 = tot_mat(1, 1) == 1;
sol = solve([eq1, eq2], [d_before, d_after]);
vpa(sol.d_before)
vpa(sol.d_after)


f_avg = f_avg / length(wavelengths_nm)
n_bk7_avg = n_bk7_avg / length(wavelengths_nm);
R1_new = 2*f_avg*(n_bk7_avg-1);
R2_new = -R1_new;
% Plot the focal length as a function of the wavelength


for i = 1:length(wavelengths_nm)
    lambda_f = wavelengths_nm(i) / 1000; % Convertion to appropriate units 
    n1 = nbk7(lambda_f);
    M1=[1 0;(1-n1)/(n1*R1_new) 1/n1];
    M2 = [1 d1; 0 1];
    M3 = [1 0; (n1-1)/1*R2_new n1/1];
    M_path = M3*M2*M1;
    focal_lengths(i) = -1 / M_path(2, 1);
end

% tot_mat = [1 d_before; 0 1]*M_path*[1 d_after; 0 1];
% eq1 = tot_mat(1, 2) == 0;
% eq2 = tot_mat(1, 1) == 1;
% sol = solve([eq1, eq2], [d_before, d_after]);
% vpa(sol.d_before)
% vpa(sol.d_after)
plot(wavelengths_nm, focal_lengths);
xlabel('Wave-Length (nm)', 'Interpreter', 'latex');
ylabel('Focal length (m)', 'Interpreter', 'latex');
title('Focal Length lens bk7', 'Interpreter', 'latex');
grid on;
hold off

%% 2 d 3



H1 = zeros(1, length(wavelengths_nm));
H2 = zeros(1, length(wavelengths_nm));

%%


H1avg =  0;
H2avg = 0;
for i = 1:length(wavelengths_nm)
    lambda_f = wavelengths_nm(i) / 1000; % Convertion to appropriate units 

    n1 = nbk7(lambda_f);
    n2 = nf2(lambda_f);
    M1=[1 0;(1-n1)/(n1*R) 1/n1];
    M2 = [1 d1; 0 1];
    M3 = [1 0; (n1-n2)/n2*-R n1/n2];
    M4= [1 d2; 0 1];
    M5=[1 0;(n2-1)/R3 n2];
    M_tot = M5*M4*M3*M2*M1;
    focal_lengths(i) = -1 / M_tot(2, 1);
    f_avg = focal_lengths(i) +f_avg;
    n_bk7_avg = n_bk7_avg + n1;
    A = M_tot(1,1);
    B = M_tot(1,2);
    C = M_tot(2,1);
    D = M_tot(2,2);
    H1(i)=-(1-D)/C;
    H1avg = H1(i) + H1avg;
    H2(i) = (1-A)/C;
    H2avg = H2(i) + H2avg;

end

H1avg = H1avg / length(wavelengths_nm)
H2avg = H2avg / length(wavelengths_nm)


figure(5)
plot(wavelengths_nm, H1, '--', 'MarkerSize', 100,'Color','b'); 
hold on;
xlabel('lambda nm');
ylabel('focal length mm');
title('focal length as function of lambda ');
xlabel('Wave-Length (nm)', 'Interpreter', 'latex');
ylabel('Focal length (m)', 'Interpreter', 'latex');
title('H1 as a function of Wavelength', 'Interpreter', 'latex');
grid on;

figure(4)
plot(wavelengths_nm, H2, '--', 'MarkerSize', 100,'Color','b'); 
hold on;
xlabel('lambda nm');
ylabel('focal length mm');
title('focal length as function of lambda ');
xlabel('Wave-Length (nm)', 'Interpreter', 'latex');
ylabel('Focal length (m)', 'Interpreter', 'latex');
title('H2 as a function of Wavelength', 'Interpreter', 'latex');
grid on;
hold off


%% function
function A = matrixA(d1)
    A = [1,d1;0,1];
end

function B = matrixB(R,n_BK7)
    B = [1,0;(1-n_BK7)/(n_BK7*R),1/n_BK7];
end

function C = matrixC(R,n_BK7)
    C = [1,0;(n_BK7-1)/(-R),n_BK7];
end

function D = matrixD(d2)
    D = [1,d2;0,1];
end


% Multiply the matrices to obtain M
function M = calculateM(d1,d2,n_BK7,R)
    A = matrixA(d1);
    B = matrixB(R,n_BK7);
    C = matrixC(R,n_BK7);
    D = matrixD(d2);
    M = D * C * B * A;
end








