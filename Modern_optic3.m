%% 1 b + c
clear 
clc
% Define radius of the circular aperture
radius = (mod(445,5)+1)*10^-2;
L = 0.2;
N = 200; % grid size
n = 1:1:N;
dx = L/N;
x_n = n*dx;
y_n = x_n;
[x1, y1] = meshgrid(linspace(-L/2, L/2 - dx , N));

% Generate circular aperture function
aperture = circ(x1, y1, radius);

% Display the aperture as an image
imagesc(x1(1,:), y1(:,1), aperture);
colormap(gray);  % Set colormap to grayscale
axis image;  % Set aspect ratio to equal
xlabel('X (m)');
ylabel('Y (m)');
title('Circular Aperture');
%%
syms x
equ = 2*x^2 == radius^2;
sol = solve(equ,x);
xy_max = abs(sol);
%% 1 d
lambda = 800e-9;
K = 2*pi/lambda;
z_0 =5*(K*(radius^2)/2);
M = 1/dx;
k = 1:1:N;
df = 1/L;
f_k = k*df;
fx=-1/(2*dx):1/L:1/(2*dx)-1/L; % Spatial frequency vector (in 1/meters)
fy = fx;
[Fx,Fy] = meshgrid(linspace(-M/2,M/2 - df,N));
aperture_dft = F(aperture);
%aperture_dft = fftshift(fft2(aperture));


% Display the aperture as an image
imagesc(Fx(1,:),Fy(:,1), abs(aperture_dft));
colormap('jet');
axis equal;
xlabel('Spatial Frequency (1/m)');
ylabel('Spatial Frequency (1/m)');
title('DFT Magnitude');
colorbar;

%%
% Display the magnitude of the DFT using surf
subplot(1, 2, 2);
surf((-N/2:N/2-1)/(L), (-N/2:N/2-1)/(L), abs(aperture_dft));
colormap('jet');
axis equal;
xlabel('Spatial Frequency (1/m)');
ylabel('Spatial Frequency (1/m)');
zlabel('Magnitude');
title('DFT Magnitude - surf');
camlight left;  % Add light from the left
lighting phong;  % Set lighting model to phong
shading interp;  % Interpolated shading for smoother appearance

%% 2 a
lambda = 800e-9;
K = 2*pi/lambda;
R = 0.01;
% franofer condition
z_0 = 5*(K*(R^2)/2);

%% 2 b + c

%frahanufer
H_fra = 1;
E_out_fra = F(aperture.*H_fra);
x_prop=lambda*z_0*fx; %axis
% Display the aperture as an image
imagesc(x_prop,x_prop,abs(E_out_fra));
colormap("jet");  % Set colormap to grayscale
axis image;  % Set aspect ratio to equal
xlabel('x2');
ylabel('y2');
title('Eout with Frahanufer');
%%
surf(abs(E_out_fra))

%% 2 d
% Choose a specific line in the image
line_index = N / 2; % Choose the middle line

% Extract the intensity values along the chosen line
intensity_profile = abs(E_out_fra(line_index, :));

% Plot the intensity profile
figure;
plot(x_prop, intensity_profile, 'k', 'LineWidth', 2);
xlabel('Position along x (\mm)');
ylabel('Intensity');
title('Intensity Profile along a Cross-section');
%%
 % Calculate intensity profile along x-axis
intensity_profile = abs(E_out_fra(N/2, :)).^2;

% Find maximum intensity
[max_intensity, max_index] = max(intensity_profile);

% Find half-maximum intensity
half_max_intensity = max_intensity / 2;

% Find positions where intensity equals half-maximum on both sides of the maximum
left_index = find(intensity_profile(1:max_index) <= half_max_intensity, 1, 'last');
right_index = find(intensity_profile(max_index:end) <= half_max_intensity, 1) + max_index - 1;

% Calculate FWHM
FWHM = x_prop(right_index) - x_prop(left_index);

% Plot intensity profile
figure;
plot(x_prop*1e6, intensity_profile, 'LineWidth', 2);
hold on;

plot([x_prop(left_index), x_prop(right_index)]*1e6, [half_max_intensity, half_max_intensity], 'r', 'LineWidth', 2);
xlabel('x (mm)');
ylabel('Intensity');
title(['|Intensity| Profile along x-axis at z = z_0']);
legend('|Intensity| Profile', 'Half Maximum');
grid on;
%% Q 3 b + c + d
%H = exp(1i*(K/(2*z_0)).*(x1(1,:).^2+y1(:,1).^2));
%[aperture_dft,x_prop] = discrete_prop(aperture,lambda,z,fx,K,x1,y1);

z1 = z_0;
z = [0.01*z1, 0.1*z1, 0.5*z1, z1, 2*z1, 10*z1];

% Create a subplot with 2 rows and 6 columns
figure;
for i = 1:3
    subplot(2, 3, i); % Specify subplot position
    H = exp(1i*(K/(2*z(i))).*(x1.^2+y1.^2));
    E_out_fre = (abs(F(aperture.*H)));
    % Display the magnitude of the Fresnel transform
    imagesc(x_prop, x_prop, abs(E_out_fre));
    colormap('jet');
    axis equal;
    xlabel('x2 mm');
    ylabel('y2 mm');
    title(['Fresnel Magnitude, z = ', num2str(z(i))]);
    colorbar;
end
for i = 4:6
    subplot(2, 3, i); % Specify subplot position
    H = exp(1i*(K/(2*z(i))).*(x1.^2+y1.^2));
    E_out_fre = (abs(F(aperture.*H)));
    % Display the magnitude of the Fresnel transform
    imagesc(x_prop, x_prop, abs(E_out_fre));
    colormap('jet');
    axis equal;
    xlabel('x2 mm');
    ylabel('y2 mm');
    title(['Fresnel Magnitude, z = ', num2str(z(i))]);
    colorbar;
end


%% check difference between Franhufer and fresnel
zc = z_0;
H = exp(1i*(K/(2*zc)).*(x1.^2+y1.^2));
E_out_fre = F(aperture.*H);
E_out_fra = F(aperture);
subplot(1,2,1);
imagesc(x_prop, x_prop, abs(E_out_fre));
subplot(1,2,2); 
imagesc(x_prop, x_prop, abs(E_out_fra));



 % Calculate intensity profile along x-axis
intensity_fre = abs(E_out_fre(N/2, :));%.^2;

% Find maximum intensity
[max_intensity, max_index] = max(intensity_fre);

% Find half-maximum intensity
half_max_intensity = max_intensity / 2;

% Find positions where intensity equals half-maximum on both sides of the maximum
left_index = find(intensity_fre(1:max_index) <= half_max_intensity, 1, 'last');
right_index = find(intensity_fre(max_index:end) <= half_max_intensity, 1) + max_index - 1;

% Calculate FWHM
FWHM = x_prop(right_index) - x_prop(left_index);

% Plot intensity profile
figure;
plot(x_prop*1e6, intensity_fre, 'LineWidth', 2);
hold on;

plot([x_prop(left_index), x_prop(right_index)]*1e6, [half_max_intensity, half_max_intensity], 'r', 'LineWidth', 2);
xlabel('x (mm)');
ylabel('Intensity');
title(['|E out| Profile along x-axis at z = 10z_0']);
legend('|E out| Profile', 'Half Maximum');
grid on;

%%
z1 = z_0;
z = [0.01*z1, 0.1*z1, 0.5*z1, z1, 2*z1, 10*z1];

% Initialize arrays to store intensity data
intensity_data = zeros(length(x_prop), 6);

for i = 1:6
    H = exp(1i*(K/(2.*z(i))).*(x1.^2+y1.^2));
    E_out_fre = F(aperture.*H);
    intensity_fre = abs(E_out_fre(N/2, :)).^2;
    intensity_data(:, i) = intensity_fre; % Store intensity data
end

% Plot all intensity profiles on the same axes
figure;
hold on;
for i = 1:6
    plot(x_prop*1e6/z_0, intensity_data(:, i), 'LineWidth', 2);
    xlabel('x (mm)');
    ylabel('Intensity');
    title(['|Intensity| Profile along x-axis']);
    legend('|Intensity| Profile');
grid on;
end
hold off;

