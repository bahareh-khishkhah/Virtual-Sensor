% Improving the Quality of Photoacoustic Images 
% Using Virtual Sensor Points
% by bahareh khishkhah

clc
clear all
close all

% =========================================================================
% SIMULATION
% =========================================================================

% load the initial pressure distribution from an image and scale
p0_magnitude = 2;
p0 = p0_magnitude * loadImage('196.bmp');

% assign the grid size and create the computational grid
PML_size = 20;              % size of the PML in grid points
Nx = 256 - 2 * PML_size;    % number of grid points in the x direction
Ny = 256 - 2 * PML_size;    % number of grid points in the y direction
x = 10e-3;                  % total grid size [m]
y = 10e-3;                  % total grid size [m]
dx = x / Nx;                % grid point spacing in the x direction [m]
dy = y / Ny;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% resize the input image to the desired number of grid points
p0 = resize(p0, [Nx, Ny]);

% smooth the initial pressure distribution and restore the magnitude
% p0 = smooth(p0, true);
p0(p0~=0)=1;% assign to the source structure
source.p0 = p0;

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

%plot the initial pressure and sensor distribution
figure;
imshow(p0);
colorbar;
% ylabel('x-position [mm]');
% xlabel('y-position [mm]');
axis image;
title('the initial pressure');

% define the frequency response of the sensor elements
center_freq = 5e6;      % [Hz]
bandwidth = 90;         % [%]
sensor.frequency_response = [center_freq, bandwidth];
% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% User input for sensor type selection
sensor_type = input('Choose sensor type (1 for Linear, 2 for Circular): ');
% User input for number of sensors
number_sensors = input('Choose number of sensors (note: the number of sensor points should not exceed the number of grid points): ');
% Sensor Configuration Based on User Input

if sensor_type == 1
    % Linear Sensor Configuration
    num_sensors = number_sensors; % Number of sensors for linear array
    sensor_size = 1; % Size of each sensor in pixels
    interval = (Nx/num_sensors)-1; % Interval between sensors in pixels

    % Initialize mask matrix for linear sensors
    sensor.mask = zeros(kgrid.Nx, kgrid.Ny); 

    % Define linear sensor positions
    for i = 1:num_sensors
        x_start = 1;
        x_end = sensor_size;
        y_start = (i - 1) * (sensor_size + interval) + 1;
        y_end = y_start + sensor_size - 1;
        sensor.mask(x_start:x_end, y_start:y_end) = 1;
    end

    % Plot initial pressure and linear sensor distribution
    figure;
    imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0 + sensor.mask);
    colorbar;
    title('Initial Pressure with Linear Sensor Distribution');
    axis image;

elseif sensor_type == 2
    % Circular Sensor Configuration
    sensor_radius = 4.5e-3; % [m]
    num_sensor_points_circular = number_sensors;

    cart_sensor_mask_circular = makeCartCircle(sensor_radius, num_sensor_points_circular, [0, 0], pi);
    
    % Assign to sensor structure
    sensor.mask = cart_sensor_mask_circular;

    % Plot initial pressure and circular sensor distribution
    figure;
    imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0 + cart2grid(kgrid, cart_sensor_mask_circular), [-1, 1]);
    colorbar;
    title('Initial Pressure with Circular Sensor Distribution');
    xlabel('x-position [mm]');
    ylabel('y-position [mm]');
    axis image;

else
    error('Invalid selection. Please choose either 1 for Linear or 2 for Circular.');
end

% set the input arguments
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, 'Smooth', false, 'PlotSim', false,};

% run the simulation for omnidirectional detector elements
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Step 1: Data preparation
% Import data
timeSeriesData =sensor_data;
% Removing intermediate sensors
timeSeriesData=timeSeriesData(1:2:end, :);

[rows, columns] = size(timeSeriesData);
num_time_steps =columns;

rows1=rows-1;
% Preallocate Kalman_Data
Kalman_Data = zeros(rows1, num_time_steps);
% Preallocate arrays for plotting
Kalman_gain = zeros(rows1, num_time_steps);
Error_covariance = zeros(rows1, num_time_steps);

for k = 1:rows1
    % Extract X and Y values from the data
    X = timeSeriesData(k,:);
    Y = timeSeriesData(k+1,:);

    % Step 2: Feature Extraction
    % Given statistical features for Time Series X
    X_mean = mean(X);
    X_std_dev = std(X);

    % Given statistical features for Time Series Y
    Y_mean = mean(Y);
    Y_std_dev = std(Y);

    % Initialize the state vector and covariance matrix
    initial_guess = (X_mean + Y_mean) / 2;
    initial_covariance = max(X_std_dev^2, Y_std_dev^2);

    state_estimate = initial_guess;
    covariance_estimate = initial_covariance;

    % Define the measurement function
    H = @(state) (X + Y) / 2;

    % Initialize sensor_measurements using statistical features
    sensor_measurements = (X + Y) / 2;

    % Measurement noise covariance
    measurement_residual_covariance =1e-2;%10^2  
    
    dt = 1;  % Time step
    t = 0:dt:(num_time_steps-1);  % Time vector
    n = length(t);
    

    % Perform data assimilation using Kalman filter
    for i = 1:n
        % Prediction step
        state_predict = state_estimate;
        covariance_predict = covariance_estimate;
        % Calculate the predicted measurements for the ensemble
        predicted_measurement = H(state_predict);

        % Update step
        kalman_gain = covariance_predict / (covariance_predict + measurement_residual_covariance);
        measurement_residual = sensor_measurements(i) - H(state_predict);

        state_estimate = state_predict + kalman_gain * measurement_residual;
        covariance_estimate = (1 - kalman_gain) * covariance_predict;
        % Calculate the predicted measurements for the ensemble
        predicted_measurement1 = H(state_estimate);
        % Calculate the Kalman gain and error covariance matrix
        K = kalman_gain;
        P = covariance_estimate;



    end
        % Store the result in Kalman_Data
        Kalman_Data(k, :) = predicted_measurement1;
        % Store the Kalman gain and error covariance matrix for plotting
        Kalman_gain(k, :) = K;
        Error_covariance(k, :) = diag(P);
        
        % Calculate the lower and upper bounds
        z_forecast = predicted_measurement1(1:end);
        C = [1 0];  % Measurement matrix
        P_est = covariance_estimate;
        R = measurement_residual_covariance;
        lower_bound = z_forecast - 2 * sqrt(C * P_est * C' + R);
        upper_bound = z_forecast + 2 * sqrt(C * P_est * C' + R);


end

% combine new data
data1 = sensor_data(1:2:end, :);
data2 = Kalman_Data;

combinedData1 = sensor_data(end,:);
data2 = [data2; combinedData1];

[rows1, columns1] = size(sensor_data);
num_time_steps =columns1;
% Create a matrix with rows*columns dimensions and fill with 0s
matrix = zeros(rows1, columns1);

% Populate even rows with values from file2.csv and odd rows with values from file1.csv
for i = 1:rows1
    if mod(i, 2) == 0
        matrix(i, :) = data2(i/2, :);
    else
        matrix(i, :) = data1((i+1)/2, :);
    end
end
% virtual sensor data with kalman fiter
NEw_Sensordata=matrix;
% reset the initial pressure
source.p0 = 0;

% with kalman data
% assign the time reversal data for the virtual case
sensor.time_reversal_boundary_data = NEw_Sensordata;

% run the time reversal reconstruction
NEW_p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
NEW_p0_recon(NEW_p0_recon < 0) = 0;
NEW_p0_recon=mat2gray(NEW_p0_recon);
NEW_p0_recon=imadjust(NEW_p0_recon);

% original sensor data for Comparison with virtual sensor data
% with half sensor points
% Sensor Configuration Based on User Input
if sensor_type == 1
    % Linear Sensor Configuration
    num_sensors = number_sensors/2; % Number of sensors for linear array
    sensor_size = 1; % Size of each sensor in pixels
    interval = (Nx/num_sensors)-1; % Interval between sensors in pixels

    % Initialize mask matrix for linear sensors
    sensor.mask = zeros(kgrid.Nx, kgrid.Ny); 

    % Define linear sensor positions
    for i = 1:num_sensors
        x_start = 1;
        x_end = sensor_size;
        y_start = (i - 1) * (sensor_size + interval) + 1;
        y_end = y_start + sensor_size - 1;
        sensor.mask(x_start:x_end, y_start:y_end) = 1;
    end

    % Plot initial pressure and linear sensor distribution
    figure;
    imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0 + sensor.mask);
    colorbar;
    title('Initial Pressure with Linear Sensor Distribution');
    axis image;

elseif sensor_type == 2
    % Circular Sensor Configuration
    sensor_radius = 4.5e-3; % [m]
    num_sensor_points_circular = number_sensors/2;

    cart_sensor_mask_circular = makeCartCircle(sensor_radius, num_sensor_points_circular, [0, 0], pi);
    
    % Assign to sensor structure
    sensor.mask = cart_sensor_mask_circular;

    % Plot initial pressure and circular sensor distribution
    figure;
    imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0 + cart2grid(kgrid, cart_sensor_mask_circular), [-1, 1]);
    colorbar;
    title('Initial Pressure with Circular Sensor Distribution');
    xlabel('x-position [mm]');
    ylabel('y-position [mm]');
    axis image;

else
    error('Invalid selection. Please choose either 1 for Linear or 2 for Circular.');
end
sensor.time_reversal_boundary_data = timeSeriesData;

% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
p0_recon(p0_recon < 0) = 0;
p0_recon=mat2gray(p0_recon);

% Calculate SSIM for the original and Kalman-filtered reconstructions
ssim_original = ssim(p0_recon, p0);
ssim_kalman = ssim(NEW_p0_recon, p0);

% Display the SSIM values
disp(['SSIM for the original sensor data reconstructed image: ', num2str(ssim_original)]);
disp(['SSIM for the virtual sensor data reconstructed image: ', num2str(ssim_kalman)]);

% Calculate the percentage improvement in SSIM
percentage_improvement = ((ssim_kalman - ssim_original) / ssim_original) * 100;

% Display the percentage improvement
disp(['Percentage Improvement in SSIM: ', num2str(percentage_improvement), '%']);

% Calculate PSNR for the original and Kalman-filtered reconstructions
psnr_original = psnr(p0_recon, p0);
psnr_kalman = psnr(NEW_p0_recon, p0);

% Display the PSNR values
disp(['PSNR for tthe original sensor data reconstructed image: ', num2str(psnr_original)]);
disp(['PSNR for the virtual sensor data reconstructed image: ', num2str(psnr_kalman)]);

% Calculate the percentage improvement in PSNR
percentage_improvement = ((psnr_kalman - psnr_original) / psnr_original) * 100;

% Display the percentage improvement
disp(['Percentage Improvement in PSNR: ', num2str(percentage_improvement), '%']);

% Set the desired figure size (in pixels)
figureWidth = 512;  % Width of the figure
figureHeight = 512; % Height of the figure
set(gcf, 'Units', 'pixels', 'Position', [100, 100, figureWidth, figureHeight]);

% Plot the reconstructed images in one figure without axes
figure; 

% Create a 1x2 subplot layout
subplot(1, 2, 1);
imagesc(p0_recon); 
title('Reconstructed Image from Original Sensor Data');
axis off; % Turn off axis

subplot(1, 2, 2);
imagesc(NEW_p0_recon); 
title('Reconstructed Image from Virtual Sensor Data');
axis off; % Turn off axis

% Display SSIM and PSNR values
disp(['SSIM for the original sensor data reconstructed image: ', num2str(ssim_original)]);
disp(['SSIM for the virtual sensor data reconstructed image: ', num2str(ssim_kalman)]);
disp(['Percentage Improvement in SSIM: ', num2str(percentage_improvement), '%']);

disp(['PSNR for the original sensor data reconstructed image: ', num2str(psnr_original)]);
disp(['PSNR for the virtual sensor data reconstructed image: ', num2str(psnr_kalman)]);
disp(['Percentage Improvement in PSNR: ', num2str(percentage_improvement), '%']);


