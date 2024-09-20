clear
clc
clf
% Check parallel pool status
poolObj = gcp('nocreate'); % Gets current pool without creating a new one if it doesn't exist
tic;  % Start timing
cmap = colororder();

% Load the fitted kernel density estimation probability density functions of the process noise
load('kdeFunction.mat');
load('kdeFunction_w.mat');

%% Simulation of a model for generating hypothetical measurements

% Process and measurement noise parameters

A_process_1 = 0;
A_measurement_1 = 0;

B_process_1 = 0.05;
B_measurement_1 = 0.05;

name = 'Logistic'; % Type of distribution of the noise (Non-Gaussian)


i=0; x(:,1)=[0;0];
for ti=tsim-2*dt
    i=i+1;

    er_process(:,i+1) = [random(name, A_process_1, B_process_1); random(name, A_process_1, 10*B_process_1)]; % Additive Laplace process noise for both variables

    er(:,i+1) = [random(name, A_measurement_1, B_measurement_1); random(name, A_measurement_1, 10*B_measurement_1)]; % Additive Laplace measurement noise for both variables


    u = 5*sin(0.2*ti);
    uc(:,i+1) = u;
    x(:,i+1) = f_model(x(:,i),wf(:,i),u,dt) + er_process(:,i+1); % wf is responsbile for injecting the faults | x is the combined state vector for I_dot and omega_dot
    % The codes for f_model and faults are not provided as do not have permission to make them public but you can use any dynamic system. The faults are injected locally but they are not necessary for an estimation problem so they are not necessary as well.
    z(:,i+1) = x(:,i+1) + er(:,i+1);
    x_true(:,i+1) = f_model(x_true(:,i),wf(:,i),u,dt);
end

X_real=x(1,:); % only I
X_real_w=x(2,:); % only I
X_real_true=x_true(1,:); % only I
Z_real=z(1,:); % only I
Z_real_w=z(2,:); % only w



subplotXmanyY_er(7,1)
plot(Z_real, 'color', cmap(3,:))
ylim([-1.5 1.5])
xlim([0 1050])
set(gca, 'FontSize', 11); % Set the font size of the current axes
xlabel('Time Step', 'FontSize', 12);
ylabel('$I (A)$', 'interpreter', 'latex', 'FontSize', 15);


prompt = "Do you want to continue (y/n)? ";  % So the code up to this point, generates a signal output and this prompt makes sure if you want to perform state estimation based on this measurement or sensor signal.
inputStr = input(prompt, 's');  % Treat the input as a string
if inputStr == 'n'
    disp('Exiting...');
    return;
end

%% Forming the importance density function



data = zeros(1, 1000);
for xii = 1:1000
    data(xii) = (1-0.5)*random(name, A_process_1, 1*B_process_1) - 0.5*random(name, A_measurement_1, 1*B_measurement_1);
end

xi_kde = linspace(min(data), max(data), 1000); % 1000 evenly spaced points over the range
[f_kde, xi_kde] = ksdensity(data);

integral_kde_q = trapz(xi_kde, f_kde);

kde_pdf_function_q = @(x) interp1(xi_kde, f_kde, x, 'linear', 'extrap');


data_w = zeros(1, 1000);
for xii = 1:1000
    data_w(xii) = (1-0.5)*random(name, A_process_1, 1*10*B_process_1) - 0.5*random(name, A_measurement_1, 1*10*B_measurement_1);
end

xi_kde_w = linspace(min(data_w), max(data_w), 1000); % 1000 evenly spaced points over the range
[f_kde_w, xi_kde_w] = ksdensity(data_w);

integral_kde_q_w = trapz(xi_kde_w, f_kde_w);

kde_pdf_function_q_w = @(x) interp1(xi_kde_w, f_kde_w, x, 'linear', 'extrap');


% Plotting the KDE's:

% x_range = linspace(min(data) - 1, max(data) + 1, 1000);  % Extend the range slightly beyond the min and max of the data
% kde_values = kde_pdf_function_q(x_range);
% x_range_w = linspace(min(data_w) - 1, max(data_w) + 1, 1000);  % Extend the range slightly beyond the min and max of the data
% kde_values_w = kde_pdf_function_q_w(x_range_w);


% figure;  % Create a new figure window
% plot(x_range, kde_values, 'LineWidth', 2);  % Plot the KDE
% title('Kernel Density Estimate (I)');
% xlabel('Data Values');
% ylabel('Density');
% grid on;  % Turn on the grid
% 
% 
% figure;  % Create a new figure window
% plot(x_range_w, kde_values_w, 'LineWidth', 2);  % Plot the KDE
% title('Kernel Density Estimate (w)');
% xlabel('Data Values');
% ylabel('Density');
% grid on;  % Turn on the grid



% Define system parameters
numParticles = 100;  % Number of particles
N_eff_ratio = 1;
N_eff_thresh = floor(numParticles * N_eff_ratio);
numModes = 9;  % Number of modes

% Initialize particles
particles = rand(1, numParticles);  % Random initial state
particles_w = rand(1, numParticles);  % Random initial state
mode = randi([1, numModes], 1, numParticles);  % Initial mode for each particle


% Initialization with specific state range
initialStateRange = [-1*0.01, 1*0.01];  % Initial state range (I)
initialStateRange_w = [-6*0.01, 6*0.01];  % Initial state range (omega or w)

% Initialize particles' states
particles = initialStateRange(1) + (initialStateRange(2) - initialStateRange(1)) * rand(1, numParticles);
particles_history = zeros(T, numParticles);
particles_history(1, :) = particles;
particles_history_resamp = zeros(T, numParticles); % Storing the computed states for each particle accounting for their resamplings. The other history variable stores the histories that are compatible with the significant weights that got resampled for getting fed into the model.
particles_history_resamp(1, :) = particles;

particles_w = initialStateRange_w(1) + (initialStateRange_w(2) - initialStateRange_w(1)) * rand(1, numParticles);
particles_history_w = zeros(T, numParticles);
particles_history_w(1, :) = particles_w;
particles_history_resamp_w = zeros(T, numParticles); % Storing the computed states for each particle accounting for their resamplings. The other history variable stores the histories that are compatible with the significant weights that got resampled for getting fed into the model.
particles_history_resamp_w(1, :) = particles_w;

% Initialize all particles to start in mode one, reflecting 100% certainty
mode = ones(1, numParticles);  % Mode one for all particles
modeHistory = zeros(T, numParticles);
modeHistory(1, :) = mode; % Store the modes for the current timestep
modeHistory_resamp = zeros(T, numParticles); % Storing the modes for each particle accounting for their resamplings. The other history variable stores the histories that are compatible with the significant weights that got resampled for getting fed into the model.
modeHistory_resamp(1, :) = mode; % Store the modes for the current timestep

% Storing the final estimated states of the particle filter:
estimatedState_history = zeros(1,T);
estimatedState_history(1) = 0;

estimatedState_history_w = zeros(1,T);
estimatedState_history_w(1) = 0;

N_eff_history = zeros(1,T);
N_eff_history(1, 1) = numParticles;
N_eff_history_w = zeros(1,T);
N_eff_history_w(1, 1) = numParticles;

% Initialize weights uniformly
weights = ones(1, numParticles) / numParticles;
parent_weights = ones(1, numParticles) / numParticles;
weights_history = zeros(T, numParticles);
weights_history(1, :) = weights;

integral_result = zeros(1, numParticles);
integral_result_w = zeros(1, numParticles);

weights_history_plt = zeros(T, numParticles);
weights_history_plt(1, :) = weights;

weights_w = ones(1, numParticles) / numParticles;
parent_weights_w = ones(1, numParticles) / numParticles;
weights_history_w = zeros(T, numParticles);
weights_history_w(1, :) = weights_w;

weights_history_w_plt = zeros(T, numParticles);
weights_history_w_plt(1, :) = weights_w;


% Main filtering loop
for t = 2:1:T  % Assume T is the number of time steps
    measurement = getMeasurementAtTime(t, Z_real);  % Function to get measurement
    measurement_w = getMeasurementAtTime(t, Z_real_w);
    prior_output_all = zeros(2, numParticles);
    %prior_output = zeros(2, numParticles);


    % Parent weights calculations:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor p = 1:numParticles


          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the importance density function for drawing the particles(p) from
        % or sampling them:

        % Parameters for the logistic distribution
        [prior1, prior2] = model_output(t, mode(p), particles_history, particles_history_w, modeHistory, p);  % I do not have permission to make the model_output function public but you ban substitute any dynamic system simulator with it with the caveat that this returns on time-step's worth of propagation for system dynamics
        % The modes refer to mode of the system at each time-step (if your dynamic system is multi-modal, if not, you can set all your time-step modes to be one value or remove it entirely from this code)
        prior_output = [prior1; prior2];


        particles(p) = (1-0.5)*(prior_output(1,:) + random(name, A_process_1, 1*B_process_1)) + 0.5*(measurement - random(name, A_measurement_1, 1*B_measurement_1));
        particles_w(p) = (1-0.5)*(prior_output(2,:) + random(name, A_process_1, 1*10*B_process_1)) + 0.5*(measurement_w - random(name, A_measurement_1, 1*10*B_measurement_1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        % the importance density function for the weight update stage

        prior_output = prior1;
        prior_output_w = prior2;




        likelihood = pdf(name, measurement, particles(p) + A_measurement_1, B_measurement_1);
        likelihood_w = pdf(name, measurement_w, particles_w(p) + A_measurement_1, 10*B_measurement_1);

        % the process or prior distribution
        prior_x_i_minus_1 = prior_output;
        prior_x_i_minus_1_w = prior_output_w;

        prior = kde_likelihood(particles(p), prior_x_i_minus_1(1,:), kde_pdf_function);
        prior_w = kde_likelihood(particles_w(p), prior_x_i_minus_1_w(1,:), kde_pdf_function_w);


        % Calculating the integrals for the parent weights (For auxiliary particle filters):

        I_n_min = -1; % Adjust based on your data
        I_n_max = 1;  % Adjust based on your data
        num_points = 200;  % More points for better accuracy
        I_n_values = linspace(I_n_min, I_n_max, num_points);

        prior_values = zeros(1, num_points);
        likelihood_values = zeros(1, num_points);
        products = zeros(1, num_points);


        for i = 1:num_points
            I_n = I_n_values(i);

            % Compute prior using KDE
            % Assuming kde_likelihood is a function you have defined or available
            prior_values(i) = kde_likelihood(I_n, prior_x_i_minus_1(1,:), kde_pdf_function);

            % Compute likelihood
            % Assuming you have a function or a way to compute the logistic PDF
            % Modify this with the correct function to compute the  logistic PDF
            likelihood_values(i) = pdf(name, measurement, I_n + A_measurement_1, B_measurement_1);

            % Compute the product of prior and likelihood
            products(i) = prior_values(i) * likelihood_values(i);
        end

        % Numerically integrate using the trapezoidal rule
        integral_result(p) = trapz(I_n_values, products);

        parent_weights(p) = weights_history(t-1, p)*integral_result(p);


        w_n_min = -10; % Adjust based on your data
        w_n_max = 10;  % Adjust based on your data
        num_points = 200;  % More points for better accuracy
        w_n_values = linspace(w_n_min, w_n_max, num_points);

        prior_values = zeros(1, num_points);
        likelihood_values = zeros(1, num_points);
        products = zeros(1, num_points);


        for i = 1:num_points
            w_n = w_n_values(i);

            % Compute prior using KDE
            % Assuming kde_likelihood is a function you have defined or available
            prior_values(i) = kde_likelihood(w_n, prior_x_i_minus_1_w(1,:), kde_pdf_function_w);

            % Compute likelihood
            % Assuming you have a function or a way to compute logistic PDF
            % Modify this with the correct function to compute logistic PDF
            likelihood_values(i) = pdf(name, measurement_w, w_n + A_measurement_1, 10*B_measurement_1);

            % Compute the product of prior and likelihood
            products(i) = prior_values(i) * likelihood_values(i);
        end

        % Numerically integrate using the trapezoidal rule
        integral_result_w(p) = trapz(w_n_values, products);

        %parent_integral = integral_result * integral_result_w;

        parent_weights_w(p) = weights_history_w(t-1, p)*integral_result_w(p);


%         % weight update
%         weights(p) = likelihood*prior/importance_density/integral_result;
% 
%         weights_w(p) = likelihood_w*prior_w/importance_density_w/integral_result_w;
%         %weights(p) = likelihood*weights_history(t-1, p);
    end

    parent_weights = parent_weights / sum(parent_weights);  % Normalize parent weights
    parent_weights_w = parent_weights_w / sum(parent_weights_w);  % Normalize parent weights
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sampling the parent trajectories based on parent weights:
    % (Presampling instead of resampling)

    [particles, mode, parent_weights, particles_history, modeHistory, integral_result] = presample(particles, mode, parent_weights, numParticles, particles_history, modeHistory, integral_result);
    [particles_w, parent_weights_w, particles_history_w, integral_result_w] = presample_w(particles_w, parent_weights_w, numParticles, particles_history_w, integral_result_w);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    % Prediction step
    parfor p = 1:numParticles
        currentState = particles(p);
        currentMode = mode(p);

 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the importance density function for drawing the particles(p) from
        % or sampling them:

        % Parameters for the logistic distribution
        [prior1, prior2] = model_output(t, mode(p), particles_history, particles_history_w, modeHistory, p);
        prior_output = [prior1; prior2];


        particles(p) = (1-0.5)*(prior_output(1,:) + random(name, A_process_1, 1*B_process_1)) + 0.5*(measurement - random(name, A_measurement_1, 1*B_measurement_1));
        particles_w(p) = (1-0.5)*(prior_output(2,:) + random(name, A_process_1, 1*10*B_process_1)) + 0.5*(measurement_w - random(name, A_measurement_1, 1*10*B_measurement_1));

        prior_output_all(:, p) = prior_output;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Mode transition (if your dynamic system is multi-modal, if not, you can set all your time-step modes to be one value or remove it entirely from this code)
        mode(p) = ;  % Sample new mode
    end


    particles_history(t, :) = particles;
    particles_history_resamp(t, :) = particles;
    particles_history_w(t, :) = particles_w;
    particles_history_resamp_w(t, :) = particles_w;
    modeHistory(t, :) = mode; % Store the modes for the current timestep
    modeHistory_resamp(t, :) = mode;





    % Measurement update (simplified for one state)
    
    %weights = zeros(1, numParticles);
    parfor p = 1:numParticles

    
        % the importance density function for the weight update stage

        prior_output = prior_output_all(1, p);
        prior_output_w = prior_output_all(2, p);

        importance_density = kde_likelihood(particles(p), (1-0.85)*prior_output(1,:) + 0.85*measurement, kde_pdf_function_q);
        importance_density_w = kde_likelihood(particles_w(p), (1-0.5)*prior_output_w(1,:) + 0.5*measurement_w, kde_pdf_function_q_w);


        likelihood = pdf(name, measurement, particles(p) + A_measurement_1, B_measurement_1);
        likelihood_w = pdf(name, measurement_w, particles_w(p) + A_measurement_1, 10*B_measurement_1);

        % the process or prior distribution
        prior_x_i_minus_1 = prior_output;
        prior_x_i_minus_1_w = prior_output_w;

        prior = kde_likelihood(particles(p), prior_x_i_minus_1(1,:), kde_pdf_function);
        prior_w = kde_likelihood(particles_w(p), prior_x_i_minus_1_w(1,:), kde_pdf_function_w);


        % weight update
        weights(p) = likelihood*prior/importance_density/integral_result(p);

        weights_w(p) = likelihood_w*prior_w/importance_density_w/integral_result_w(p);
        %weights(p) = likelihood*weights_history(t-1, p);
    end

    weights = weights / sum(weights);  % Normalize weights
    weights_w = weights_w / sum(weights_w);  % Normalize weights
    

    weights_history_plt(t, :) = weights;
    weights_history_w_plt(t, :) = weights_w;

     % Optional resampling step can be added here
     N_eff = 1/(sum(weights.^2, 'all'));
     N_eff_w = 1/(sum(weights_w.^2, 'all'));
%     if N_eff > N_eff_thresh
%         continue
%     else
%         [particles, mode, weights, particles_history, particles_history_w, modeHistory, weights_w] = resample(particles, mode, weights, numParticles, particles_history, particles_history_w, modeHistory, weights_w);
%     end

    weights_history(t, :) = weights;
    weights_history_w(t, :) = weights_w;

    % Estimate state
    estimatedState = sum(particles .* weights);
    estimatedState_history(t) = estimatedState;

    estimatedState_w = sum(particles_w .* weights_w);
    estimatedState_history_w(t) = estimatedState_w;

    N_eff_history(1, t) = N_eff;
    N_eff_history_w(1, t) = N_eff_w;


    % Continue to next time step...
    fprintf('Timestep: %d\n', t);
end

elapsedTimeSeconds = toc;  % End timing



subplotXmanyY_er(7,2)
plot(estimatedState_history, 'color', cmap(3,:))
hold on
plot(X_real, 'color', cmap(4,:))
hold on
plot(Z_real, 'color', cmap(5,:))
legend('Estimate', 'True state', 'Measurement');
ylim([-1.5 1.5])
xlim([0 1050])
xlabel('Time Step', 'FontSize', 12);
ylabel('$I (A)$', 'interpreter', 'latex', 'FontSize', 8);

subplotXmanyY_er(7,3)
plot(estimatedState_history_w, 'color', cmap(3,:))
hold on
plot(X_real_w, 'color', cmap(4,:))
hold on
plot(Z_real_w, 'color', cmap(5,:))
legend('Estimate', 'True state', 'Measurement');
ylim([-10 10])
xlim([0 1050])
xlabel('Time Step', 'FontSize', 12);
ylabel('$Omega (rad/s)$', 'interpreter', 'latex', 'FontSize', 8);

subplotXmanyY_er(7,4)
for i = 1:numParticles
    ww = weights_history(:, i)';
    %loglog(weights_history_plt(:, i), 'color', cmap(1,:))
    plot(weights_history(:, i), 'color', cmap(1,:))
    hold on
end
xlabel('Time Step', 'FontSize', 12);
ylabel('Weights (I)', 'interpreter', 'latex', 'FontSize', 8);

subplotXmanyY_er(7,5)
for i = 1:numParticles
    ww = weights_history_w(:, i)';
    %loglog(weights_history_w_plt(:, i), 'color', cmap(1,:))
    plot(weights_history_w(:, i), 'color', cmap(1,:))
    hold on
end
xlabel('Time Step', 'FontSize', 12);
ylabel('Weights (omega)', 'interpreter', 'latex', 'FontSize', 8);

subplotXmanyY_er(7,6)
plot(N_eff_history, 'color', cmap(3,:))
xlabel('Time Step', 'FontSize', 12);
ylabel('Number of effective particles', 'interpreter', 'latex', 'FontSize', 8);

subplotXmanyY_er(7,7)
plot(N_eff_history_w, 'color', cmap(3,:))
xlabel('Time Step', 'FontSize', 12);
ylabel('Number of effective particles', 'interpreter', 'latex', 'FontSize', 8);

rmse = sqrt(mean((estimatedState_history - X_real).^2));
fprintf('RMSE of Estimation (RW Current): %f \n', rmse)
estimatedState_history_imputed = fillmissing(estimatedState_history, 'linear');  % Linearly interpolate missing values
rmse = sqrt(mean((estimatedState_history_imputed - X_real).^2));
fprintf('RMSE of Estimation (Imputed) (RW Current): %f \n', rmse)
rmse = sqrt(mean((Z_real - X_real).^2));
fprintf('RMSE of the measurement signal compared to the true state (RW Current): %f \n', rmse)


rmse = sqrt(mean((estimatedState_history_w - X_real_w).^2));
fprintf('RMSE of Estimation (RW Omega): %f \n', rmse)
estimatedState_history_imputed_w = fillmissing(estimatedState_history_w, 'linear');  % Linearly interpolate missing values
rmse = sqrt(mean((estimatedState_history_imputed_w - X_real_w).^2));
fprintf('RMSE of Estimation (Imputed) (RW Omega): %f \n', rmse)
rmse = sqrt(mean((Z_real_w - X_real_w).^2));
fprintf('RMSE of the measurement signal compared to the true state (RW Omega): %f \n', rmse)

fprintf('Elapsed time: %.2f minutes\n', elapsedTimeSeconds / 60);


%% Custom functions:

function measurement = getMeasurementAtTime(t, measurements)
% Check if t is within the range of measurements
if t >= 1 && t <= length(measurements)
    measurement = measurements(t);
else
    error('Time t is out of the range of available measurements');
end
end

function [x, m, w, p_h, p_h_w, m_h, w_w]=resample(x, m, w, N, p_h, p_h_w, m_h, w_w)
% Multinomial sampling with Ripley's method
u=cumprod(rand(1,N).^(1./[N:-1:1]));
u=fliplr(u);
wc=cumsum(w);
k=1;
for i=1:N
    while(wc(k)<u(i))
        k=k+1;
    end
    ind(i)=k;
end
x=x(:,ind);
m=m(:,ind);
p_h=p_h(:,ind);
p_h_w = p_h_w(:,ind);
m_h=m_h(:,ind);
w=ones(1,N)./N;
w_w=ones(1,N)./N;
end

% Presampling is only used for auxiliary particle filters
function [x, m, w, p_h, m_h, int_res]=presample(x, m, w, N, p_h, m_h, int_res)
% Multinomial sampling with Ripley's method
u=cumprod(rand(1,N).^(1./[N:-1:1]));
u=fliplr(u);
wc=cumsum(w);
k=1;
for i=1:N
    while(wc(k)<u(i))
        k=k+1;
    end
    ind(i)=k;
end
x=x(:,ind);
m=m(:,ind);
p_h=p_h(:,ind);
m_h=m_h(:,ind);
int_res = int_res(:,ind);
w=ones(1,N)./N;
end

function [x, w, p_h, int_res]=presample_w(x, w, N, p_h, int_res)
% Multinomial sampling with Ripley's method
u=cumprod(rand(1,N).^(1./[N:-1:1]));
u=fliplr(u);
wc=cumsum(w);
k=1;
for i=1:N
    while(wc(k)<u(i))
        k=k+1;
    end
    ind(i)=k;
end
x=x(:,ind);
p_h=p_h(:,ind);
int_res = int_res(:,ind);
w=ones(1,N)./N;
end
