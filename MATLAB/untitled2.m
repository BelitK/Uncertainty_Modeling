% Define the state transition function
function x_next = stateTransitionFcn(x, u)
    A = [1 0; 0 1]; % State transition matrix
    B = [1; 1];     % Control input matrix
    x_next = A * x + B * u; % Update state
end

% Define the measurement function
function z = measurementFcn(x)
    C = [1 0]; % Measurement matrix
    z = C * x; % Output measurement
end

% Initial state and covariance
initialState = [0; 0]; % Initial state vector
initialCovariance = eye(2); % Initial covariance matrix

% Create the extended Kalman filter object
ekf = extendedKalmanFilter(@stateTransitionFcn, @measurementFcn, initialState);

% Set process and measurement noise
ekf.ProcessNoise = 0.1 * eye(2); % Process noise covariance
ekf.MeasurementNoise = 0.1;       % Measurement noise covariance

% Simulate the EKF
numSteps = 100;
u = rand(numSteps, 1); % Random control inputs
measurements = zeros(numSteps, 1); % Placeholder for measurements
estimates = zeros(numSteps, 2); % Store estimates for plotting

for k = 1:numSteps
    % Simulate a measurement
    trueState = stateTransitionFcn(initialState, u(k));
    measurements(k) = measurementFcn(trueState) + randn() * sqrt(ekf.MeasurementNoise);

    % Predict the next state with control input
    predict(ekf, u(k)); % Pass the control input as an argument

    % Correct the state with the measurement
    correct(ekf, measurements(k));

    % Store the estimated state
    estimates(k, :) = ekf.State';

    % Update the initial state for the next iteration
    initialState = ekf.State;
end

% Plot the results
figure;
subplot(2, 1, 1);
plot(1:numSteps, measurements, 'r', 'DisplayName', 'Measurements');
hold on;
plot(1:numSteps, estimates(:, 1), 'b', 'DisplayName', 'Estimated State 1');
plot(1:numSteps, estimates(:, 2), 'g', 'DisplayName', 'Estimated State 2');
xlabel('Time Step');
ylabel('State Value');
title('State Estimation');
legend;
grid on;

subplot(2, 1, 2);
plot(1:numSteps, measurements - estimates(:, 1), 'k', 'DisplayName', 'Measurement Residuals');
xlabel('Time Step');
ylabel('Residual');
title('Measurement Residuals');
legend;
grid on;

% Display the final estimated state
disp('Estimated State:');
disp(ekf.State);
