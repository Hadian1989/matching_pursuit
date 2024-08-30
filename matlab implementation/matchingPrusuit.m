clc;
close all;
clear all;
%% Step 1: Initialization
dim = 200;
numOfMeasurmentOfOrigSignal = 30;
numOfCoefficient = 6;
posisionOfCoffs=  sort(randsample(40, numOfCoefficient));
valueOfCoeffs = -0.8 + (1.6) .* rand(1, numOfCoefficient);

% Create a sparse vector
idealX_CoeffSparseArray = sparse(posisionOfCoffs, 1, valueOfCoeffs, dim, 1);

z_values = linspace(-7, 7, dim); %(1x200)
f_original = zeros(1, dim);
for i = 1:dim
    % Compute g(z) for each point in z_values
    f_original(i) = compute_g(z_values(i),numOfCoefficient,valueOfCoeffs,posisionOfCoffs);
end

measuringPositionsOfOrigSignal = sort(randperm(dim, numOfMeasurmentOfOrigSignal)); 
%  30x200
AtomSparseMatrix = sparse(1:numOfMeasurmentOfOrigSignal, measuringPositionsOfOrigSignal, 1, numOfMeasurmentOfOrigSignal, dim);

% Matrix of change of basis from Hermite to canonical
Phi_HermitePolyMatrix = zeros(dim, dim);
for i = 1:dim
    % Basis matrix (Hermite polynomials)
    Phi_HermitePolyMatrix(i, :) = arrayfun(@(x) hermitFunc(i-1, x), z_values); 
end
% Measurement matrix (30x200)
A_ProjectionMatrix = AtomSparseMatrix * transpose(Phi_HermitePolyMatrix);
f_hat = f_original(measuringPositionsOfOrigSignal)';

%% Step 2: Orthogonal Matching Pursuit Algorithm

% Initialization
% Error threshold
epsilon0 = 1e-10;
% Initial guess (all coefficients zero) 200x1
gaussed_waiting_factors = zeros(dim, 1);
% Initial residual is the measurement vector
residual = f_hat; 
% Complement of the support (all indices)
indexOfCompOfSupport = 1:dim;
% Calculate the error for all columns of A that are not already in the support
epsilon = arrayfun(@(j) norm((A_ProjectionMatrix(:, j)' * residual / norm(A_ProjectionMatrix(:, j))^2)* A_ProjectionMatrix(:, j) - residual)^2, indexOfCompOfSupport);
% Select the one with the smallest error
[~, min_idx] = min(epsilon);
j_MinValueErrorIndex = indexOfCompOfSupport(min_idx);
% Update the error and complement of the support
epsilon = epsilon(j_MinValueErrorIndex);
% Update the complement of the support
indexOfCompOfSupport(indexOfCompOfSupport == j_MinValueErrorIndex) = [];
% Update the support
support = setdiff(1:dim, indexOfCompOfSupport); 
% Extract the columns of A corresponding to the support set
A_OfSupport = A_ProjectionMatrix(:, support);
% Find the best fit for the new estimate of x (least squares solution)
gaussed_waiting_factors(support) = (A_OfSupport' * A_OfSupport) \ (A_OfSupport' * f_hat);
% Update the residuals
residual = f_hat - A_ProjectionMatrix * gaussed_waiting_factors;
tmp = gaussed_waiting_factors;
% To store results
result_waiting_factors = {};


% Orthogonal Matching Pursuit
while norm(residual) > epsilon0 
    % Calculate error vector epsilon
    epsilon = arrayfun(@(j) norm((A_ProjectionMatrix(:, j)' * residual / norm(A_ProjectionMatrix(:, j))^2) * A_ProjectionMatrix(:, j) - residual)^2, indexOfCompOfSupport);
    % Select the index with the minimum error (best fit)
    [~, min_idx] = min(epsilon);
    % Update epsilon
    epsilon = epsilon(min_idx); 
    % Remove the selected index from the complement set
    indexOfCompOfSupport(min_idx) = []; 
    % Update the support set
    support = setdiff(1:dim, indexOfCompOfSupport);
    % Update the solution for the new support set
    A_OfSupport = A_ProjectionMatrix(:, support);
    % Least squares solution
    gaussed_waiting_factors(support) = (A_OfSupport' * A_OfSupport) \ (A_OfSupport' * f_hat);
    % Update the residuals
    residual = f_hat - A_ProjectionMatrix * gaussed_waiting_factors;    
    result_waiting_factors{end+1} = gaussed_waiting_factors;
end
result_waiting_factors=[tmp,result_waiting_factors];

%% Step 3: Generate the animation
% % Define video dimensions

% Initialize video writer
videoFileName = 'matching_pursuit_animation.mp4';  % Video file name
frameRate = 6;  % Frame rate
v = VideoWriter(videoFileName, 'MPEG-4');  % Create a video writer object
v.FrameRate = frameRate;  % Set frame rate
open(v);  % Open the video file for writing

p0 = cell(1, numOfMeasurmentOfOrigSignal);

for k = 1:numOfMeasurmentOfOrigSignal
    fig = figure('Visible', 'off', 'Color', 'w');
     % Adjust figure size: [left, bottom, width, height]
    set(fig, 'Position', [100, 100, 1200, 600]); 
    % Create the first subplot for the line plot
    subplot(1, 2, 1); 
    plot(f_original, 'k', 'LineWidth', 2);
    hold on;
    plot(measuringPositionsOfOrigSignal(1:k), f_hat(1:k), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    ylim([-1, 1]);
    legend('Original Function', 'Measurements');
    title('Line Plot');
    axis off;
     % Second plot (horizontal line with shaded area)
    subplot(1, 2, 2); 
    stem(1:50, zeros(1, 50), 'm', 'MarkerSize', 5, 'MarkerFaceColor', 'm'); 
    ylim([-1, 1]);
    xlim([0,50]);
    xlabel('Element of the basis', 'FontWeight', 'bold');
    ylabel('Coefficient', 'FontWeight', 'bold');
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontWeight', 'bold'); 
    p0{k} = gcf;  % Save the current figure handle into the cell array
    % Capture the current frame and write it to the video
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Close the figure to free up memory
    close(gcf);
end
tau_values = 0:0.1:1;  % Tau values for the second loop

% Second Loop: Create the second set of figures and add them to the video
for tau = tau_values
    % Create a new figure for each frame
    fig = figure('Visible', 'off', 'Color', 'w');
    set(fig, 'Position', [100, 100, 1200, 600]);
    subplot(1, 2, 1);
    plot(f_original, 'k', 'LineWidth', 1.5);
    hold on;
    evolving_signal = (1 - tau)*0 + tau*transpose(Phi_HermitePolyMatrix)*result_waiting_factors{1};
    plot(evolving_signal, 'Color', [1, 0.5, 0] ,'LineWidth', 1.5);  % Plot evolving signal
    plot(measuringPositionsOfOrigSignal, f_hat, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    ylim([-1, 1]);  % Set y-axis limits
    axis off;  % Hide axes
    
    % Second subplot (Right): Plot the evolution of coefficients
    subplot(1, 2, 2);  % Create second subplot
    coeffs = tau * result_waiting_factors{1};  % Compute the evolving coefficients
    h = stem(1:200, coeffs, 'm', 'MarkerSize', 6, 'MarkerFaceColor', 'm');  % Plot coefficients
    % Ensure the color is consistent
    set(h, 'Color', 'm');  % Set color for the stem lines
    set(h, 'MarkerEdgeColor', 'm');  % Set color for marker edges
    set(h, 'MarkerFaceColor', 'm');  % Set color for marker faces
    ylim([-1, 1]);  % Set y-axis limits
    xlim([0,50]);
    xlabel('Element of the basis', 'FontWeight', 'bold');
    ylabel('Coefficient', 'FontWeight', 'bold');
    title('Coefficients');  % Add title
    
    % Capture the current frame and write it to the video
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Close the figure to free up memory
    close(gcf);
end

for k = 2:length(support) - 1
    for tau = tau_values
        % Create a new figure for each frame
        fig = figure('Visible', 'off', 'Color', 'w');
        set(fig, 'Position', [100, 100, 1200, 600]);  % Adjust figure size: [left, bottom, width, height]
        % First subplot: Plot df and evolving signal
        subplot(1, 2, 1);  % Create the first subplot
        plot(f_original, 'k' ,'LineWidth', 1.5);  % Plot df as gray line
        hold on;
        evolving_signal_p2 = (1 - tau) * Phi_HermitePolyMatrix' * result_waiting_factors{k - 1} ...
                           + tau * Phi_HermitePolyMatrix' * result_waiting_factors{k};
        plot(evolving_signal_p2 ,'Color', [1, 0.5, 0], 'LineWidth', 1.5);  % Plot evolving signal
        plot(measuringPositionsOfOrigSignal, f_hat, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        ylim([-1, 1]);  % Set y-axis limits
        axis off;  % Hide axes
        % Second subplot: Plot coefficients
        subplot(1, 2, 2);  % Create the second subplot
        coeffs_p2 = (1 - tau) * result_waiting_factors{k - 1} + tau * result_waiting_factors{k};
        h2 = stem(1:200, coeffs_p2, 'm', 'MarkerSize', 6, 'MarkerFaceColor', 'm');
        % Ensure the color is consistent
        set(h2, 'Color', 'm');  % Set color for the stem lines
        set(h2, 'MarkerEdgeColor', 'm');  % Set color for marker edges
        set(h2, 'MarkerFaceColor', 'm');  % Set color for marker faces
        ylim([-1, 1]);  % Set y-axis limits
        xlim([0,50]);
        xlabel('Element of the basis', 'FontWeight', 'bold');
        ylabel('Coefficient', 'FontWeight', 'bold');
        title('Coefficients');  % Add title
        
        % Capture the current frame and write it to the video
        frame = getframe(gcf);
        writeVideo(v, frame);
        
        % Close the figure to free up memory
        close(gcf);
    end
end

% Finalize and close the video file
close(v);

function g_val = compute_g(z,nb,c,nzc)
g_val = 0;
for j = 1:nb
    g_val = g_val + c(j) * hermitFunc(nzc(j), z);  % Sum of terms
end
end