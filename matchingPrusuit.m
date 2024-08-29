clc;
close all;
clear all;

dim = 200;
numberOfMeasurmentOfOriginalSignal = 30;
numberOfCoefficient = 6;
posisionOfCoffecients=  sort(randsample(40, numberOfCoefficient));
valueOfCoefficients = -0.8 + (1.6) .* rand(1, numberOfCoefficient);

% Create a sparse vector
idealX_CoefficientSparseArray = sparse(posisionOfCoffecients, 1, valueOfCoefficients, dim, 1);  % Create sparse vector

z_values = linspace(-7, 7, dim); %(1x200)
f_original = zeros(1, dim);
for i = 1:dim
    f_original(i) = compute_g(z_values(i),numberOfCoefficient,valueOfCoefficients,posisionOfCoffecients);  % Compute g(z) for each point in z_values
end

measuringPositionsOfOriginalSignal = sort(randperm(dim, numberOfMeasurmentOfOriginalSignal)); 
AtomSparseMatrix = sparse(1:numberOfMeasurmentOfOriginalSignal, measuringPositionsOfOriginalSignal, 1, numberOfMeasurmentOfOriginalSignal, dim); % 30x200

% Matrix of change of basis from Hermite to canonical
Phi_HermitePolynomialsMatrix = zeros(dim, dim);
for i = 1:dim
    Phi_HermitePolynomialsMatrix(i, :) = arrayfun(@(x) hermitFunc(i-1, x), z_values); % Basis matrix (Hermite polynomials)
end

A = AtomSparseMatrix * transpose(Phi_HermitePolynomialsMatrix);    % Measurement matrix (30x200)
f_hat = f_original(measuringPositionsOfOriginalSignal)';

%% Step 2: Orthogonal Matching Pursuit Algorithm

% Initialization
epsilon0 = 1e-10;    % Error threshold
gauss_waiting_factors = zeros(dim, 1);   % Initial guess (all coefficients zero) 200x1
residual = f_hat;         % Initial residual is the measurement vector 
indexOfComplementOfSupport = 1:dim;          % Complement of the support (all indices)

% Calculate the error for all columns of A that are not already in the support
epsilon = arrayfun(@(j) norm((A(:, j)' * residual / norm(A(:, j))^2)* A(:, j) - residual)^2, indexOfComplementOfSupport); %error

% Select the one with the smallest error
[~, min_idx] = min(epsilon);
j_MinValueErrorIndex = indexOfComplementOfSupport(min_idx);
% Update the error and complement of the support
epsilon = epsilon(j_MinValueErrorIndex);
indexOfComplementOfSupport(indexOfComplementOfSupport == j_MinValueErrorIndex) = []; % Update the complement of the support
support = setdiff(1:dim, indexOfComplementOfSupport); % Update the support

% Extract the columns of A corresponding to the support set
A_OfSupport = A(:, support);
% Find the best fit for the new estimate of x (least squares solution)
gauss_waiting_factors(support) = (A_OfSupport' * A_OfSupport) \ (A_OfSupport' * f_hat);
% Update the residuals
residual = f_hat - A * gauss_waiting_factors;
tmp=gauss_waiting_factors;
result_waiting_factors = {};      % To store results


% Orthogonal Matching Pursuit
while norm(residual) > epsilon0 % Continue until the residual is small enough
    % Calculate error vector epsilon
    epsilon = arrayfun(@(j) norm((A(:, j)' * residual / norm(A(:, j))^2) * A(:, j) - residual)^2, indexOfComplementOfSupport);
    % Select the index with the minimum error (best fit)
    [~, min_idx] = min(epsilon);
    % index_min_epsilon = nS(min_idx);     % Selected index
    epsilon = epsilon(min_idx); % Update epsilon

    % Update the complement and support sets
    indexOfComplementOfSupport(min_idx) = []; % Remove the selected index from the complement set
    support = setdiff(1:dim, indexOfComplementOfSupport);  % Update the support set

    % Update the solution for the new support set
    A_OfSupport = A(:, support);
    gauss_waiting_factors(support) = (A_OfSupport' * A_OfSupport) \ (A_OfSupport' * f_hat); % Least squares solution
    residual = f_hat - A * gauss_waiting_factors;    % Update the residuals
    result_waiting_factors{end+1} = gauss_waiting_factors;
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

p0 = cell(1, numberOfMeasurmentOfOriginalSignal);

for k = 1:numberOfMeasurmentOfOriginalSignal
    fig = figure('Visible', 'off');
    set(fig, 'Position', [100, 100, 1200, 600]);  % Adjust figure size: [left, bottom, width, height]
    % Create the first subplot for the line plot
    subplot(1, 2, 1);  % 1 row, 2 columns, first subplot
    plot(f_original, 'k', 'LineWidth', 2);
    hold on;
    plot(measuringPositionsOfOriginalSignal(1:k), f_hat(1:k), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');  % Overlay red points
    ylim([-1, 1]);
    legend('Original Function', 'Measurements');
    title('Line Plot');
    axis off;

     % Second plot (horizontal line with shaded area)
    subplot(1, 2, 2);  % Second subplot (right side)
    stem(1:50, zeros(1, 50), 'm', 'MarkerSize', 5, 'MarkerFaceColor', 'm');  % Horizontal line at y=0
    ylim([-1, 1]);
    xlim([0,50]);
    xlabel('Element of the basis', 'FontWeight', 'bold');
    ylabel('Coefficient', 'FontWeight', 'bold');
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontWeight', 'bold');  % Set label style
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
    fig = figure('Visible', 'off');
    set(fig, 'Position', [100, 100, 1200, 600]);  % Adjust figure size: [left, bottom, width, height]
    
    % First subplot (Left): Plot df and evolving signal
    subplot(1, 2, 1);  % Create first subplot
    plot(f_original, 'k', 'LineWidth', 1.5);  % Plot df as gray line
    hold on;
    evolving_signal = (1 - tau)*0 + tau*transpose(Phi_HermitePolynomialsMatrix)*result_waiting_factors{1};  % Evolving signal
    plot(evolving_signal, 'Color', [1, 0.5, 0] ,'LineWidth', 1.5);  % Plot evolving signal
    plot(measuringPositionsOfOriginalSignal, f_hat, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');  % Plot points overlay
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
        fig = figure('Visible', 'off');
        set(fig, 'Position', [100, 100, 1200, 600]);  % Adjust figure size: [left, bottom, width, height]
        % First subplot: Plot df and evolving signal
        subplot(1, 2, 1);  % Create the first subplot
        plot(f_original, 'k' ,'LineWidth', 1.5);  % Plot df as gray line
        hold on;
        evolving_signal_p2 = (1 - tau) * Phi_HermitePolynomialsMatrix' * result_waiting_factors{k - 1} + tau * Phi_HermitePolynomialsMatrix' * result_waiting_factors{k};
        plot(evolving_signal_p2 ,'Color', [1, 0.5, 0], 'LineWidth', 1.5);  % Plot evolving signal
        plot(measuringPositionsOfOriginalSignal, f_hat, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');  % Plot points overlay
        ylim([-1, 1]);  % Set y-axis limits
        axis off;  % Hide axes
        % Second subplot: Plot coefficients
        subplot(1, 2, 2);  % Create the second subplot
        coeffs_p2 = (1 - tau) * result_waiting_factors{k - 1} + tau * result_waiting_factors{k};
        h2 = stem(1:200, coeffs_p2, 'm', 'MarkerSize', 6, 'MarkerFaceColor', 'm');  % Stem plot of coefficients
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