%% Step 3: Generate the animation
% Initialize video writer
videoFileName = 'matching_pursuit_animation.mp4';  % Video file name
frameRate = 6;  % Frame rate
v = VideoWriter(videoFileName, 'MPEG-4');  % Create a video writer object
v.FrameRate = frameRate;  % Set frame rate
open(v);  % Open the video file for writing

p0 = cell(1, numberOfMeasurmentOfOriginalSignal);

for k = 1:numberOfMeasurmentOfOriginalSignal
    fig = figure('Visible', 'off', 'Color', 'w');
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
    fig = figure('Visible', 'off', 'Color', 'w');
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
        fig = figure('Visible', 'off', 'Color', 'w');
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
