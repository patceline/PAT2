%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes and plots the numerical error
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical error is computed (difference between numerical error and data
% used as input). The maximal error over the timespan is computed as well
% as the error in the last time step (for comparisons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [max_error_apo,last_error] = ...
          plot_error(position,actual_position)

% Error vector is the same size as the output
m = length(position);
% Initialize error vector
error = zeros(m,1);

% Compute the error. The sign is changed when 
% position(i)<actual_position(i) so that the oscillatory behavior can be
% seen (if there is one)
for i = 1:length(position)
    if position(i)>actual_position(i)
        error(i) = norm(position(i)-actual_position(i));
    elseif position(i)<actual_position(i)
        error(i) = -norm(position(i)-actual_position(i));
    end
end     

% Compute maximal error
max_error_apo = max(abs(error));
% Compute error in the last time step
last_error = abs(error(end));

% Plot the error (in norm)
figure;
hold on;
plot(error);
set(gca,'FontSize',14)
title('Numerical error in the position')
xlabel('Time [days]');
set(gca,'FontSize',14)
ylabel('Error [km]')

% Plot the error (in each component)
figure;
plot(position-actual_position);
set(gca,'FontSize',14)
title('Numerical error in the position')
legend('x','y','z')
xlabel('Time [days]');
set(gca,'FontSize',14)
ylabel('Error [km]')


end