%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the acceleration for body i
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the acceleration of body i given the relativistic equations of
% motion. The vector y_tilde containing the positions and velocities of 
% the other bodies is needed here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_dot_dot_i] = rel_acceleration(i,y_tilde)
global mu n beta gamma c 

% Compute the constants used in the equations of motion
% These are all terms in 1/c^2
[cst1,cst2,cst3,cst4,cst5,cst6,~,cst8,...
 cst9,cst10,~] = get_constants(beta,gamma,c);

% Initialize the acceleraiton of body i (3x1 vector)
r_dot_dot_i = zeros(3,1);

% Position and velocity of body i, extracted from the state vector y_tilde
r_i_tilde = y_tilde(3*(i-1)+1:3*(i-1)+3,1);
v_i_tilde = y_tilde(3*n+3*(i-1)+1:3*n+3*(i-1)+3,1);

% Computes the first sum in the equations of motion
first_sum = sum1(i,y_tilde,r_i_tilde);

% Loop over all bodies. Compute the contribution that each body has on body
% to the acceleration of body i
for j = 1:n
    if j == i
        r_dot_dot_i = r_dot_dot_i + 0;
    else
        
    % Position and velocity of body j, extracted from the state vector
    r_j_tilde = y_tilde(3*(j-1)+1:3*(j-1)+3,1);
    v_j_tilde = y_tilde(3*n+3*(j-1)+1:3*n+3*(j-1)+3,1);
    % mu of the body j
    mu_j = mu(j);
    
    % Compute the second sum in the equations of motion
    second_sum = sum2(j,y_tilde,r_j_tilde);
    
    beta_coef = coef(i);
    
    % Contribution of body j to the acceleration of body i
    r_dot_dot_i = r_dot_dot_i + ...
         mu_j*(r_j_tilde-r_i_tilde)/...
         norm(r_i_tilde-r_j_tilde)^3*...
         (...
         1-beta_coef-...
         cst1*first_sum-...
         cst2*second_sum+...
         cst3*dot(v_i_tilde,v_i_tilde)+...
         cst4*dot(v_j_tilde,v_j_tilde)-...
         cst5*dot(v_i_tilde,v_j_tilde)-...
         cst6*(dot(r_i_tilde-r_j_tilde,v_j_tilde)...
         /norm(r_i_tilde-r_j_tilde))^2)+...
         cst8*mu_j/norm(r_i_tilde-r_j_tilde)^3*(...
         dot(r_i_tilde-r_j_tilde,cst9*v_i_tilde-cst10*v_j_tilde)*...
         (v_i_tilde-v_j_tilde)...
         );
    
    end
end

end