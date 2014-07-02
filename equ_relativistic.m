%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the right hand side of the equations of motion
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the right hand side (RHS) of the equations of motion, i.e. the
% acceleration. The loop over i corresponds to looping over all bodies.
% The acceleration for each body is stored in a vector r_dot_dot.
% The right hand side contains the velocities and accelerations (derivative
% of the state vector (r,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RHS = equ_relativistic(t,y_tilde)
global n %fileID
% The equation to solve is M(t,y)y' = f(t,y) = RHS

for i = 1:n 
    [gamma_r_dot_dot_i] = rel_acceleration(i,y_tilde);
    r_dot_dot(3*(i-1)+1:3*(i-1)+3,1) = gamma_r_dot_dot_i;
end


RHS = [y_tilde(3*n+1:6*n);r_dot_dot];

%fprintf(fileID,'%.14f\n',t);

% Comment this out if you don't want to see the time. Otherwise, this
% allows you to know how far along the integration you are
t

end