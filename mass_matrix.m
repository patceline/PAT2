function M = mass_matrix(~,y)
global n
global gamma c mu
% make it sparse!!!
non_trivial_M = zeros(3*n,3*n);

for i = 1:n
    for j = 1:n
        if i == j
        non_trivial_M(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3) = eye(3);
        else
        r_i = y(3*(i-1)+1:3*(i-1)+3,1);
        r_j = y(3*(j-1)+1:3*(j-1)+3,1);

        %check this
        non_trivial_M(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3) = ...
            1/(2*c^2)/norm(r_j-r_i)^3*(r_j-r_i)*(r_j-r_i)'...
            -(3+4*gamma)/c^2*mu(j)/norm(r_j-r_i);
        end
    end
end

zero_matrix = zeros(3*n,3*n);
M = sparse([eye(3*n) zero_matrix; ...
     zero_matrix non_trivial_M]);




end
