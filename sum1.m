function first_sum = sum1(i,y_tilde,r_i_tilde)
global n mu

first_sum = 0;

for k = 1:n
    if k == i 
        continue;
     else
        r_k_tilde = y_tilde(3*(k-1)+1:3*(k-1)+3,1);
        mu_k = mu(k);
        first_sum = first_sum + mu_k/norm(r_i_tilde-r_k_tilde);
     end
end


end