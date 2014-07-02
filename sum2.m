function second_sun = sum2(j,y_tilde,r_j_tilde)
global n mu 
second_sun = 0;
for k = 1:n
    if k == j 
         continue;
    else
         mu_k = mu(k);
         r_k_tilde = y_tilde(3*(k-1)+1:3*(k-1)+3,1);
         second_sun = second_sun + mu_k/norm(r_j_tilde-r_k_tilde);
     end
end   

end