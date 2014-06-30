h_values = [0.2, 0.1, 0.05, 0.02]
avg_errors = zeros(length(h_values))
for i = 1 : length(h_values)
    xs = 0 : h_values(i) : 1
    ys = 0 : h_values(i) : 1
    
    //Z = get_z(xs, ys, kiyono)
    //Z_cg = cg()
    Z = get_z(xs, ys, sample8)
    Z_cg = cg(all_four, h_values(i), 10^-10, 0, 1, 0, 1, sample7, sample6, sample5, sample4)
    
    //Eliminate bounderies
    Z_test = Z(2:(length(xs) - 1), 2:(length(ys) - 1))
    Z_cg_test = Z_cg((2:length(xs) - 1), 2:(length(ys) - 1))
    
    //disp(Z_test)
    
    error_matrix = map_matrix(Z_test - Z_cg_test, abs)
    avg_error = sum(error_matrix / (size(Z_test, 'c') * size(Z_test, 'r')))
    
    avg_errors(i) = avg_error
    disp(h_values(i), avg_error)    
end

plot2d(h_values, avg_errors)
