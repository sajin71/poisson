function A = make_coefficient_matrix(h, num_of_x, num_of_y)
    //xs = xl : h : xr
    //ys = yl : h : yr
    
    //num_of_x = size(xs, 'c')
    //num_of_y = size(ys, 'c')
    
    A = zeros(num_of_x * num_of_y, num_of_x * num_of_y)
    for u_x_index = 1: num_of_x
        for u_y_index = 1: num_of_y
            row =  zeros(1, num_of_x * num_of_y)
            for i = 1 : num_of_x
                for j = 1: num_of_y
                    if i == u_x_index - 1 & j == u_y_index | i == u_x_index + 1 & j == u_y_index | i == u_x_index & j == u_y_index - 1 | i == u_x_index & j == u_y_index + 1 then
                        row(num_of_x * (i - 1) + j) = 1 / (h^2)
                     elseif ((i == u_x_index) & (j == u_y_index)) then
                         row(num_of_x * (i - 1) + j) = - 4 / (h^2)
                     end     
                end
            end
            A(((u_x_index - 1) * num_of_x + u_y_index),:) = row
        end
    end
endfunction


function b = make_b_vector(f, h, xl, xr, yl, yr, xl_f, xr_f, yl_f, yr_f)
    xs = xl : h : xr
    ys = yl : h : yr
    
    num_of_x = size(xs, 'c') - 2
    num_of_y = size(ys, 'c') - 2
    
    b = zeros(num_of_x * num_of_y, 1)
    //disp(xs, ys)
    for j = 1: num_of_y
        for i = 1: num_of_x
           b(i + num_of_y * (j - 1)) = b(i+ num_of_y * (j - 1)) - f(xs(i + 1), ys(j + 1))
           
           if i == 1 then
               //disp(i + num_of_y * (j- 1))
               b(i + num_of_y * (j - 1)) = b(i + num_of_y * (j - 1)) - xl_f(xs(1), ys(j + 1)) / (h ^ 2)
           end
           
           if j == 1 then
                b(i + num_of_y * (j - 1)) = b(i + num_of_y * (j - 1)) - yl_f(xs(i + 1), ys(1)) / (h ^ 2)
           end
           
           if i == num_of_x then
               b(i + num_of_y * (j - 1)) = b(i + num_of_y * (j - 1)) - xr_f(xs(num_of_x), ys(j + 1)) / (h ^ 2)
           end
           
           if j == num_of_y then
                b(i + num_of_y * (j - 1)) = b(i + num_of_y * (j - 1)) - yr_f(xs(i + 1), ys(num_of_y)) / (h ^ 2)
           end
           
        end
    end
endfunction


function Z = make_z_matrix(u, xs, ys, xl_f, xr_f, yl_f, yr_f)
    //disp(size(u))
    
    num_of_x = size(xs, 'c')
    num_of_y = size(ys, 'c')
    
    //disp(num_of_x, num_of_y)
    for i = 1: num_of_x
        for j = 1: num_of_y
            if i == 1 then
                Z(1, j) = xl_f(xs(1), ys(j))
            elseif j ==1 then
                Z(i, 1) = yl_f(xs(i), ys(1))
            elseif i == num_of_x then
                Z(num_of_x, j) = xr_f(xs(num_of_x), ys(j))
            elseif j == num_of_y then
                Z(i, num_of_y) = yr_f(xs(i), ys(num_of_y))
            else
                //disp(i, j)
                //disp((i - 1) + (num_of_y - 2) * (j - 1))
                Z(i, j) = u((i - 1) + (num_of_y - 2) * (j - 2))
            end  
        end 
    end
endfunction


function u = adaptability(f, xs, ys)
    u = zeros((length(xs)  - 2) * (length(ys) -2))
    for i = 1 : (length(ys) -2)
        for j = 1: (length(xs) -2)
            u(j + (i - 1)* (length(ys) -2)) = f(xs(j + 1), ys(i + 1))
        end
    end
endfunction

function Z = cg(f, h, ep, xl, xr, yl, yr, xl_f, xr_f, yl_f, yr_f)
    xs = xl : h : xr
    ys = yl : h : yr
    
    num_of_x = size(xs, 'c') - 2
    num_of_y = size(ys, 'c') - 2
    
    A = make_coefficient_matrix(h, num_of_x, num_of_y)
    b = make_b_vector(f, h, xl, xr, yl, yr, xl_f, xr_f, yl_f, yr_f)
    
    u = zeros(num_of_x * num_of_y, 1)
    next_u = zeros(num_of_x * num_of_y, 1)
    r = zeros(num_of_x * num_of_y, 1)
    next_r = zeros(num_of_x * num_of_y, 1)
    
    //adaptability test for kiyono
    correct_u = adaptability(kiyono, xs, ys)
    //disp(correct_u)
    //disp(length(xs), length(ys))
    //disp(size(A))
    //disp(size(correct_u))
    //disp(size(b))
    //disp(length(xs), length(ys))
    adaptability_error= map_matrix(A * correct_u - b, abs)
    //disp(adaptability_error)
    avg_error = sum(adaptability_error) / (size(adaptability_error, 'c') * size(adaptability_error, 'r'))
    disp(avg_error)
    //adaptablility test end
    
    p = b - A * u
    next_r = b - A * u
    i = 0 //counter
    alpha = 0
    beeta = 0
    while(norm(A * u - b) / norm(b) > ep)
        i  = i + 1
        //disp(p)
        //disp(r)
        u = next_u
        r = next_r
        alpha = (r' * r) / (p' * A * p)
        next_r = r - alpha * (A * p)
        next_u = u + alpha * p
        beeta = (next_r' * next_r) / (r' * r)
        p = next_r + beeta * p
        //disp(alpha)
     end
     
     Z = make_z_matrix(u, xs, ys, xl_f, xr_f, yl_f, yr_f)
     disp(i)
     //plot3d(xs, ys,  Z)
     //scf(14)
     //clf(14)
     //xset("colormap",jetcolormap(64))
     //surfxs, ys, Z)
 

endfunction
