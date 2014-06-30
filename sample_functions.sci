N = 100

function z = kiyono_imp(x, y)
    if x == 0 & y == 0 then
        z = 0
    else
        z = (x * y) / (x^2 + y^2)
    end
endfunction

function z = kiyono_imp_f(x, y)
    if (x == 0 & y == 0) then
        z = 0
    else
        z = - 12 * x * y / (x^2 + y^2)^2 + 8 * x * y^3 / (x^2 + y^2)^3 + 8 * x^3 * y / (x^2 + y^2)^3
        //z = 2 * x^2 / (x^2 + y^2)
        //z = z - 20 * x^2 * y^2 / (x^2 + y^2)^2
        //z = z + 2 * y^2 / (x^2 + y^2)
        //z = z + 8 * x^2 * y^4 / (x^2 + y^2)^3
        //z = z + 8 * x^4 * y^2 / (x^2 + y^2)^3
    end
endfunction

function z = kiyono_imp_x1_bound(x, y)
    z = (1 * y) / (1 + y^2) 
endfunction

function z = kiyono_imp_xminus1_bound(x, y)
    z = -1 * kiyono_imp_x1_bound(x, y)
endfunction

function z = kiyono_imp_y1_bound(x, y)
    z = kiyono_imp_x1_bound(y, x)
endfunction

function z = kiyono_imp_yminus1_bound(x, y)
    z = -1 * kiyono_imp_y1_bound(x, y)
endfunction

function z = kiyono(x, y)
    if x == 0 & y == 0 then
        z = 0
    else
        z = ((x ^ 3 )* y) / (x^2 + y^2)
    end
endfunction

function z = kiyono_x1_bound(x, y)
    z = (1 * y) / (1 + y^2)
endfunction

function z = kiyono_xminus1_bound(x, y)
    z = (-1 * y) / (1 + y^2)
endfunction

function z = kiyono_y1_bound(x, y)
    z = (1 * x^3) / (1 + x^2)
endfunction

function z = kiyono_yminus1_bound(x, y)
    z = (-1 * x^3) / (1 + x^2)
endfunction

function z = kiyono_f(x, y)
    if x== 0 & y == 0 then
        z = 0
    else
        z = 6 * x * y / (x^2 + y^2)
        z = z  + 8 * (x^5) * y / (x^2 + y^2)^3
        z = z  - 20 * (x^3) * y / (x^2 + y^2)^2
        z = z + 8 * (x^3) * (y^3) / (x^2 + y^2)^3
        
        //flip sign
        z = -z
    end
endfunction

function z = sample1(x, y)
    if y <= 0.5 then
        z = 1
    else
        z = -1
    end
endfunction

function z = sample4(x, y)
    z = 1 + x^2
endfunction

function z = sample5(x, y)
    z = x^2
endfunction

function z = sample6(x, y)
    z = 1 + y^2
endfunction

function z = sample7(x, y)
    z = y^2
endfunction

function z = sample8(x, y)
    z = x^2 + y^2
endfunction

function z = all_one(x, y)
    z = 1
endfunction

function z = all_four(x, y)
    z = -4
endfunction

function z = sample3(x, y)
    if x <= 0.5  then
        z = x
    else
        z = x - 1
    end
endfunction

function z = all_zero(x, y)
    z = 0
endfunction

function z = sample2(x, y)
    z = min(x, 1- x)
endfunction

function z = sin_pi(x, y)
    z = sin(%pi * x)
endfunction

function z = laplace(x, y)
    z = 0
    for m = 1 : N
        mth = (2 / %pi) * (floor(1 - (-1)^m)) / (m * sinh(m * %pi)) * sinh(m * %pi * y) * sin(m * %pi * x)
        z = z + mth
    end 
endfunction

function Z = get_laplace_strict(h)
    xs = 0 : h : 1
    ys = 0 : h : 1
    
    for i = 1 : length(xs)
        for j = 1 : length(ys)
            Z(i, j) = laplace(xs(i), ys(j))
        end
    end
endfunction
function z = poisson(x, y)
    //g = 1
    z = 0
    for m = 1 : N
        mth = 0//(4 / (%pi)^2) * (2 / ((2 * m + 1) %pi)) * (2 / ((2 * m ))) * sin((2 * m * 1) * %pi * x) * cos((2 * m + 1) * %pi * y)
        z = z + mth
    end
endfunction

function Z = get_z(xs, ys, f)
    for i = 1 : length(xs)
        for j = 1 : length(ys)
            Z(i, j) = f(xs(i), ys(j))
        end
    end
endfunction

function Z = map_matrix(X, f)
    for i = 1 : size(X, 'r')
        for j = 1: size(X, 'c')
            Z(i, j) = f(X(i, j))
        end
    end
endfunction
