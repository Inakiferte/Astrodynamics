function [c] = c_function (i,z)
% Function that computes Stumpff function ci € [0,1,2,3].
if i == 0  % C_0(z)

    if z > 0
        c = cos(sqrt(z));
    elseif z < 0
        c = cosh(sqrt(-z));
    else
        c = 1;
    end
end

if i == 1  % C_1(z)

    if z > 0
        c = sin(sqrt(z)) / sqrt(z);
    elseif z < 0 
        c = sinh(sqrt(-z)) / sqrt(-z);
    else
        c = 1;
    end
end

if i == 2  % C_2(z)

     c = (1-c_function(0,z))/z;
end

if i == 3  % C_3(z)

     c = (1-c_function(1,z))/z;
end

if i == 4 % C_4(z)

    c = (0.5 - c_function(2,z))/z;
end

if i == 5 % C_5(z)

    c = ((1/6) - c_function(3,z))/z;
end

end