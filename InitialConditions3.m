function [rho,u,p] = InitialConditions3(x)  %膨胀激波问题
    if x<0.5
        rho = 5; 
        u = sqrt(1.4);
        p = 29;
    else
        rho = 1; 
        u = 5*sqrt(1.4);
        p = 1;
    end
end 