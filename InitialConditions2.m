function [rho,u,p] = InitialConditions2(x)  %激波管问题
    if x<0.5
        rho = 1; 
        u = 0;
        p = 1;
    else
        rho = 0.125; 
        u = 0;
        p = 0.1;
    end
end 