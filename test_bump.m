function f = test_bump(r)
    a = 1/2;
    
    if(r > a)
        f = 0;
    else
        f = exp(-1/(a^2 - r^2));
    end
end