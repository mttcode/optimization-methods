function y = f2(x)
    global f2Calls;
    f2Calls = f2Calls + 1;
    %y = (x(1)-1).^2 + 10.*(x(1)-x(2).^2).^2;
    %y = (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2)^2 - 7)^2;
    y = (x(1)-1)^2 + 100*(x(2)-x(1)^2)^2;
end
