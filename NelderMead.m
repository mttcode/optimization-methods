clc;
clear;

disp('Zadaj nasledujuce hodnoty: ');
preklop = 1;
predlz  = 2;
skrat   = 0.5;
reduk   = 0.5;

x0      = input('Startovaci bod (zvycajne [-1, 1]): ');

e = 1e-50;
n = length(x0);

%Osetri ak zadany vektor x0 je vertikalny
if length(x0(1,:)) ~= length(x0)
    x0=x0';
end

% vytvorenie simplexu
x = preklop*eye(n+1, n);
for i=1:n+1
    x(i, :) = x(i, :) + x0;
end;
fData = zeros(n+1, 1);

global f2Calls;
f2Calls = 0;
iter = 0;
minGx = Inf;
maxGx = -Inf;
minGy = Inf;
maxGy = -Inf;
while (iter < 1000)
    iter = iter + 1;
    
    minGx = min([minGx min(x(:, 1))]);
    maxGx = max([maxGx max(x(:, 1))]);
    minGy = min([minGy min(x(:, 2))]);
    maxGy = max([maxGy max(x(:, 2))]);
    
    % vypocet hodnot fcie    
    for i=1:n+1
        fData(i) = f2(x(i, :));
    end;

    % usporiadame body
    for i=1:n+1
        for j=1:n
            if fData(j+1) < fData(j)
                temp       = fData(j);
                fData(j)   = fData(j+1);
                fData(j+1) = temp;
            
                temp      = x(j, :);
                x(j, :)   = x(j+1, :);
                x(j+1, :) = temp;            
            end;
        end;
    end;

    % vypocet a vyznacenie taziska
    x0  = sum(x(1:n, :)) ./ n;
    x0z = f2(x0);
  
    % dosiahli sme urcenu presnost ? (vsetky body su max. e daleko od taziska)
    quit = 1;
    for i=1:n+1
        if norm(x0 - x(i, :)) > e
            quit = 0;
        end;
    end;
    if quit > .5
        break;
    end;
    
    % reflekcia
    xr  = x0 + preklop.*(x0 - x(n+1, :));
    xrf = f2(xr);
    if (fData(1) <= xrf) && (xrf < fData(n))
        fData(n+1) = xrf;
        x(n+1, :)  = xr;
        continue;
    end;
    
    % expanzia
    if xrf < fData(1)
        xe  = x0 + predlz.*(x0 - x(n+1, :));
        xef = f2(xe);
        
        if xef < xrf
            x(n+1, :)  = xe;
            fData(n+1) = xef;            
        else
            x(n+1, :)  = xr;
            fData(n+1) = xrf;            
        end;
        continue;
    end;
    
    % kontrakcia
    if xrf >= fData(n)
        xc  = x(n+1, :) + skrat.*(x0 - x(n+1, :));
        xcf = f2(xc);
        
        if xcf < fData(n+1)
            x(n+1, :)  = xc;
            fData(n+1) = xcf; 
            continue;
        end;        
    end;
        
    % redukcia
    for i=2:n+1
        x(i, :) = x(1, :) + reduk.*(x(i, :) - x(1, :));
    end;      
end;

% Vypis
format long;
fprintf('\nPozadovane minimum sa naslo v bode [%d, %d] s hodnotou %d.\n', x0(1), x0(2), f2(x0))
fprintf('Pocet iteracii bol %d a pocet volani funkcie %d.\n', iter, f2Calls)
