function [X,F,Iters] = newton_opt(N, X, gradToler, XToler, MaxIter, myFx)

bGoOn = true;
Iters = 0;

while bGoOn

  Iters = Iters + 1;
  if Iters > MaxIter
    break;
  end

  g = FirstDerivatives(X, N, myFx);
  fnorm = norm(g);
  if fnorm < gradToler 
    break;
  end 
  J = SecondDerivatives(X, N, myFx);
  DeltaX = g / J;
  
  X = X - DeltaX;
  
  bStop = true;
  for i=1:N
    if abs(DeltaX(i)) > XToler(i)
      bStop = false;
    end
  end  
  
  bGoOn = ~bStop; 
  
end

F = feval(myFx, X);


function FirstDerivX = FirstDerivatives(X, N, myFx)

for iVar=1:N  
  xt = X(iVar);
  h = 0.01 * (1 + abs(xt));
  X(iVar) = xt + h;
  fp = feval(myFx, X);
  X(iVar) = xt - h;
  fm = feval(myFx, X);
  X(iVar) = xt;
  FirstDerivX(iVar) = (fp - fm) / 2 / h;    
end
end

function SecondDerivX = SecondDerivatives(X, N, myFx)

for i=1:N  
  for j=1:N
    if i == j
      f0 = feval(myFx, X);
      xt = X(i);
      hx = 0.01 * (1 + abs(xt));
      X(i) = xt + hx;
      fp = feval(myFx, X);
      X(i) = xt - hx;
      fm = feval(myFx, X);
      X(i) = xt;
      y = (fp - 2 * f0 + fm) / hx ^ 2;
    else
      xt = X(i);
      yt = X(j);
      hx = 0.01 * (1 + abs(xt));
      hy = 0.01 * (1 + abs(yt));

      X(i) = xt + hx;
      X(j) = yt + hy;
      fpp = feval(myFx, X);

      X(i) = xt - hx;
      X(j) = yt - hy;
      fmm = feval(myFx, X);

      X(i) = xt + hx;
      X(j) = yt - hy;
      fpm = feval(myFx, X);

      X(i) = xt - hx;
      X(j) = yt + hy;
      fmp = feval(myFx, X);
      X(i) = xt;
      X(j) = yt;
      y = (fpp - fmp - fpm + fmm) / (4 * hx * hy);
    end  
    SecondDerivX(i,j) = y;
  end
end
end

end
