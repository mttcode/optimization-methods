function [X,LastF,Iters]=conjgrad(N,X,Eps_Fx,MaxIter,myFx)

initStep = 0.1;
minStep = 0.000001;
LastF = feval(myFx, X);
[dfnorm,Deriv] = getgradients(X, N, myFx);

lambda = 0;
lambda = linsearch(X, N, lambda, initStep, minStep, Deriv, myFx);
X = X + lambda * Deriv;

bGoOn = true;
Iters = 0;

while bGoOn

  Iters = Iters + 1;
  if Iters > MaxIter
    break;
  end

  dfnormOld = dfnorm;
  DerivOld = Deriv;
  [dfnorm,Deriv] = getgradients(X, N, myFx);
  Deriv = (dfnorm / dfnormOld)^2 * DerivOld - Deriv;
  if dfnorm < Eps_Fx
    break;
  end
  lambda = 0;
  lambda = linsearch(X, N, lambda, initStep, minStep, Deriv, myFx);
  X = X + lambda * Deriv;
  F = feval(myFx, X);
  if abs(F - LastF) < Eps_Fx
    bGoOn = false;
  else
    LastF = F;
  end

end

LastF = feval(myFx, X);

function y = myFxEx(N, X, DeltaX, lambda, myFx)

  X = X + lambda * DeltaX;
  y = feval(myFx, X);

function [fnorm,Deriv] = getgradients(X, N, myFx)

  for i=1:N
    xx = X(i);
    h = 0.01 * (1 + abs(xx));
    X(i) = xx + h;
    Fp = feval(myFx, X);
    X(i) = xx - h;
    Fm = feval(myFx, X);
    X(i) = xx;
    Deriv(i) = (Fp - Fm) / 2 / h;
  end
  fnorm = norm(Deriv);

function lambda = linsearch(X, N, lambda, initStep, minStep, D, myFx)

  f1 = myFxEx(N, X, D, lambda, myFx);
  while initStep > minStep
    f2 = myFxEx(N, X, D, lambda + initStep, myFx)  ;
    if f2 < f1
      f1 = f2;
      lambda = lambda + initStep;
    else
      f2 = myFxEx(N, X, D, lambda - initStep, myFx);
      if f2 < f1
        f1 = f2;
        lambda = lambda - initStep;
      else
        initStep = initStep / 10;
      end
    end
  end
