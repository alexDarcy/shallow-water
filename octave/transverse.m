% Compute transverse flux
function [Zminus, Zplus] = transverse(isX, signA, qL, qR)
  g = 9.81;
  hL = qL(1);
  uL = qL(2)/qL(1);
  vL = qL(3)/qL(1);
  hR = qR(1);
  uR = qR(2)/qR(1);
  vR = qR(3)/qR(1);

  hBar = 0.5*(hL+hR);
  uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
  cTilde = sqrt(g*hBar);

  r1 = zeros(3,1);
  r2 = zeros(3,1);
  r3 = zeros(3,1);
  alpha = zeros(3,1);
  lambda = zeros(3,1);
  G = zeros(3,1);

  r1(1) = 1;
  r1(2) = uTilde;
  r1(3) = vTilde-cTilde;
  r2(1) = 0;
  r2(2) = -1;
  r2(3) = 0;
  r3(1) = 1;
  r3(2) = uTilde;
  r3(3) = vTilde+cTilde;

  delta(1) = qR(1)-qL(1);
  delta(2) = qR(2)-qL(2);
  delta(3) = qR(3)-qL(3);

  alpha(1) = ((vTilde+cTilde)*delta(1)-delta(3))/(2*cTilde);
  alpha(2) = uTilde*delta(1)-delta(2);
  alpha(3) = (-(vTilde-cTilde)*delta(1)+delta(3))/(2*cTilde);
  lambda(1) = vTilde-cTilde;
  lambda(2) = vTilde;
  lambda(3) = vTilde+cTilde;

  A = zeros(3,1);
  if (sign(lambda(1)) == signA)
    A = lambda(1)*alpha(1)*r1;
  end
  if (sign(lambda(2)) == signA)
    A = A + lambda(2)*alpha(2)*r2;
  end
  if (sign(lambda(3)) == signA)
    A = A + lambda(3)*alpha(3)*r3;
  end

  Zplus = zeros(3,1);
  Zminus = zeros(3,1);
  beta1 = computeBeta(isX,A,uTilde,vTilde,cTilde);
  if (lambda(1) < 0)
    Zminus = beta1(1)*lambda(1)*r1;
  else
    Zplus = beta1(1)*lambda(1)*r1;
  end

  if (lambda(2) < 0)
    Zminus = Zminus+beta1(2)*lambda(2)*r2;
  else
    Zplus = Zplus+beta1(2)*lambda(2)*r2;
  end

  if (lambda(3) < 0)
    Zminus = Zminus+beta1(3)*lambda(3)*r3;
  else
    Zplus = Zplus+beta1(3)*lambda(3)*r3;
  end
end

function coef = computeBeta(isX, A, u, v, c)
  coef = zeros(3,1);
  if (isX == 1)
    coef(1) = ((u+c)*A(1)-A(2))/(2*c);
    coef(2) = -v*A(1)+A(3);
    coef(3) = (-(u-c)*A(1)+A(2))/(2*c);
  else
    coef(1) = ((v+c)*A(1)-A(3))/(2*c);
    coef(2) = u*A(1)-A(2);
    coef(3) = (-(v-c)*A(1)+A(3))/(2*c);
  end
end


