% Solves the Riemann problem on each axis 

function F = riemann(isX,i, j, q1, q2, q3)
  if (isX == 1)
    F = riemannX(i, j, q1, q2, q3);
  else
    F = riemannY(i, j, q1, q2, q3);
  end
end

function F = riemannX(i, j, q1, q2, q3)
  g = 9.81;
  hL = q1(i,j);
  uL = q2(i,j)/q1(i,j);
  vL = q3(i,j)/q1(i,j);
  hR = q1(i+1,j);
  uR = q2(i+1,j)/q1(i+1,j);
  vR = q3(i+1,j)/q1(i+1,j);
  hBar = 0.5*(hL+hR);
  uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  %if (abs(uR-uL) < 1e-9)
  vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
  %else
  %  aL = hL*(uTilde-uL);
  %  aR = hR*(uR-uTilde);
  %  vTilde = (aL*vL+aR*vR)/(aL+aR);
  %end
  cTilde = sqrt(g*hBar);

  r1(1) = 1;
  r1(2) = uTilde-cTilde;
  r1(3) = vTilde;
  r2(1) = 0;
  r2(2) = 0;
  r2(3) = 1;
  r3(1) = 1;
  r3(2) = uTilde+cTilde;
  r3(3) = vTilde;

  delta(1) = q1(i+1,j)-q1(i,j);
  delta(2) = q2(i+1,j)-q2(i,j);
  delta(3) = q3(i+1,j)-q3(i,j);
  alpha(1) = ((uTilde+cTilde)*delta(1)-delta(2))/(2*cTilde);
  alpha(2) = -vTilde*delta(1)+delta(3);
  alpha(3) = (-(uTilde-cTilde)*delta(1)+delta(2))/(2*cTilde);
  lambda(1) = uTilde-cTilde;
  lambda(2) = uTilde;
  lambda(3) = uTilde+cTilde;
  w = 0.5*(phi(lambda(1))*alpha(1)*r1 + phi(lambda(2))*alpha(2)*r2 +phi(lambda(3))*alpha(3)*r3);

  F(1) = 0.5*(q2(i,j)+q2(i+1,j));
  F(2) = 0.5*(q2(i)^2/q1(i)+0.5*g*q1(i)^2 + q2(i+1)^2/q1(i+1)+0.5*g*q1(i+1)^2);
  F(3) = 0.5*(q2(i,j)*q3(i,j)/q1(i,j)+q2(i+1,j)*q3(i+1,j)/q1(i+1,j)) ;
  F = F - w;

end

function G = riemannY(i, j, q1, q2, q3)
  g = 9.81;
  hL = q1(i,j);
  uL = q2(i,j)/q1(i,j);
  vL = q3(i,j)/q1(i,j);
  hR = q1(i,j+1);
  uR = q2(i,j+1)/q1(i,j+1);
  vR = q3(i,j+1)/q1(i,j+1);
  hBar = 0.5*(hL+hR);
  uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  %if (abs(uR-uL) < 1e-9)
    vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
  %else
  %  aL = hL*(uTilde-uL);
  %  aR = hR*(uR-uTilde);
  %  vTilde = (aL*vL+aR*vR)/(aL+aR);
  %end
 
  cTilde = sqrt(g*hBar);


  r1(1) = 1;
  r1(2) = uTilde;
  r1(3) = vTilde-cTilde;
  r2(1) = 0;
  r2(2) = -1;
  r2(3) = 0;
  r3(1) = 1;
  r3(2) = uTilde;
  r3(3) = vTilde+cTilde;

  delta(1) = q1(i,j+1)-q1(i,j);
  delta(2) = q2(i,j+1)-q2(i,j);
  delta(3) = q3(i,j+1)-q3(i,j);
  alpha(1) = ((vTilde+cTilde)*delta(1)-delta(3))/(2*cTilde);
  alpha(2) = uTilde*delta(1)-delta(2);
  alpha(3) = (-(vTilde-cTilde)*delta(1)+delta(3))/(2*cTilde);
  lambda(1) = vTilde-cTilde;
  lambda(2) = vTilde;
  lambda(3) = vTilde+cTilde;
  w = 0.5*(phi(lambda(1))*alpha(1)*r1 + phi(lambda(2))*alpha(2)*r2 +phi(lambda(3))*alpha(3)*r3);


  G(1) = 0.5*(q3(i,j)+q3(i,j+1));
  
  G(2) = 0.5*(q2(i,j)*q3(i,j)/q1(i,j)+q2(i,j+1)*q3(i,j+1)/q1(i,j+1)) ;

  G(3) = 0.5*(q3(i,j)^2/q1(i,j)+0.5*g*q1(i,j)^2 + q3(i,j+1)^2/q1(i,j+1)+0.5*g*q1(i,j+1)^2);

  G = G - w;
end

% Harten entropy fix
function z = phi(lambda)
  % empirical value
  epsilon = 2;
  if (abs(lambda) >= epsilon)
    z = abs(lambda);
  else
    z = (lambda^2 + epsilon^2)/(2*epsilon);
 end
end

