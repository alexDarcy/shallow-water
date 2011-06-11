% Solves the Riemann problem on each axis 

function F  = riemann(isX,qL,qR)
  if (isX == 1)
    F  = riemannX(qL,qR);
  else
    F  = riemannY(qL,qR);
  end
end
    
function F  = riemannX(qL,qR)
  g = 9.81;
  hL = qL(1);
  uL = qL(2)/qL(1);
  vL = qL(3)/qL(1);
  hR = qR(1);
  uR = qR(2)/qR(1);
  vR = qR(3)/qR(1);

  F = zeros(3,1);

  hBar = 0.5*(hL+hR);
  uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
  cTilde = sqrt(g*hBar);

  r1 = zeros(3,1);
  r2 = zeros(3,1);
  r3 = zeros(3,1);
  alpha = zeros(3,1);
  lambda = zeros(3,1);
  r1(1) = 1;
  r1(2) = uTilde-cTilde;
  r1(3) = vTilde;
  r2(1) = 0;
  r2(2) = 0;
  r2(3) = 1;
  r3(1) = 1;
  r3(2) = uTilde+cTilde;
  r3(3) = vTilde;

  delta(1) = qR(1)-qL(1);
  delta(2) = qR(2)-qL(2);
  delta(3) = qR(3)-qL(3);
  alpha(1) = ((uTilde+cTilde)*delta(1)-delta(2))/(2*cTilde);
  alpha(2) = -vTilde*delta(1)+delta(3);
  alpha(3) = (-(uTilde-cTilde)*delta(1)+delta(2))/(2*cTilde);
  lambda(1) = uTilde-cTilde;
  lambda(2) = uTilde;
  lambda(3) = uTilde+cTilde;
  w = 0.5*(phi(lambda(1))*alpha(1)*r1 + phi(lambda(2))*alpha(2)*r2 +phi(lambda(3))*alpha(3)*r3);

  F(1) = 0.5*(qL(2)+qR(2));
  F(2) = 0.5*(qL(2)^2/qL(1)+0.5*g*qL(1)^2 + qR(2)^2/qR(1)+0.5*g*qR(1)^2);
  F(3) = 0.5*(qL(2)*qL(3)/qL(1)+qR(2)*qR(3)/qR(1)) ;
  F = F - w;

end

function G = riemannY(qL,qR)
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
  G = zeros(3,1);

  r1 = zeros(3,1);
  r2 = zeros(3,1);
  r3 = zeros(3,1);
  alpha = zeros(3,1);
  lambda = zeros(3,1);

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
  w = 0.5*(phi(lambda(1))*alpha(1)*r1 + phi(lambda(2))*alpha(2)*r2 +phi(lambda(3))*alpha(3)*r3);


  G(1) = 0.5*(qL(3)+qR(3));
  
  G(2) = 0.5*(qL(2)*qL(3)/qL(1)+qR(2)*qR(3)/qR(1)) ;

  G(3) = 0.5*(qL(3)^2/qL(1)+0.5*g*qL(1)^2 + qR(3)^2/qR(1)+0.5*g*qR(1)^2);

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

