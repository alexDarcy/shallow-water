function shallow2D
  Lx = 10;
  Ly = 10;
  n1 = 25;
  n2 = 25;
  dt = 0.01;
  g = 9.81;
  % Initial disk
  xC = Lx/2;
  yC = Ly/2;
  radius = 0.1;
  stepX = Lx/n1;
  stepY = Ly/n2;
  h1 = 1;
  h0 = 0.5;

  % Create mesh
  x = 0:stepX:Lx;
  y = 0:stepY:Ly;
  %z = zeros(nX+1,nY+1)
  nX = length(x);
  nY = length(y);

  for i=1:nX
    for j=1:nY
      %dist = (x(i)-xC)^2 + (y(j)-yC)^2;
      %dist = sqrt(dist);
      %if (dist <= radius )
      if (x(i) < Lx/2)
        h(i,j) = h1;
      else
        h(i,j) = h0;
      end
      u(i,j) = 0;
      v(i,j) = 0;
    end
  end

  surf(x,y,h);
  view(110,10);
  drawnow;

  q1 = h;
  for i=1:nX
    for j=1:nY
      q2(i,j) = u(i,j)*h(i,j);
      q3(i,j) = v(i,j)*h(i,j);
    end
  end

  for t=0:dt:150*dt
    for i=1:nX-1
      %for j=1:nY-1
      j = 1;

        hL = q1(i,j);
        uL = q2(i,j)/q1(i,j);
        vL = q3(i,j)/q1(i,j);
        hR = q1(i+1,j);
        uR = q2(i+1,j)/q1(i+1,j);
        vR = q3(i+1,j)/q1(i+1,j);

        hBar = 0.5*(hL+hR);
        uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
        vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
        cTilde = sqrt(g*hBar);
        FTmp = riemannX(uTilde, vTilde, cTilde, i, j, q1, q2, q3);
        F1(i) = FTmp(1);
        F2(i) = FTmp(2);
        F3(i) = FTmp(3);

      %end
    end


    for i=2:nX-1
      for j=1:nY
        q1(i,j) = q1(i,j) - dt/stepX*(F1(i)-F1(i-1));
        q2(i,j) = q2(i,j) - dt/stepX*(F2(i)-F2(i-1));
        q3(i,j) = q3(i,j) - dt/stepX*(F3(i)-F3(i-1));
      end
    end

    surf(x,y,q1);
    view(110,10);
    drawnow;
  end
  %surf(x,y,q1);
end

function F = riemannX(uTilde, vTilde, cTilde, i, j, q1, q2, q3)
  g = 9.81;
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
  alpha1 = ((uTilde+cTilde)*delta(1)-delta(2))/(2*cTilde);
  alpha2 = -vTilde*delta(1)+delta(3);
  alpha3 = (-(uTilde-cTilde)*delta(1)+delta(2))/(2*cTilde);
  lambda1 = uTilde-cTilde;
  lambda2 = uTilde;
  lambda3 = uTilde+cTilde;

  F(1) = 0.5*(q2(i,j)+q2(i+1,j));
  F(1) = F(1) - 0.5*(phi(lambda1)*alpha1*r1(1)
  + phi(lambda2)*alpha2*r2(1) +phi(lambda3)*alpha3*r3(1));

  F(2) = 0.5*(q2(i)^2/q1(i)+0.5*g*q1(i)^2 + q2(i+1)^2/q1(i+1)+0.5*g*q1(i+1)^2);
  F(2) = F(2) - 0.5*(phi(lambda1)*alpha1*r1(2)
  + phi(lambda2)*alpha2*r2(2) + phi(lambda3)*alpha3*r3(2));

  F(3) = 0.5*(q2(i,j)*q3(i,j)/q1(i,j)+q2(i+1,j)*q3(i+1,j)/q1(i+1,j)) ;
  F(3) = F(3) - 0.5*(phi(lambda1)*alpha1*r1(3)
  + phi(lambda2)*alpha2*r2(3) + phi(lambda3)*alpha3*r3(3));
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

