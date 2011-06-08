function shallow2D(toPrint=0)
  Lx = 10;
  Ly = 10;
  n1 = 25;
  n2 = 25;
  dt = 0.01;
  g = 9.81;
  % Initial disk
  xC = Lx/2;
  yC = Ly/2;
  radius = 3;
  stepX = Lx/n1;
  stepY = Ly/n2;
  h1 = 1;
  h0 = 0.5;

  % Create mesh with ghost cells
  x = 0:stepX:Lx+stepX;
  y = 0:stepY:Ly+stepY;
  %z = zeros(nX+1,nY+1)
  nX = length(x);
  nY = length(y);
  % For visu
  limits = [0;Lx;0;Lx;0;1.5];

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

  q1 = h;
  for i=1:nX
    for j=1:nY
      q2(i,j) = u(i,j)*h(i,j);
      q3(i,j) = v(i,j)*h(i,j);
    end
  end

  xCut = x(1:nX-1);
  yCut = y(1:nY-1);
  hCut = q1(1:nX-1,1:nY-1);
  surf(xCut,yCut,hCut);
  view(110,10);
  axis(limits);
  if (toPrint > 0)
    print(strcat('movie/f',num2str(0),'.jpg'),'-djpg','-r200');
  else
    drawnow;
  end



  for k=1:150
    for i=1:nX-1
      for j=1:nY-1
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

        GTmp = riemannY(uTilde, vTilde, cTilde, i, j, q1, q2, q3);
        G1(i) = GTmp(1);
        G2(i) = GTmp(2);
        G3(i) = GTmp(3);

      end
    end


    for i=2:nX-1
      for j=2:nY-1
        q1(i,j) = q1(i,j) - dt/stepX*(F1(i)-F1(i-1));
        q1(i,j) = q1(i,j) - dt/stepY*(G1(i)-G1(i-1));

        q2(i,j) = q2(i,j) - dt/stepX*(F2(i)-F2(i-1));
        q2(i,j) = q2(i,j) - dt/stepY*(G2(i)-G2(i-1));

        q3(i,j) = q3(i,j) - dt/stepX*(F3(i)-F3(i-1));
        q3(i,j) = q3(i,j) - dt/stepY*(G3(i)-G3(i-1));
      end
    end

    hCut = q1(1:nX-1,1:nY-1);

    surf(xCut,yCut,hCut);
    axis(limits);
    view(110,10);
    if (toPrint > 0)
      print(strcat('movie/f',num2str(k),'.jpg'),'-djpg','-r200');
    else
      drawnow;
    end
  end
end

% Solves the Riemann problem on x
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

% Solves the Riemann problem on y
function G = riemannY(uTilde, vTilde, cTilde, i, j, q1, q2, q3)
  g = 9.81;
  r1(1) = 1;
  r1(2) = uTilde-cTilde;
  r1(3) = vTilde;
  r2(1) = 0;
  r2(2) = -1;
  r2(3) = 0;
  r3(1) = 1;
  r3(2) = uTilde;
  r3(3) = vTilde+cTilde;

  delta(1) = q1(i+1,j)-q1(i,j);
  delta(2) = q2(i+1,j)-q2(i,j);
  delta(3) = q3(i+1,j)-q3(i,j);
  alpha1 = ((vTilde+cTilde)*delta(1)-delta(3))/(2*cTilde);
  alpha2 = uTilde*delta(1)-delta(3);
  alpha3 = (-(vTilde-cTilde)*delta(1)+delta(3))/(2*cTilde);
  lambda1 = vTilde-cTilde;
  lambda2 = vTilde;
  lambda3 = vTilde+cTilde;

  G(1) = 0.5*(q3(i,j)+q3(i+1,j));
  G(1) = G(1) - 0.5*(phi(lambda1)*alpha1*r1(1)
  + phi(lambda2)*alpha2*r2(1) +phi(lambda3)*alpha3*r3(1));

  G(2) = 0.5*(q2(i,j)*q3(i,j)/q1(i,j)+q2(i+1,j)*q3(i+1,j)/q1(i+1,j)) ;
  G(2) = G(2) - 0.5*(phi(lambda1)*alpha1*r1(2)
  + phi(lambda2)*alpha2*r2(2) + phi(lambda3)*alpha3*r3(2));

  G(3) = 0.5*(q3(i)^2/q1(i)+0.5*g*q1(i)^2 + q3(i+1)^2/q1(i+1)+0.5*g*q1(i+1)^2);
  G(3) = G(3) - 0.5*(phi(lambda1)*alpha1*r1(3)
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

