function shallow2D
  Lx = 1;
  Ly = 1;
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
    end
  end

  % Create mesh for visu
  xT = [];
  yT = [];
  hT = [];
  tri = [];
  k = 1;
  for j=1:nY-1
    for i=1:nX-1
      xT = [xT;x(i);x(i+1);x(i)];
      yT = [yT;y(j);y(j);y(j+1)];
      hT = [hT;h(i,j); h(i+1,j); h(i,j+1)];
      tri = [tri;k,k+1,k+2];

      xT = [xT;x(i+1);x(i+1);x(i)];
      yT = [yT;y(j);y(j+1);y(j+1)];
      hT = [hT;h(i+1,j); h(i+1,j+1); h(i,j+1)];
      tri = [tri;k+3,k+4,k+5];
      k += 6;
    end
  end

  trimesh(tri, xT,yT,hT);

  q1 = h;
  for i=1:nX
    for i=1:nY
      q2(i,j) = u(i,j)*h(i,j);
    end
  end

  for t=0:dt:dt
    for i=1:nX-1
      %for j=1:nY-1
      j = 1;
        hL = q1(i,j);
        uL = q2(i,j)/q1(i,j);
        hR = q1(i+1,j);
        uR = q2(i+1,j)/q1(i+1,j);

        hBar = 0.5*(hL+hR);
        uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
        cTilde = sqrt(g*hBar);

        r1(1) = 1;
        r1(2) = uTilde-cTilde;
        r2(1) = 1;
        r2(2) = uTilde+cTilde;
        delta(1) = q1(i+1,j)-q1(i,j);
        delta(2) = q2(i+1,j)-q2(i,j);
        alpha1 = ((uTilde+cTilde)*delta(1)-delta(2))/(2*cTilde);
        alpha2 = (-(uTilde-cTilde)*delta(1)+delta(2))/(2*cTilde);
        lambda1 = uTilde-cTilde;
        lambda2 = uTilde+cTilde;

        F1(i) = 0.5*(q2(i,j)+q2(i+1,j));
        F1(i) = F1(i) - 0.5*(phi(lambda1)*alpha1*r1(1)+ phi(lambda2)*alpha2*r2(1));

        F2(i) = 0.5*(q2(i)^2/q1(i)+0.5*g*q1(i)^2 + q2(i+1)^2/q1(i+1)+0.5*g*q1(i+1)^2);
        F2(i) = F2(i) - 0.5*(phi(lambda1)*alpha1*r1(2)+ phi(lambda2)*alpha2*r2(2));
      end
    end


    for i=2:nX-1
      for j=2:nY-1
        q1(i,j) = q1(i,j) - dt/stepX*(F1(i)-F1(i-1));
        q2(i,j) = q2(i,j) - dt/stepX*(F2(i)-F2(i-1));
      end
    end

    % Update mesh
    k = 1;
    for j=1:nY-1
      for i=1:nX-1
        hT(k) = q1(i,j); 
        hT(k+1) = q1(i+1,j); 
        hT(k+2) = q1(i,j+1);

        hT(k+3) = q1(i+1,j); 
        hT(k+4) = q1(i+1,j+1); 
        hT(k+5) = q1(i,j+1);
        k += 6;
      end
    end

    trimesh(tri, xT,yT,hT);
    drawnow;
  end
end

% Harten entropy fix
function z = phi(lambda)
  % empirical value
  epsilon = 2;
  if (abs(lambda) >= epsilon)
    z = abs(lambda);
  else
    z = (lambda^2 +
    silon^2)/(2*epsilon);
  end
end

