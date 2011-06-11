% 1D shallow water
% Approximate roe solver for Godunov's method with entropy fix

function shallow1D_roe
  Lx = 10;
  n = 100;
  dx = Lx/n;
  dt = 0.01;
  g = 9.81;

  % Initialisation
  for i=1:n
    x(i) = dx/2 + (i-1)*dx;
  end

  for i=1:n
    if (abs(x(i)-5) < 3)
      h(i) = 1;
    else
      h(i) = 0.5;
    end
    u(i) = 0;
  end

  for i=1:n
    q1(i) = h(i);
    q2(i) = u(i)*h(i);
  end

  for t=0:dt:100*dt
    for i=1:n-1
      hL = q1(i);
      uL = q2(i)/q1(i);
      hR = q1(i+1);
      uR = q2(i+1)/q1(i+1);

      hBar = 0.5*(hL+hR);
      uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
      cTilde = sqrt(g*hBar);

      r1(1) = 1;
      r1(2) = uTilde-cTilde;
      r2(1) = 1;
      r2(2) = uTilde+cTilde;
      delta(1) = q1(i+1)-q1(i);
      delta(2) = q2(i+1)-q2(i);
      alpha1 = ((uTilde+cTilde)*delta(1)-delta(2))/(2*cTilde);
      alpha2 = (-(uTilde-cTilde)*delta(1)+delta(2))/(2*cTilde);
      lambda1 = uTilde-cTilde;
      lambda2 = uTilde+cTilde;

      F1(i) = 0.5*(q2(i)+q2(i+1));
      F1(i) = F1(i) - 0.5*(phi(lambda1)*alpha1*r1(1)+ phi(lambda2)*alpha2*r2(1));

      F2(i) = 0.5*(q2(i)^2/q1(i)+0.5*g*q1(i)^2 + q2(i+1)^2/q1(i+1)+0.5*g*q1(i+1)^2);
      F2(i) = F2(i) - 0.5*(phi(lambda1)*alpha1*r1(2)+ phi(lambda2)*alpha2*r2(2));
    end


    for i=2:n-1
      q1(i) = q1(i) - dt/dx*(F1(i)-F1(i-1));
      q2(i) = q2(i) - dt/dx*(F2(i)-F2(i-1));
    end
    plot(x,q1);
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
    z = (lambda^2 + epsilon^2)/(2*epsilon);
  end
end

