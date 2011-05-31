function shallow1D
  clear all;
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
    if (x(i) <= Lx/2)
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
    printf('debut \n');
    for i=1:n-1
      hL = q1(i);
      uL = q2(i)/q1(i);
      hR = q1(i+1);
      uR = q2(i+1)/q1(i+1);


      hstar = (uL-uR+2*(sqrt(g*hL)+sqrt(g*hR)))^2/(16*g);
      for k=0:20
        if (hstar <= hL)
          FL = uL+2*sqrt(g*hL)-2*sqrt(g*hstar);
          FLprime = - sqrt(g/hstar);
        else
          FL = uL-(hstar-hL)*sqrt(g*(hstar+hL)/(2*hstar*hL));
          delta = sqrt(1/hL+1/hstar);
          %FLprime = sqrt(g/2)*(-delta + (hstar-hL)/(hstar^2*2*delta));
          FLprime = - g*(hL^2+hL*hstar+2*hstar^2)/(2*sqrt(2)*sqrt(g*hL*hstar^3*(hL+hstar)));
        end
        if (hstar <= hR)
          FR = uR-2*sqrt(g*hR)+2*sqrt(g*hstar);
          FRprime = sqrt(g/hstar);
        else
          FR = uR+(hstar-hR)*sqrt(g*(hstar+hR)/(2*hstar*hR));
          delta = sqrt(1/hR+1/hstar);
          %FRprime = sqrt(g/2)*(delta - (hstar-hR)/(hstar^2*2*delta));
          FRprime = g*(hR^2+hR*hstar+2*hstar^2)/(2*sqrt(2)*sqrt(g*hstar*hR^3*(hR+hstar)));
        end
        G = FL-FR;
        Gprime = FLprime-FRprime;
        hstar = hstar - G/Gprime;

      end
      hstar

      if (hstar <= hL) 
        ustar = uL+2*sqrt(g*hL)-2*sqrt(g*hstar);
      else 
        ustar = uL-(hstar-hL)*sqrt(g*(hstar+hL)/(2*hstar*hL));
      end

      if (hstar > hL) % shock wave 
        if (abs(hstar -hL) > 1e-6)
          sLeft1 = (hstar*ustar-hL*uL)/(hstar-hL);
        else
          sLeft1 = uL -sqrt(g*hL);
        end
        sLeft2 = sLeft1 ;
      else % rarefaction wave
        sLeft1 = uL-sqrt(g*hL);
        sLeft2 = ustar-sqrt(g*hstar);
      end
      if (hstar > hR) % shock wave 
        if (abs(hstar -hR) > 1e-6)
          sRight1 = (hstar*ustar-hR*uR)/(hstar-hR);
        else
          sRight1 = uR +sqrt(g*hR);
        end
        sRight2 = sRight1 ;
      else % rarefaction wave
        sRight1 = ustar+sqrt(g*hstar);
        sRight2 = uR+sqrt(g*hR);
      end

      % Goudunov's method
      if (sLeft1 > 0)
        hTmp = hL;
        uTmp = uL;
      elseif (sRight2 < 0)
        hTmp = hR;
        uTmp = uR;
      elseif (sLeft1 < 0 && sLeft2 > 0)
        hTmp = (uL+2*sqrt(g*hL))^2/(9*g);
        uTmp = uL+2*sqrt(g*hL)-2*sqrt(g*hTmp);
      elseif (sRight1 < 0 && sRight2 > 0)
        hTmp = (uR-2*sqrt(g*hR))^2/(9*g);
        uTmp = uL-2*sqrt(g*hR)+2*sqrt(g*hTmp);
      else
        hTmp = hstar;
        uTmp = hstar*ustar;
      end

      F1(i) = hTmp*uTmp;
      F2(i) = hTmp*uTmp^2+0.5*g*hTmp^2;
    end

    for i=2:n-1
      q1(i) = q1(i) - dt/dx*(F1(i)-F1(i-1));
      q2(i) = q2(i) - dt/dx*(F2(i)-F2(i-1));
    end
    q1
    plot(x,q1);
    drawnow;
  end

end
