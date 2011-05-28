function shallow1D
  clear all;
  g= 9.81;
  dt = 0.01;

  n = 90;
  L = 10;
  dx = L/n;
  for j=1:n
    x(j) = (j-1)*dx + dx/2;
  end
  % Init;
  n1 = int32(n/2),
  state = 0;
  for i=1:n
   % if (x(i) <= L/2)
   %   h(i) = 1;
   % else
      h(i) = 0.5;
    %end
  end
  u = zeros(1,n);
  q1 = h;
  q2 = u.*h;


  dh = 0.02;
  for t=dt:dt:2
    if (h(1) >= 1)
      state = 1;
    end
    if (h(1) <= 0.2)
      state = 0;
    end

    if (state == 1)
      h(1) = h(1) - dh;
    else 
      h(1) = h(1) + dh;
    end

    for i=1:n-1
      uL = u(i);
      uR = u(i+1);
      hL = h(i);
      hR = h(i+1);

      hstar = (uL-uR + 2*sqrt(g)*(sqrt(hL)+sqrt(hR)))^2/(16*g);
      for k=1:20
        if (hstar <= hL)
          FL = uL+2*sqrt(g*hL)-2*sqrt(g*hstar);
          FLprime = - sqrt(g/hstar);
          ondeL = 0;
        else
          FL = uL-(hstar-hL)*sqrt(g*(hstar+hL)/(2*hstar*hL));
          FLprime = - g*(hL^2+hL*hstar+2*hstar^2)/(2*sqrt(2)*sqrt(g*hL*hstar^3*(hL+hstar)));
          ondeL = 1;
        end
        if (hstar <= hR)
          FR = uR-2*sqrt(g*hR)+2*sqrt(g*hstar);
          FRprime = sqrt(g/hstar);
          ondeR = 0;
        else
          FR = uR+(hstar-hR)*sqrt(g*(hstar+hR)/(2*hstar*hR));
          FRprime = g*(hR^2+hR*hstar+2*hstar^2)/(2*sqrt(2)*sqrt(g*hstar*hR^3*(hR+hstar)));
          ondeR = 1;
        end
        G = FL - FR;
        Gprime = FLprime - FRprime;
        hstar = hstar - G/Gprime;
      end

      % shock wave
      speedLL = 0;
      speedRR = 0;
      if (ondeL == 1)
        ustar = uL-(hstar-hL)*sqrt(g*(hstar+hL)/(2*hstar*hL));
        if (abs(hL - hstar) > 1e-6)
          speedLL = (hL*uL-hstar*ustar)/(hL-hstar);
        else 
          speedLL = uL - sqrt(g*hL);
        end
        speedLR = speedLL;
      else
        ustar = uL+2*sqrt(g*hL)-2*sqrt(g*hstar);
        speedLL = uL - sqrt(g*hL);
        speedLR = ustar - sqrt(g*hstar);
      end
 
      if (ondeR == 1)
        ustar = uR+(hstar-hR)*sqrt(g*(hstar+hR)/(2*hstar*hR));
        if (abs(hR - hstar) > 1e-6)
          speedRR = (hstar*ustar-hR*uR)/(hstar-hR);
        else
          speedRR = uR + sqrt(g*hR);
        end
        speedRL = speedRR;
      else
        ustar = uR-2*sqrt(g*hR)+2*sqrt(g*hstar);
        speedRR = uR + sqrt(g*hR);
        speedRL = ustar + sqrt(g*hstar);
      end

      if (speedLL > 0)
        hk = hL;
        uk = uL;
      elseif (speedRR < 0)
        hk = hR;
        uk = uR;
      elseif (speedLL < 0 && speedLR > 0)
        hk = (uL + 2*sqrt(g*hL))/(9*g);
        uk = uL + 2*sqrt(g*hL) - 2*sqrt(g*hk);
      elseif (speedRL < 0 && speedRR > 0)
        hk = (uR - 2*sqrt(g*hR))/(9*g);
        uk = uR - 2*sqrt(g*hR) + 2*sqrt(g*hk);
      else
        hk = hstar;
        uk = ustar;
      end

      F1(i) = hk*uk;
      F2(i) = hk*uk^2 + 0.5*g*hk^2;
    end
    for i=2:n-1
      q1(i) = q1(i) - dt/dx*(F1(i)-F1(i-1));
      q2(i) = q2(i) - dt/dx*(F2(i)-F2(i-1));
    end
    for i=2:n-1
      h(i) = q1(i) ;
      u(i) = q2(i)/q1(i);
    end


    drawnow;
    plot(x,h);
    axis([0,L,0,1]);
    pause(0.02);
  end
end

