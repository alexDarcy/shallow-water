function shallow1D
  Lx = 10;
  n = 100;
  dx = Lx/n;
  dt = 0.01;
  g = 9.81;

  % Initialisation
  for i=1:n
    x(i) = i*dx/2;
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

  for t=0:dt:dt
    for i=1:n-1
      hL = q1(i);
      uL = q2(i)/q1(i);
      hR = q1(i+1);
      uR = q2(i+1)/q1(i+1);


      hstar = (uL-uR+2*(sqrt(g*hL)+sqrt(g*hR))^2)/(16*g);
      % TODO

      if (hstar > hL) 
        ustar = uL+2*sqrt(g*hL)-2*sqrt(g*hstar);
      else 
        ustar = uL-(hstar-hL)*sqrt(g*(hstar+hL)/(2*hstar*hL));
      end

      if (hstar > hL) % shock wave 
        sLeft1 = (hstar*ustar-hL*uL)/(hstar-hL);
        sLeft2 = sLeft1 ;
      else % rarefaction wave
        sLeft1 = uL-sqrt(g*hL);
        sLeft2 = ustar-sqrt(g*hstar);
      end
      if (hstar > hR) % shock wave 
        sRight1 = (hstar*ustar-hL*uL)/(hstar-hL);
        sRight2 = sRight1 ;
      else % rarefaction wave
        sRight1 = ustar-sqrt(g*hstar);
        sRight2 = uR-sqrt(g*hR);
      end

      if (sLeft1 > 0)
        q1(i) = hL;
        q2(i) = hL*uL;
      elseif (sRight1 > 0)
        q1(i) = hR;
        q2(i) = hR*uR;
      elseif (sLeft1 < 0 && sLeft2 > 0)
        hTmp = (uL+2*sqrt(g*hL))^2/(9*g);
        uTmp = uL+2*sqrt(g*hL)-2*sqrt(g*hTmp);
      elseif (sRight1 < 0 && sRight2 > 0)
        hTmp = (uR-2*sqrt(g*hR))^2/(9*g);
        uTmp = uL-2*sqrt(g*hL)+2*sqrt(g*hTmp);
      else
        q1(i) = hstar;
        q2(i) = hstar*ustar;
      end;


    end
    drawnow;
    plot(x,q1);
  end

end
