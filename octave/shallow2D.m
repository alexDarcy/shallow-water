% Finite volume solver for 2D shallow water equations
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
  x = -stepX:stepX:Lx+stepX;
  y = -stepY:stepY:Ly+stepY;
  %z = zeros(nX+1,nY+1)
  nX = length(x);
  nY = length(y);
  % For visu
  limits = [0;Lx;0;Ly;0;1.5];

  for i=1:nX
    for j=1:nY
      %dist = (x(i)-xC)^2 + (y(j)-yC)^2;
      %dist = sqrt(dist);
      %if (dist <= radius )
      %if (abs(x(i)-xC) <= 3 && abs(y(j)-yC) <= 3)
      if (x(i) < Lx/2 )
      %if (x(i) < Lx/2 || y(j) < Ly/2)
        h(i,j) = h1;
      else
        h(i,j) = h0;
      end
    end
    u(i,j) = 0;
    v(i,j) = 0;
  end

  q1 = h;
  for i=1:nX
    for j=1:nY
      q2(i,j) = u(i,j)*h(i,j);
      q3(i,j) = v(i,j)*h(i,j);
    end
  end

  xCut = x(2:nX-1);
  yCut = y(2:nY-1);
  hCut = q1(2:nX-1,2:nY-1);
  surf(xCut,yCut,hCut);
  %surf(x,y,h);
  view(110,10);
  axis(limits);
  if (toPrint > 0)
    print(strcat('movie/f',num2str(0),'.jpg'),'-djpg','-r200');
  else
    drawnow;
  end



  for k=1:100
    for i=1:nX-1
      for j=1:nY-1
        % riemann x
        FTmp = riemann(1,i, j, q1, q2, q3);
        F1(i) = FTmp(1);
        F2(i) = FTmp(2);
        F3(i) = FTmp(3);

        GTmp = riemann(0,i, j, q1, q2, q3);
        G1(j) = GTmp(1);
        G2(j) = GTmp(2);
        G3(j) = GTmp(3);

      end
    end


    for i=1:nX-1
      for j=1:nY-1
        if (i > 1)
          q1(i,j) = q1(i,j) - dt/stepX*(F1(i)-F1(i-1));
          q2(i,j) = q2(i,j) - dt/stepX*(F2(i)-F2(i-1));
          q3(i,j) = q3(i,j) - dt/stepX*(F3(i)-F3(i-1));
        end

        if (j > 1)
          q1(i,j) = q1(i,j) - dt/stepY*(G1(j)-G1(j-1));
          q2(i,j) = q2(i,j) - dt/stepY*(G2(j)-G2(j-1));
          q3(i,j) = q3(i,j) - dt/stepY*(G3(j)-G3(j-1));
        end
      end
    end

    boundary(q1,q2,q3,nX,nY);

    hCut = q1(2:nX-1,2:nY-1);
    surf(xCut,yCut,hCut);
    %surf(x,y,h);
    axis(limits);
    view(110,10);
    if (toPrint > 0)
      print(strcat('movie/f',num2str(k),'.jpg'),'-djpg','-r200');
    else
      drawnow;
    end
  end
end

% interpolate boundary conditions
function boundary(q1,q2,q3,nX,nY)
  for j=1:nY
    q1(1,j) = 2*q1(2,j) - q1(3,j);
    q2(1,j) = 2*q2(2,j) - q2(3,j);
    q3(1,j) = 2*q3(2,j) - q3(3,j);

    q1(nX,j) = 2*q1(nX-1,j) - q1(nX-2,j);
    q2(nX,j) = 2*q2(nX-1,j) - q2(nX-2,j);
    q3(nX,j) = 2*q3(nX-1,j) - q3(nX-2,j);
  end
  for i=1:nX
    q1(i,1) = 2*q1(i,2) - q1(i,3);
    q2(i,1) = 2*q2(i,2) - q2(i,3);
    q3(i,1) = 2*q3(i,2) - q3(i,3);

    q1(i,nY) = 2*q1(i,nY-1) - q1(i,nY-2);
    q2(i,nY) = 2*q2(i,nY-1) - q2(i,nY-2);
    q3(i,nY) = 2*q3(i,nY-1) - q3(i,nY-2);
  end
end



% Solves the Riemann problem on x
% and add transverse flux

