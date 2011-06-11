% Finite volume solver for 2D shallow water equations
function shallow2D(toPrint=0)
  Lx = 9;
  Ly = 9;
  n1 = 20;
  n2 = 20;
  dt = 0.01;
  g = 9.81;
  % Initial disk
  xC = Lx/2; 
  yC = Ly/2;
  radius = 1;
  dx = Lx/n1;
  dy = Ly/n2;
  h1 = 1;
  h0 = 0.5;

  % Create mesh with ghost cells
  x = -dy:dx:Lx+dx;
  y = -dy:dy:Ly+dy;
  nX = length(x);
  nY = length(y);
  % For visu
  limits = [0;Lx;0;Ly;0;1];

  for i=1:nX
    for j=1:nY
      dist = (x(i)-xC)^2 + (y(j)-yC)^2;
      dist = sqrt(dist);
      if (dist <= radius )
        %if (x(i) > Lx/2 )
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
  if (toPrint > 0)
    f = figure('Visible','off');
  end
  surf(xCut,yCut,hCut);
  %surf(x,y,h);
  view(110,10);
  axis(limits);
  if (toPrint > 0)
    print(strcat('movie/f',num2str(0),'.jpg'),'-djpg','-r200');
  else
    drawnow;
  end

  qL = zeros(3,1);
  qRX = zeros(3,1);
  qRY = zeros(3,1);
  lMax = 0;
  for k=1:100
    for i=1:nX-1
      for j=1:nY-1
        % riemann x
        qL(1) = q1(i,j);
        qL(2) = q2(i,j);
        qL(3) = q3(i,j);
        qRX(1) = q1(i+1,j);
        qRX(2) = q2(i+1,j);
        qRX(3) = q3(i+1,j);
        qRY(1) = q1(i,j+1);
        qRY(2) = q2(i,j+1);
        qRY(3) = q3(i,j+1);

        [FTmp,l1] = riemann(1,qL,qRX);
        [GTmp,l2] = riemann(0,qL,qRY);
        lMax = max(lMax,max(l1,l2));

        F1(i,j) = FTmp(1);
        F2(i,j) = FTmp(2);
        F3(i,j) = FTmp(3);

        G1(i,j) = GTmp(1);
        G2(i,j) = GTmp(2);
        G3(i,j) = GTmp(3);
      end
    end

    froude = max(l1,l2)*dt/dx;
    printf('froude %f \n',froude);

    for i=1:nX-1
      for j=1:nY-1
        if (i > 1)
          q1(i,j) = q1(i,j) - dt/dx*(F1(i,j)-F1(i-1,j));%+F1tilde(i,j)-F1tilde(i-1,j));
          q2(i,j) = q2(i,j) - dt/dx*(F2(i,j)-F2(i-1,j));%+F2tilde(i,j)-F2tilde(i-1,j));
          h3(i,j) = q3(i,j) - dt/dx*(F3(i,j)-F3(i-1,j));%+F3tilde(i,j)-F3tilde(i-1,j));
        end

        if (j > 1)
          q1(i,j) = q1(i,j) - dt/dy*(G1(i,j)-G1(i,j-1));%+G1tilde(i,j)-G1tilde(i,j-1));
          q2(i,j) = q2(i,j) - dt/dy*(G2(i,j)-G2(i,j-1));%+G2tilde(i,j)-G2tilde(i,j-1));
          q3(i,j) = q3(i,j) - dt/dy*(G3(i,j)-G3(i,j-1));%+G3tilde(i,j)-G3tilde(i,j-1));
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

%G1tilde = zeros(nX,nY);
%    G2tilde = zeros(nX,nY);
%    G3tilde = zeros(nX,nY);
%    F1tilde = zeros(nX,nY);
%    F2tilde = zeros(nX,nY);
%    F3tilde = zeros(nX,nY);



%[BrAr, BlAr] = transverse(1,1,qL,qRX);
%[BrAl, BlAl] = transverse(1,0,qL,qRX);
%G1tilde(i+1,j+1) = G1tilde(i+1,j+1) - dt/(2*dx)*BrAr(1);
%G2tilde(i+1,j+1) = G2tilde(i+1,j+1) - dt/(2*dx)*BrAr(2);
%G3tilde(i+1,j+1) = G3tilde(i+1,j+1) - dt/(2*dx)*BrAr(3);

%G1tilde(i,j+1) = G1tilde(i,j+1) - dt/(2*dx)*BrAl(1);
%G2tilde(i,j+1) = G2tilde(i,j+1) - dt/(2*dx)*BrAl(2);
%G3tilde(i,j+1) = G3tilde(i,j+1) - dt/(2*dx)*BrAl(3);

%G1tilde(i+1,j) = G1tilde(i+1,j) - dt/(2*dx)*BlAr(1);
%G2tilde(i+1,j) = G2tilde(i+1,j) - dt/(2*dx)*BlAr(2);
%G3tilde(i+1,j) = G3tilde(i+1,j) - dt/(2*dx)*BlAr(3);

%G1tilde(i,j) = G1tilde(i,j) - dt/(2*dx)*BlAl(1);
%G2tilde(i,j) = G2tilde(i,j) - dt/(2*dx)*BlAl(2);
%G3tilde(i,j) = G3tilde(i,j) - dt/(2*dx)*BlAl(3);

%[BrAr, BlAr] = transverse(0,1,qL,qRY);
%[BrAl, BlAl] = transverse(0,0,qL,qRY);
%F1tilde(i+1,j+1) = F1tilde(i+1,j+1) - dt/(2*dy)*BrAr(1);
%F2tilde(i+1,j+1) = F2tilde(i+1,j+1) - dt/(2*dy)*BrAr(2);
%F3tilde(i+1,j+1) = F3tilde(i+1,j+1) - dt/(2*dy)*BrAr(3);

%F1tilde(i+1,j) = F1tilde(i+1,j) - dt/(2*dy)*BrAl(1);
%F2tilde(i+1,j) = F2tilde(i+1,j) - dt/(2*dy)*BrAl(2);
%F3tilde(i+1,j) = F3tilde(i+1,j) - dt/(2*dy)*BrAl(3);

%F1tilde(i+1,j) = F1tilde(i+1,j) - dt/(2*dy)*BlAr(1);
%F2tilde(i+1,j) = F2tilde(i+1,j) - dt/(2*dy)*BlAr(2);
%F3tilde(i+1,j) = F3tilde(i+1,j) - dt/(2*dy)*BlAr(3);

%G1tilde1(i+1,j) = G1tilde(i+1,j) - dt/(2*dy)*BlAl(1);
%G2tilde2(i+1,j) = G2tilde(i+1,j) - dt/(2*dy)*BlAl(2);
%G3tilde3(i+1,j) = G3tilde(i+1,j) - dt/(2*dy)*BlAl(3);


