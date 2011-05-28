function shallow2D
  Lx = 1;
  Ly = 1;
  nX = 3;
  nY = 3;
  stepX = Lx/nX;
  stepY = Ly/nY;
  
  % Create mesh
  x = 0:stepX:Lx;
  y = 0:stepY:Ly;
  for j=1:length(y)
    for i=1:length(x)
      z(i,j) = sin(x(i))*sin(y(i));
    end
  end
  
  % Create mesh for visu
  xT = [];
  yT = [];
  zT = [];
  tri = [];
  k = 1;
  for j=1:length(y)-1
    for i=1:length(x)-1
      xT = [xT;x(i);x(i+1);x(i)];
      yT = [yT;y(j);y(j);y(j+1)];
      zT = [zT;z(i,j); z(i+1,j); z(i,j+1)];
      tri = [tri;k,k+1,k+2];

      xT = [xT;x(i+1);x(i+1);x(i)];
      yT = [yT;y(j);y(j+1);y(j+1)];
      zT = [zT;z(i+1,j); z(i+1,j+1); z(i,j+1)];
      tri = [tri;k+3,k+4,k+5];
      k += 6;
    end
  end

  trimesh(tri, xT,yT,zT);
end
