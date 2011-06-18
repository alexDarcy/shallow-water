#include "Solver.h"

Solver::Solver():g(9.81f){}

Solver::Solver(float LX,float nX, float LZ,float nZ):g(9.81f),
  nbX(nX),nbZ(nZ),
  Lx(LX),Lz(LZ)
{
  dx = Lx / (nbX-1);
  dz = Lz / (nbZ-1);
  nbPoints = nbX*nbZ;
  nbQuads = (nbX-1)*(nbZ-1);

  float h1 = 1;
  float h0 = 0.5;
  float xC = Lx/2;
  float zC = Lz/2;
  float radius = Lx/4;
  float dist;
  int i = 0;
  dt = 0.02;
  h = new float[nbPoints];
  u = new float[nbPoints];
  v = new float[nbPoints];
  F1 = new float[nbPoints];
  F2 = new float[nbPoints];
  F3 = new float[nbPoints];
  G1 = new float[nbPoints];
  G2 = new float[nbPoints];
  G3 = new float[nbPoints];


  for (float z = 0; z <= Lz ; z += dz)
  {
    for (float x = 0; x <= Lx ; x += dx)
    {
      
      dist = (x-xC)*(x-xC) + (z-zC)*(z-zC);
      dist = sqrt(dist);
      if (dist <= radius ) // water drop
        h[i] = h1;
      else
        h[i] = h0;

      u[i] = 0;
      v[i] = 0;
      
      i++;
    }
  }
  q1 = new float[nbPoints];
  q2 = new float[nbPoints];
  q3 = new float[nbPoints];
  for (int i =0; i < nbPoints; i++)
  {
    q1[i] = h[i];
    q2[i] = u[i]*h[i];
    q3[i] = v[i]*h[i];
  }


}

Solver::~Solver(){
}

void Solver::run(float t)
{
  t++;
  Vector3 qL, qRX, qRY;
  for (int j =0; j < nbZ-1; j++)
  {
    for (int i =0; i < nbX-1; i++)
    {
      qL(1) = q1[i+j*nbX];
      qL(2) = q2[i+j*nbX];
      qL(3) = q3[i+j*nbX];
      qRX(1) = q1[i+1+j*nbX];
      qRX(2) = q2[i+1+j*nbX];
      qRX(3) = q3[i+1+j*nbX];
      qRY(1) = q1[i+j*nbX+nbX];
      qRY(2) = q2[i+j*nbX+nbX];
      qRY(3) = q3[i+j*nbX+nbX];

      Vector3 FTmp = riemannX(qL,qRX);
      Vector3 GTmp = riemannY(qL,qRY);

      F1[i+j*nbX] = FTmp(1);
      F2[i+j*nbX] = FTmp(2);
      F3[i+j*nbX] = FTmp(3);

      G1[i+j*nbX] = GTmp(1);
      G2[i+j*nbX] = GTmp(2);
      G3[i+j*nbX] = GTmp(3);
    }
  }

  for (int j =0; j < nbZ-1; j++)
  {
    for (int i =0; i < nbX-1; i++)
    {

      if (i > 0)
      {
        // if ((abs(y(j)-5) < 1) || i < limit || i > limit+1) % dam break
        q1[i+j*nbX] = q1[i+j*nbX] - dt/dx*(F1[i+j*nbX]-F1[i-1+j*nbX]);//%+F1tilde(i,j)-F1tilde(i-1,j));
        q2[i+j*nbX] = q2[i+j*nbX] - dt/dx*(F2[i+j*nbX]-F2[i-1+j*nbX]);//%+F2tilde(i,j)-F2tilde(i-1,j));
        q3[i+j*nbX] = q3[i+j*nbX] - dt/dx*(F3[i+j*nbX]-F3[i-1+j*nbX]);//%+F3tilde(i,j)-F3tilde(i-1,j));
      }

      if (j > 0)
      {
        //% if ((abs(y(j)-5) < 1) || i < limit || i > limit+1) % dam break
        q1[i+j*nbX] = q1[i+j*nbX] - dt/dz*(G1[i+j*nbX]-G1[i+j*nbX-nbX]);//%+G1tilde(i,j)-G1tilde(i,j-1));
        q2[i+j*nbX] = q2[i+j*nbX] - dt/dz*(G2[i+j*nbX]-G2[i+j*nbX-nbX]);//%+G2tilde(i,j)-G2tilde(i,j-1));
        q3[i+j*nbX] = q3[i+j*nbX] - dt/dz*(G3[i+j*nbX]-G3[i+j*nbX-nbX]);//%+G3tilde(i,j)-G3tilde(i,j-1));
      }
    }
  }
}


Vector3 Solver::riemannX(Vector3& qL,Vector3& qR)
{
  float hL = qL(1);
  float uL = qL(2)/qL(1);
  float vL = qL(3)/qL(1);
  float hR = qR(1);
  float uR = qR(2)/qR(1);
  float vR = qR(3)/qR(1);

  float hBar = 0.5*(hL+hR);
  float uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  float vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
  float cTilde = sqrt(g*hBar);

  Vector3 r1 (1, uTilde-cTilde, vTilde);
  Vector3 r2 (0, 0, 1);
  Vector3 r3 (1, uTilde+cTilde, vTilde);

  Vector3 alpha, lambda;
  Vector3 delta = qR;
  delta -= qL;
  Vector3 w;
  alpha(1) = ((uTilde+cTilde)*delta(1)-delta(2))/(2*cTilde);
  alpha(2) = -vTilde*delta(1)+delta(3);
  alpha(3) = (-(uTilde-cTilde)*delta(1)+delta(2))/(2*cTilde);
  lambda(1) = uTilde-cTilde;
  lambda(2) = uTilde;
  lambda(3) = uTilde+cTilde;
  w = r1*phi(lambda(1))*alpha(1);
  w += r2*phi(lambda(2))*alpha(2);
  w += r3*phi(lambda(3))*alpha(3);
  w *= 0.5;

  Vector3 F;
  F(1) = 0.5*(qL(2)+qR(2));
  F(2) = 0.5*(qL(2)*qL(2)/qL(1)+0.5*g*qL(1)*qL(1) + qR(2)*qR(2)/qR(1)+0.5*g*qR(1)*qR(1));
  F(3) = 0.5*(qL(2)*qL(3)/qL(1)+qR(2)*qR(3)/qR(1)) ;
  F = F - w;
  //lambdaMax = max(lambda);
  return F;
}

Vector3 Solver::riemannY(Vector3& qL,Vector3& qR)
{
  float hL = qL(1);
  float uL = qL(2)/qL(1);
  float vL = qL(3)/qL(1);
  float hR = qR(1);
  float uR = qR(2)/qR(1);
  float vR = qR(3)/qR(1);
  float hBar = 0.5*(hL+hR);
  float uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  float vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));

  float cTilde = sqrt(g*hBar);

  Vector3 r1 (1, uTilde, vTilde-cTilde);
  Vector3 r2 (0, -1, 0);
  Vector3 r3 (1, uTilde, vTilde+cTilde);

  Vector3 delta = qR;
  delta -= qL;

  Vector3 alpha;
  Vector3 lambda;
  Vector3 w;
  alpha(1) = ((vTilde+cTilde)*delta(1)-delta(3))/(2*cTilde);
  alpha(2) = uTilde*delta(1)-delta(2);
  alpha(3) = (-(vTilde-cTilde)*delta(1)+delta(3))/(2*cTilde);
  lambda(1) = vTilde-cTilde;
  lambda(2) = vTilde;
  lambda(3) = vTilde+cTilde;
  w = r1*phi(lambda(1))*alpha(1);
  w += r2*phi(lambda(2))*alpha(2);
  w += r3*phi(lambda(3))*alpha(3);
  w *= 0.5;



  Vector3 G;
  G(1) = 0.5*(qL(3)+qR(3));
  G(2) = 0.5*(qL(2)*qL(3)/qL(1)+qR(2)*qR(3)/qR(1)) ;
  G(3) = 0.5*(qL(3)*qL(3)/qL(1)+0.5*g*qL(1)*qL(1) + qR(3)*qR(3)/qR(1)+0.5*g*qR(1)*qR(1));

  G = G - w;
  return G;

  //lambdaMax = max(lambda);
}

// Harten entropy fix
float Solver::phi(float lambda)
{
  //% empirical value
  //  %epsilon = 2;
  //%if (abs(lambda) >= epsilon)
  return abs(lambda);
  //%else
  //  %   z = (lambda^2 + epsilon^2)/(2*epsilon);
  //%end
}
