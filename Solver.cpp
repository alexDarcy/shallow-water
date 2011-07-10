#include "Solver.h"

Solver::Solver():g(9.81f){}

Solver::Solver(float LX,float nX, float LZ,float nZ):g(9.81f)
{
  isDamBreak = true;
  isWaterDrop = false;//true;
  isTsunami= false;

  dx = LX / (nX-1);
  dz = LZ / (nZ-1);

  /* Ghost cells */
  nbX = nX+2;
  nbZ = nZ+2;
  Lx = LX+dx; 
  Lz = LZ+dz; 
  nbPoints = nbX*nbZ;
  nbQuads = (nbX-1)*(nbZ-1);

  if (isDamBreak)
    h1 = 2;
  else
    h1 = 3;

  h0 = 0.5;
  xC1 = LX/3;
  zC1 = Lz/3;
  xC2 = 2*LX/3;
  zC2 = 2*Lz/3;
  radius = LX/10;
  float dist1;
  float dist2;
  int i = 0;
  dt = 0.007;
  F1 = new float[nbPoints];
  F2 = new float[nbPoints];
  F3 = new float[nbPoints];
  G1 = new float[nbPoints];
  G2 = new float[nbPoints];
  G3 = new float[nbPoints];
  /* Dam properties */
  float damL = Lx/4;
  float damW = 2*dx;
  float holeL = Lz/2;
  float holeW = 5*dz;
  damLimit = damL / dx;
  damWidth = damW / dx;

  holeLimit = holeL / dz;
  holeWidth = holeW / dz;
  bool cond;

  q1 = new float[nbPoints];
  q2 = new float[nbPoints];
  q3 = new float[nbPoints];

  for (float z = -dz ; z <= Lz ; z += dz)
  {
    for (float x = -dx; x <= Lx ; x += dx)
    {
      if (isDamBreak)
      {
        //cout << i << " " << damLimit << " " << damWidth << endl;
        cond = (i % nbX) <= damLimit+damWidth;
      }
      else if (isWaterDrop)
      {
        dist1 = (x-xC1)*(x-xC1) + (z-zC1)*(z-zC1);
        dist1 = sqrt(dist1);
        dist2 = (x-xC2)*(x-xC2) + (z-zC2)*(z-zC2);
        dist2 = sqrt(dist2);

        cond =dist1 <= radius;// || dist2 <= radius ; 
      }
      else 
        cond = false;

      if (cond) 
      {
        q1[i] = h1;
      }
      else
      {
        q1[i] = h0;
      }

      q2[i] = 0;
      q3[i] = 0;

      i++;
    }
  }
}

Solver::~Solver(){
    delete[] q1;
    delete[] q2;
    delete[] q3;
 
}

/* Compute flux and update Q */
void Solver::run(float t)
{
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

  bool cond;
  bool cond1;
  bool cond2;
  for (int j =0; j < nbZ-1; j++)
  {
    for (int i =0; i < nbX-1; i++)
    {
      // dam break
      if (isDamBreak)
      {
        cond = (abs(j-holeLimit) <= holeWidth*0.5); 
        cond1 = i > damLimit + damWidth+1 ; // bottom of the dam
        cond2 = i < damLimit - damWidth; 
        cond = cond || cond1 || cond2;
      }
      else 
        cond = true;

      if (cond) 
      {
        if (i > 0 && j > 0)
        {
          q1[i+j*nbX] = q1[i+j*nbX] - dt/dx*(F1[i+j*nbX]-F1[i-1+j*nbX]);
          q2[i+j*nbX] = q2[i+j*nbX] - dt/dx*(F2[i+j*nbX]-F2[i-1+j*nbX]);
          q3[i+j*nbX] = q3[i+j*nbX] - dt/dx*(F3[i+j*nbX]-F3[i-1+j*nbX]);

          q1[i+j*nbX] = q1[i+j*nbX] - dt/dz*(G1[i+j*nbX]-G1[i+(j-1)*nbX]);
          q2[i+j*nbX] = q2[i+j*nbX] - dt/dz*(G2[i+j*nbX]-G2[i+(j-1)*nbX]);
          q3[i+j*nbX] = q3[i+j*nbX] - dt/dz*(G3[i+j*nbX]-G3[i+(j-1)*nbX]);
        }
      }
    }
  }

  boundary();
  if (isWaterDrop)
    addDrops(t);
}

// Set boundary conditions
void Solver::boundary()
{
  int j0;
  int i0;
  float dh = (h1-h0)/10;
  float up = 1.f;

  float mean = 0;
  
  for (int j = 0; j < nbZ; j++)
    mean += q1[j*nbX];

  mean /= nbZ;

  if (isTsunami)
  {
    if (up > 0 && mean >= h1)
      up = -1.f;
    else if (up < 0 && mean <= h0)
      up = 1.f;
    
    for (int j = 0; j < nbZ; j++)
    {
        q1[j*nbX] += dh*up;
    }
  }
  else if (isWaterDrop)
  {
    /* Reflection on the boundary */
    for (int j = 1; j < nbZ-1; j++)
    {
      i0 = 1;
      q2[i0+j*nbX] = -q2[i0+1+j*nbX];
      i0 = nbX-2;
      q2[i0+j*nbX] = -q2[i0-1+j*nbX];
    }
    for (int i = 1; i < nbX-1; i++)
    {
      j0 = 1;
      q2[i+j0*nbX] = -q2[i+(j0+1)*nbX];
      j0 = nbZ-2;
      q2[i+j0*nbX] = -q2[i+(j0-1)*nbX];
    }
  }
  else
  {
    i0 = nbX-1;
    for (int j = 0; j < nbZ; j++)
    {
      q1[i0+j*nbX] = q1[i0+j*nbX-1];
      q2[i0+j*nbX] = q2[i0+j*nbX-1];
      q3[i0+j*nbX] = q3[i0+j*nbX-1];
    }
  }
}

void Solver::addDrops(float t)
{
  float latency = 2.f;
  float precision = 1e-2f;
  float tmp = t - (int) t;
  //cout << t << " / " << tmp << endl;
  if (tmp < precision)
  {
    int iC = xC1/dx;
    int jC = zC1/dz;
    int iWidth = radius/dx;
    int jWidth = radius/dz;
    for (int i = iC-iWidth/2; i<= iC+iWidth/2; i++)
      for (int j = jC-jWidth/2; j<= jC+jWidth/2; j++)
      {
        q1[i+j*nbX] = h1;
        q2[i+j*nbX] = 0;
        q3[i+j*nbX] = 0;
      }
  }
}

/* Solve the riemann problem on X */
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

/* Solve the riemann problem on Y */
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
  // empirical value
  float epsilon = 2;
  if (abs(lambda) >= epsilon)
    return abs(lambda);
  else
    return(lambda*lambda + epsilon*epsilon)/(2*epsilon);
}
