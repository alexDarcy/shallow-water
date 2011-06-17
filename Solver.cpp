#include "Solver.h"

Solver::Solver():g(9.81f){}

Solver::~Solver(){}
    
Vector3 Solver::riemannX(float* qL,float* qR)
{
  float hL = qL[1];
  float uL = qL[2]/qL[1];
  float vL = qL[3]/qL[1];
  float hR = qR[1];
  float uR = qR[2]/qR[1];
  float vR = qR[3]/qR[1];

  float hBar = 0.5*(hL+hR);
  float uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  float vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
  float cTilde = sqrt(g*hBar);

  Vector3 r1 (1, uTilde-cTilde, vTilde);
  Vector3 r2 (0, 0, 1);
  Vector3 r3 (1, uTilde+cTilde, vTilde);

  Vector3 alpha, lambda;
  Vector3 delta(qR[0],qR[1],qR[2]);
  for (int i = 0; i < 3; i++)
    delta(i) -= qL[i];
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
  F(1) = 0.5*(qL[2]+qR[2]);
  F(2) = 0.5*(qL[2]*qL[2]/qL[1]+0.5*g*qL[1]*qL[1] + qR[2]*qR[2]/qR[1]+0.5*g*qR[1]*qR[1]);
  F(3) = 0.5*(qL[2]*qL[3]/qL[1]+qR[2]*qR[3]/qR[1]) ;
  F = F - w;
  //lambdaMax = max(lambda);
  return F;
}

Vector3 Solver::riemannY(float* qL,float* qR)
{
  float hL = qL[1];
  float uL = qL[2]/qL[1];
  float vL = qL[3]/qL[1];
  float hR = qR[1];
  float uR = qR[2]/qR[1];
  float vR = qR[3]/qR[1];
  float hBar = 0.5*(hL+hR);
  float uTilde = (sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
  float vTilde = (sqrt(hL)*vL+sqrt(hR)*vR)/(sqrt(hL)+sqrt(hR));
 
  float cTilde = sqrt(g*hBar);

  Vector3 r1 (1, uTilde, vTilde-cTilde);
  Vector3 r2 (0, -1, 0);
  Vector3 r3 (1, uTilde, vTilde+cTilde);

  Vector3 delta(qR[0],qR[1],qR[2]);
  for (int i = 0; i < 3; i++)
    delta(i) -= qL[i];

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
  G(1) = 0.5*(qL[3]+qR[3]);
  G(2) = 0.5*(qL[2]*qL[3]/qL[1]+qR[2]*qR[3]/qR[1]) ;
  G(3) = 0.5*(qL[3]*qL[3]/qL[1]+0.5*g*qL[1]*qL[1] + qR[3]*qR[3]/qR[1]+0.5*g*qR[1]*qR[1]);

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

