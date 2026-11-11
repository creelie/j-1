#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* Theta(E) for g''+F g=0, F=E^2-1/4+K^2 cosh^2 t + 1/(4 sinh^2 t), regular at 0 */
static double K,E2;
static double F(double t){double s=sinh(t),c=cosh(t);return E2-0.25+K*K*c*c+0.25/(s*s);}
static double Fp(double t){double s=sinh(t),c=cosh(t);return 2*K*K*s*c-0.5*c/(s*s*s);}
static double rhs(double t,double ph){double f=F(t),k=sqrt(f),kp=Fp(t)/(2*k);return k+kp/(2*k)*sin(2*ph);}
int main(int argc,char**argv){
  double lam=atof(argv[1]); double E=atof(argv[2]); int odd=argc>3?atoi(argv[3]):0;
  K=2*M_PI*lam*lam; E2=E*E; double mu=-E2;
  double t0=1e-4; double x0=lam*cosh(t0);
  /* Frobenius in x: f=1+U1 (x-lam)+U2(x-lam)^2 ; ((x^2-l^2)f')'+(K'x^2-mu)f=0, K'=4pi^2 lam^2 */
  double Kp=4*M_PI*M_PI*lam*lam; double h=x0-lam;
  /* p=2 l h + h^2 ; (p f')' = p' f' + p f''. coefficients: order h^0: 2l U1 + (Kp l^2-mu)=0 */
  double U1=(mu-Kp*lam*lam)/(2*lam);
  /* order h^1: (p f')' : p=2l h+h^2, f'=U1+2U2 h, p f' = 2l U1 h + (U1+4 l U2) h^2 ; derivative: 2lU1 + 2(U1+4lU2) h
     (Kp x^2 - mu) f = (Kp l^2 - mu) + [2 Kp l + (Kp l^2-mu) U1] h ; sum h^1: 2(U1+4 l U2)+2Kp l+(Kp l^2-mu)U1=0 */
  double U2=-(2*U1+2*Kp*lam+(Kp*lam*lam-mu)*U1)/(8*lam);
  double f=1+U1*h+U2*h*h, fx=U1+2*U2*h;
  double s=sinh(t0),c=cosh(t0);
  double g=sqrt(lam*s)*f; double gp=0.5/sqrt(lam*s)*lam*c*f+sqrt(lam*s)*fx*lam*s;
  double k=sqrt(F(t0)); double ph=atan2(k*g,gp);
  double T=acosh(3e4/K); if(T<3)T=3;
  double t=t0;
  while(t<T){ double kk=sqrt(F(t)); double dt=0.02/kk; if(dt>1e-2)dt=1e-2; if(t+dt>T)dt=T-t;
    /* RK4 */
    double a1=rhs(t,ph),a2=rhs(t+dt/2,ph+dt/2*a1),a3=rhs(t+dt/2,ph+dt/2*a2),a4=rhs(t+dt,ph+dt*a3);
    ph+=dt/6*(a1+2*a2+2*a3+a4); t+=dt; }
  /* tail int_T^inf (k - K sinh) dt by substitution, simple adaptive Simpson on [T,T+40] */
  double tail=0; int n=20000; double L=25, hh=L/n;
  for(int i=0;i<=n;i++){double tt=T+i*hh; double sh=sinh(tt); double v=(E2-0.25+K*K+0.25/(sh*sh))/(sqrt(F(tt))+K*sh); double w=(i==0||i==n)?1:(i%2?4:2); tail+=w*v;}
  tail*=hh/3;
  double Th=ph+tail-K*cosh(T);
  /* I(E)=int_0^inf (sqrt(E2+K^2c^2)-K c) */
  double I=0; n=2; L=60; hh=L/n;
  for(int i=0;i<=n;i++){double tt=i*hh; double cc=cosh(tt); double v=E2/(sqrt(E2+K*K*cc*cc)+K*cc); double w=(i==0||i==n)?1:(i%2?4:2); I+=w*v;}
  I*=hh/3;
  double pred=E*log(4*E/K)-E;
  printf("%.6f %.10f %.10f %.10f %.10f\n",E,Th,I+M_PI/4,pred+M_PI/4,Th-(I+M_PI/4));
  return 0;}
