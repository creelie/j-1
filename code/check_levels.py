import subprocess, math
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
BIN=os.path.join(HERE,'theta')
lam=2**0.5
def Th(E): return float(subprocess.check_output([BIN,repr(lam),repr(E)]).split()[1])
def theta_RS(t): # Riemann-Siegel theta asymptotic
    return t/2*math.log(t/(2*math.pi))-t/2-math.pi/8+1/(48*t)
E=0.5; prev=Th(E); out=[]
while E<25:
    E2=E+0.25; v=Th(E2)
    if math.floor(prev/math.pi)!=math.floor(v/math.pi):
        m=math.floor(v/math.pi); lo,hi=E,E2
        for _ in range(18):
            mid=(lo+hi)/2
            if Th(mid)/math.pi<m: lo=mid
            else: hi=mid
        t=lo+hi
        out.append((m,t,theta_RS(t)/math.pi))
    E,prev=E2,v
print("Theta(0.5)/pi",Th(0.5)/math.pi)
for m,t,r in out: print(m,round(t,6),round(r,6))
