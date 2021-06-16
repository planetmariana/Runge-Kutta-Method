import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

Carga_masa=-1.76E11
B=5.0E-5
lam=-1.0E-11
epsilon=9.0E-12

def x_prime(x,t,y,v_x,v_y):
	return v_x


def y_prime(x,t,y,v_x,v_y):
	return v_y

def vx(x,t,y,v_x,v_y):
	return (Carga_masa*B*v_y) + ((lam*Carga_masa*x)/(2.0*np.pi*epsilon*(x**2.0 + y**2.0)))

def vy(x,t,y,v_x,v_y):
	return  (-Carga_masa*B*v_x) + ((lam*Carga_masa*y)/(2.0*np.pi*epsilon*(x**2.0 + y**2.0)))


t0=0.0
t1=2E-5
puntos=2000
delta=(t1-t0)/puntos
t,x,y,dvx,dvy=np.zeros(puntos),np.zeros(puntos),np.zeros(puntos),np.zeros(puntos),np.zeros(puntos)


def RungeK(x1,tb,y1,vx1,vy1):
	k1x=x_prime(x1,tb,y1,vx1,vy1)
	k1y=y_prime(x1,tb,y1,vx1,vy1)
	k1vx=vx(x1,tb,y1,vx1,vy1)
	k1vy=vy(x1,tb,y1,vx1,vy1)



	k2x=x_prime(x1 + k1x*delta*0.5, tb + delta*0.5, y1 + k1y*delta*0.5,vx1 + k1vx*delta*0.5, vy1 + k1vy*delta*0.5)
	k2y=y_prime(x1 + k1x*delta*0.5, tb + delta*0.5, y1 + k1y*delta*0.5,vx1 + k1vx*delta*0.5, vy1 + k1vy*delta*0.5)
	k2vx=vx(x1 + k1x*delta*0.5, tb + delta*0.5, y1 + k1y*delta*0.5,vx1 + k1vx*delta*0.5, vy1 + k1vy*delta*0.5)
	k2vy=vy(x1 + k1x*delta*0.5, tb + delta*0.5, y1 + k1y*delta*0.5,vx1 + k1vx*delta*0.5, vy1 + k1vy*delta*0.5)



	k3x=x_prime(x1 + k2x*delta*0.5, tb + delta*0.5, y1 + k2y*delta*0.5,vx1 + k2vx*delta*0.5, vy1 + k2vy*delta*0.5)
	k3y=y_prime(x1 + k2x*delta*0.5, tb + delta*0.5, y1 + k2y*delta*0.5,vx1 + k2vx*delta*0.5, vy1 + k2vy*delta*0.5)
	k3vx=vx(x1 + k2x*delta*0.5, tb + delta*0.5, y1 + k2y*delta*0.5,vx1 + k2vx*delta*0.5, vy1 + k2vy*delta*0.5)
	k3vy=vy(x1 + k2x*delta*0.5, tb + delta*0.5, y1 + k2y*delta*0.5,vx1 + k2vx*delta*0.5, vy1 + k2vy*delta*0.5)



	k4x=x_prime(x1 + k3x*delta, tb + delta, y1 + k3y*delta,vx1 + k3vx*delta, vy1 + k3vy*delta)
	k4y=y_prime(x1 + k3x*delta, tb + delta, y1 + k3y*delta,vx1 + k3vx*delta, vy1 + k3vy*delta)
	k4vx=vx(x1 + k3x*delta, tb + delta, y1 + k3y*delta,vx1 + k3vx*delta, vy1 + k3vy*delta)
	k4vy=vy(x1 + k3x*delta, tb + delta, y1 + k3y*delta,vx1 + k3vx*delta, vy1 + k3vy*delta)
		


	promedio_x_prime=(k1x+2.0*k2x+2.0*k3x+k4x)/6.0
	x2=x1 +promedio_x_prime*delta
	promedio_y_prime=(k1y+2.0*k2y+2.0*k3y+k4y)/6.0
	y2=y1+promedio_y_prime*delta
	promedio_vx=(k1vx+2.0*k2vx+2.0*k3vx+k4vx)/6.0
	vx2=vx1+promedio_vx*delta
	promedio_vy=(k1vy+2.0*k2vy+2.0*k3vy+k4vy)/6.0
	vy2=vy1+promedio_vy*delta
	tb2=tb+delta




	return x2,tb2,y2,vx2,vy2

t[0],x[0],y[0],dvx[0],dvy[0]=t0,2.0,2.0,0.0,-3000

for i in range(1,2000):
	x[i],t[i],y[i],dvx[i],dvy[i]=RungeK(x[i-1],t[i-1],y[i-1],dvx[i-1],dvy[i-1])


#Siguiendo el enunciado, se guarda la grafica en un .pdf
plt.plot(y,x,c="DeepSkyBlue",label="Runge-Kutta")
plt.legend(loc=0)
plt.ylabel("y(m)")
plt.xlabel("x(m)")
plt.savefig("Grafica-Runge-Kutta.pdf")
plt.close()






































