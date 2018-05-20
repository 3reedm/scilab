y0=1  
x0=0
x=0:0.1:1
deff("[ydot]=f(x,y)","ydot=2*x^2+2*y")  
y=ode("adams",y0,x0,x,f)
plot2d(x,y,5)
plot2d(x,y,-5)
