function[P]=NEWTON(X,Y,x)
 P=Y(1)
 n=length(X)
 for i=2:n 
  A=0
  for j=1:i
   w=1
   for k=[1:j-1,j+1:i]
    w=w*(X(j)-X(k))
   end
   A=A+Y(j)/w
  end
  _X=1
  for l=1:i-1
   _X=_X*(x-X(l))
  end
  P=P+A*_X
 end
endfunction

//Исходные данные
X=[0,%pi/4,%pi/2,3*%pi/2]
Y=[1,1.4,1,-1]
x1=[%pi/6,%pi/3,%pi]
x2=[-3*%pi/2:%pi/8:3*%pi/2]

//Функция для 3 точек
n=length(x1)
for i=1:n
  y1(i)=NEWTON(X,Y,x1(i))
end

//Для отрезка
n=length(x2)
for i=1:n
 y2(i)=NEWTON(X,Y,x2(i))
end

//Отображение графиков
plot(x2,sin(x2),x2,y2)
//plot2d(x2,sin(x2),-3)
//plot2d(x2,y2,-2)
//plot2d(x1,y1,-4)
