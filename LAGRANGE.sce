function[L]=LAGRANGE(X,Y,x)
 L=0
 n=length(X)
 for i=1:n
  p=1
  for j=[1:i-1,i+1:n]
   a=x-X(j)
   b=X(i)-X(j)
   p=p*a/b
  end
  L=L+Y(i)*p
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
  y1(i)=LAGRANGE(X,Y,x1(i))
end

//Для отрезка
n=length(x2)
for i=1:n
 y2(i)=LAGRANGE(X,Y,x2(i))
end

//Отображение графиков
plot2d(x2,sin(x2),-3)
plot2d(x2,y2,-2)
plot2d(x1,y1,-4)
