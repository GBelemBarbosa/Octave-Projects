res=[0.980, 0.988]
d=0.001
u_r1=d/(2*sqrt(3))

for r=res
  u_r2=(0.005*r+2*d)/(2*sqrt(3))
  u_r=sqrt(u_r1^2+u_r2^2)
endfor

I=[2.68, 3.34, 4.10, 4.82, 5.53]
d=0.01
u_i1=d/(2*sqrt(3))

for i=I
  u_i2=(0.005*i+3*d)/(2*sqrt(3));
  u_i=sqrt(u_i1^2+u_i2^2)
endfor

V=[3.810, 4.664, 5.621]
d=0.001
u_v1=d/(2*sqrt(3))

for v=V
  u_v2=(0.003*v+2*d)/(2*sqrt(3));
  u_v=sqrt(u_v1^2+u_v2^2)
endfor

V=[6.56, 7.49]
d=0.01
u_v1=d/(2*sqrt(3))

for v=V
  u_v2=(0.003*v+2*d)/(2*sqrt(3));
  u_v=sqrt(u_v1^2+u_v2^2)
endfor