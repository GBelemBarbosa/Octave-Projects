m=3
n=2
alfa=1
theta=zeros(3,1)
X=[ones(m,1),[0.1,0.3;0.6,-0.2;-1.2,-3.1]]
Y=[1;0;0]

h = @(x) 1 ./(1+e.^(-x))

eps=0.001;
it_max=10000;
i=0;
do 
  err=h(X.'*theta)-Y;
  theta-=(alfa/m).*X.'*err;
  i++;
until abs(err)<eps || i>it_max

X.'*theta
h(X.'*theta)
