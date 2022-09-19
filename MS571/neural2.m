m=3
n=2
L=3 #number of layers (both ends included)
N=[n,2,1] #number of nodes (not including bias node) of each layer (both ends included)
alfa=1
A=cell(1, L);
A(1)=[0.1,0.3;0.6,-0.2;-1.2,-3.1];
Theta=cell(1, L);
for i=2:L
  Theta(i-1)=rand(N(i), N(i-1)+1);
end
Theta(L)=ones(1, N(L-1)+1)
Delta=cell(1, L-1);
Ds=cell(1, m);
D=cell(1, L-1);
Theta
Y=[1;0;0]

h = @(x) 1 ./(1+e.^(-x))

eps=0.001;
it_max=0;
k=0;

do 
  for l=1:L-1
    A(l+1)=[h([ones(m, 1), cell2mat(A(l))]*cell2mat(Theta(l)).')];
  endfor
  d=cell2mat(A(L))-Y;
  for l=L:-1:2
    D=zeros(N(l), N(l-1)+1);
    Al=cell2mat(A(l));
    Al1=cell2mat(A(l-1));
    for i=1:m
      dAux=(cell2mat(Theta(l))*d).*Al(i,:).*(1-Al(i,:));
      D+=dAux.'*[1, Al1(i,:)];
    end 
    Delta(l-1)=D;   
    d
  endfor
  D
  k++;
until abs(err)<eps || k>it_max