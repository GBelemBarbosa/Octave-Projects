m=3
n=2
L=3 #number of layers (both ends included)
N=[n,2,1] #number of nodes (not including bias node) of each layer (both ends included)
alfa=1
A=cell(1, L);
x0=[0.1,0.3;0.6,-0.2;-1.2,-3.1];
Theta=cell(1, L);
Delta=cell(1, L);
for i=2:L
  Theta(i-1)=rand(N(i), N(i-1)+1);
  Delta(i-1)=zeros(N(i), N(i-1)+1);
end
Theta(L)=ones(1, N(L-1)+1)
Theta
Delta
Y=[1;0;0]

h = @(x) 1 ./(1+e.^(-x))

eps=0.001;
it_max=0;
k=0;

do 
  for i=1:m
    A=cell(1, L);
    A(1)=[1, x0(i, :)];
    for l=1:L-1
      A(l+1)=[ones(1, l!=L-1), h(cell2mat(A(l))*cell2mat(Theta(l)).')];
    end
    d=cell2mat(A(L))-Y(i)
    Delta(L-1)=cell2mat(Delta(L-1))+d.*cell2mat(A(l))
    d=((cell2mat(Theta(l))*d.').*Al.*(1-Al))
    for l=L-2:-1:1
      Delta(l)=cell2mat(Delta(l))+d(2:end).'*cell2mat(A(l));
      Al=cell2mat(A(l));
      d=((cell2mat(Theta(l))*d.').*Al.*(1-Al))
    endfor
  endfor
  k++;
until abs(err)<eps || k>it_max
