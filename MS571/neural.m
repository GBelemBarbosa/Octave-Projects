m=3
n=3
global L=3 #number of layers (both ends included)
global N=[n-1,2,1] #number of nodes (not including bias node) of each layer (both ends included)
Y=[1;0;1]
x0=[0.1,0.3;0.6,-0.2;-1.2,-3.1];
alfa=1
eps=0.001;
it_max=1;

h = @(x) 1 ./(1+e.^(-x));

function x0=calc(x0, Theta, h)
  global L N;
  for l=1:L-1
    x0=h([ones(N(l)+1, 1), x0]*cell2mat(Theta(l)).');
  end
end  

Theta=cell(1, L);
for i=2:L
  Theta(i-1)=rand(N(i), N(i-1)+1);
end
Theta(L)=ones(1, N(L-1)+1);

k=0;
err=0;
do 
  Delta=cell(1, L-1);
  for l=2:L
    Delta(l-1)=zeros(N(l), N(l-1)+1);
  end
  err=0;
  for i=1:m
    A=cell(1, L);
    A(1)=x0(i, :);
    for l=1:L-1
      A(l+1)=h([1, cell2mat(A(l))]*cell2mat(Theta(l)).');
    end
    err=max(err, max(abs(cell2mat(A(L))-Y(i))));
    d=cell2mat(A(L))-Y(i);
    err=d;
    Al=cell2mat(A(L-1));
    Delta(L-1)
    d.*[1, Al]
    Delta(L-1)=cell2mat(Delta(L-1))+d.*[1, Al];
    d=((cell2mat(Theta(L-1))(2:end).*d).*Al.*(1-Al));
    for l=L-2:-1:1
      Al=cell2mat(A(l));
      d, [1, Al]
      Delta(l)=cell2mat(Delta(l))+d.'*[1, Al];
      d=((cell2mat(Theta(l))(:, 2:end)*d.').*Al.*(1-Al));
    endfor
  endfor
  for l=1:L-1
    Theta(l)=cell2mat(Theta(l))-(alfa/m).*cell2mat(Delta(l));
  endfor
  k++;
until abs(err)<eps || k>it_max
