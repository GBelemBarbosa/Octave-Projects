n=11;
m=21;
Xspace=linspace(0,1,m);
A=vander(Xspace);
A=A(:,n:end);
Yspace=transpose(exp(sin(6.*Xspace)));

ATA=A'*A; 
R=chol(ATA); 
Q=A*inv(R); 
c=R\(Q'*Yspace); 

err=(Q*R-A).^2;
MSE=sum(err(:))/numel(A)

err=(Yspace-A*c).^2;
MSE=sum(err(:))/m

plot(Xspace, err)
title('Cholesky')
xlabel('x')
ylabel('[y(xi)-(A*c)(i)]^2')
figure()

plot(Xspace, Yspace, Xspace, A*c)
title('y=exp(sin(6x)) vs. Cholesky')
xlabel('x')
ylabel('y')

cond(ATA)
cond(R)