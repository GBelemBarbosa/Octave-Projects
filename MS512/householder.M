n=11;
m=21;
Xspace=linspace(0,1,m);
A=vander(Xspace);
A=A(:,n:end);
Yspace=transpose(exp(sin(6.*Xspace)));

[Q,R]=qr(A);
c=R\(Q'*Yspace);

err=(Q*R-A).^2;
MSE=sum(err(:))/numel(A)

err=(Yspace-A*c).^2;
MSE=sum(err(:))/m

plot(Xspace, err)
title('Householder')
xlabel('x')
ylabel('[y(xi)-(A*c)(i)]^2')
figure()

plot(Xspace, Yspace, Xspace, A*c)
title('y=exp(sin(6x)) vs. Householder')
xlabel('x')
ylabel('y')

cond(R)