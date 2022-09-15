n=11;
m=21;
Xspace=linspace(0,1,m);
A=vander(Xspace);
A=A(:,n:end);
Yspace=transpose(exp(sin(6.*Xspace)));
R=zeros(n);
Q=zeros(m,n);

for j=1:n
    v=A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j);
        v-=R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    if (R(j,j)==0)
        break
    end
    Q(:,j)=v/R(j,j);
end
c=R\(Q'*Yspace);

err=(Q*R-A).^2;
MSE=sum(err(:))/numel(A)

err=(Yspace-A*c).^2;
MSE=sum(err(:))/m

plot(Xspace, err)
title('GS')
xlabel('x')
ylabel('[y(xi)-(A*c)(i)]^2')
figure()

plot(Xspace, Yspace, Xspace, Acopia*c)
title('y=exp(sin(6x)) vs. GS')
xlabel('x')
ylabel('y')

cond(R)