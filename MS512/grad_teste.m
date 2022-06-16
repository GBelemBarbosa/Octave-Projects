c=1000;
G=c*diag(1:10);
G(1,1)=1
iterations=zeros(5,1);
scale=[0.1, 0.3, 0.5, 0.75, 1];
num_iter=100;
for t=1:num_iter
    t
    xk_0=transpose(rand(1,10));
    for i=1:10
        do
            aux=randi([-10,10]);
        until aux!=0
        xk_0(i)*=aux;
    endfor
    grad_0=G*xk_0;
    normgrad_0=dot(grad_0,grad_0);
    for i=1:5
        j=0;
        xk=xk_0;
        grad=grad_0;
        normgrad=normgrad_0;
        while (normgrad>eps)
            lambda=normgrad/(transpose(grad)*G*grad);
            xk-=scale(i)*lambda*grad;
            grad=G*xk;
            normgrad=dot(grad,grad);
            j++;
        endwhile
        iterations(i)+=j;
    endfor
endfor
hf=figure;
xlabel("Proporcao de \\lambda");
ylabel("Media de iteracoes");
hold on;
for (i=1:5)  
    scatter(scale(i), iterations(i)/num_iter, "filled");
endfor
print (hf, "grad1000.png");
hold off;
