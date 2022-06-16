num_iter=10;
function grad_x=gradiente(x)
    grad_x=[x(1)+2*x(1)^3-2*x(1)*x(2)-1; x(2)-x(1)^2];
endfunction
function hessian_x=hessian(x)
    hessian_x=[2*(x(1)^2-x(2))+4*x(1)^2+1, -2*x(1); -2*x(1), 1];
endfunction
f=@(x) ((x(1)^2-x(2))^2+(1-x(1))^2)/2;
tau=[0.92,0.7,0.8,0.75,0.85];
iterations=zeros(5,1);
iterations_lambda=zeros(5,1);
for t=1:num_iter
    t
    xk_0=transpose(rand(1,2));
    for i=1:2
        do
            aux=randi([-10,10]);
        until aux!=0
        xk_0(i)*=aux;
    endfor
    grad_0=gradiente(xk_0);
    normgrad_0=dot(grad_0,grad_0);
    for i=1:5
        i
        j=0;
        p=0;
        xk=xk_0;
        grad=grad_0;
        normgrad=normgrad_0;
        while (normgrad>eps && j<300)
            lambda=normgrad/(transpose(grad)*hessian(xk)*grad);
            f_xk=f(xk);
            aux=0.5*normgrad;          
            while (f_xk-f(xk-lambda*grad)<=aux*lambda)
                lambda*=tau(i);
                p++;
            endwhile
            xk-=lambda*grad;
            grad=gradiente(xk);
            normgrad=dot(grad,grad);
            j++;
        endwhile
        iterations(i)+=j
        iterations_lambda(i)+=p
    endfor
endfor
hf=figure;
xlabel("Proporcao de \lambda");
ylabel("Media de iteracoes");
hold on;
for (i=1:5)  
    scatter(tau(i), iterations(i)/num_iter, "filled");
endfor
print (hf, "grad.png");
hold off;
