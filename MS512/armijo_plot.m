num_iter=1000;
function grad_x=gradiente(x)
    grad_x=[x(1)+2*x(1)^3-2*x(1)*x(2)-1; x(2)-x(1)^2];
endfunction
function hessian_x=hessian(x)
    hessian_x=[2*(x(1)^2-x(2))+4*x(1)^2+1, -2*x(1); -2*x(1), 1];
endfunction
f=@(x) ((x(1)^2-x(2))^2+(1-x(1))^2)/2;
iterations=0;
not_min=0;
#hf=figure()
#hold on;
tic();
for t=1:num_iter
    xk=transpose(rand(1,2));
    for i=1:2
        do
            aux=randi([-10,10]);
        until aux!=0
        xk(i)*=aux;
    endfor
    xk_0=xk;
    j=0;
    do 
        grad=gradiente(xk);
        normgrad=dot(grad,grad);
        lambda=normgrad/(transpose(grad)*hessian(xk)*grad);
        f_xk=f(xk);
        alpha_normgrad=0.5*normgrad;
        xk1=xk-lambda*grad;
        while (f_xk-f(xk1)<=alpha_normgrad*lambda)
            lambda*=0.8;
            xk1=xk-lambda*grad;
        endwhile
        xk=xk1;         
        j++;
    until (normgrad<=eps || lambda<=eps)
    if (normgrad<=eps)
        #scatter(xk_0(1), xk_0(2), "r");
        iterations+=j;
    else
        #scatter(xk_0(1), xk_0(2), "b");
        #scatter(xk(1), xk(2), "g");
        not_min++;
    endif
endfor
toc()
iterations
not_min
#print (hf, "armijo.png");
#hold off;
