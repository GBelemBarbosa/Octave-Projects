num_iter=100;
function grad_x=gradiente(x)
    grad_x=[x(1)+2*x(1)^3-2*x(1)*x(2)-1; x(2)-x(1)^2];
endfunction
laplacian=@(x) 6*x(1)^2-2*x(2)+2;
f=@(x) ((x(1)^2-x(2))^2+(1-x(1))^2)/2;
iterations=0;
for t=1:num_iter
    t
    xk=transpose(rand(1,2));
    for i=1:2
        do
            aux=randi([-10,10]);
        until aux!=0
        xk(i)*=aux;
    endfor
    j=0;
    do 
        grad=gradiente(xk);
        normgrad=dot(grad,grad);
        xk-=grad/laplacian(xk);           
        j++;
    until (normgrad<=eps || j>10000)
    if (normgrad<=eps)
        iterations+=j;
    endif
endfor
iterations