num_iter=1000;
function grad_x=gradiente(x)
    grad_x=[x(1)+2*x(1)^3-2*x(1)*x(2)-1; x(2)-x(1)^2];
endfunction
laplacian=@(x) 6*x(1)^2-2*x(2)+2;
f=@(x) ((x(1)^2-x(2))^2+(1-x(1))^2)/2;
iterations=0;
not_min=0;
#hf=figure();
#hold on;
tic()
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
        xk-=grad/laplacian(xk);           
        j++;
    until (normgrad<=eps || j==10000)
    if (normgrad<=eps)
        #scatter(xk_0(1), xk_0(2), "r");
        iterations+=j;
    else
        #scatter(xk_0(1), xk_0(2), "b");    
        not_min++;
    endif    
endfor
toc()
iterations
not_min
#print (hf, "newton.png");
#hold off;