A=[2  ,-1  ,-1 , -1,   0 ,  1;
   1 , -1  , 3   ,0 , -1  , 0]
b=[-3;-2];
custos=[6;4;-4;0;0];

function print_sol(sol, hf)
    points=[];
    o=rows(sol);
    for i=1:columns(sol)
        if any(sol(o/2+1:o, i)==1)
            j=find(sol(o/2+1:o, i)==1)
            if any(sol(o/2+1:o, i)==2)
                p=find(sol(o/2+1:o, i)==2)
                scatter(sol(j, i), sol(p, i), "filled");
                points=[points, [sol(j, i); sol(p, i)]]
                text(sol(j, i)+0.1, sol(p, i)+0.1, strcat("X", num2str(i-1)));
            else
                scatter(sol(j, i), 0, "filled");
                points=[points, [sol(j, i); 0]]
                text(sol(j, i)+0.1, 0.1, strcat("X", num2str(i-1)));
            endif
        elseif any(sol(o/2+1:o, i)==2)
            j=find(sol(o/2+1:o, i)==2)
            scatter(0, sol(j, i), "filled");
            points=[points, [0; sol(j, i)]]
            text(0.1, sol(j, i)+0.1, strcat("X", num2str(i-1)));
        else
            scatter(0, 0, "filled");
            points=[points, [0; 0]]
            text(0.1, 0.1, strcat("X", num2str(i-1)));
        endif
    endfor
    points;
    for i=1:columns(points)-1
        quiver(points(1, i), points(2, i), points(1, i+1)-points(1, i), points(2, i+1)-points(2, i));
    endfor
endfunction

function sol=Big_M(A, custos, b)
    iB=[];
    iN=[];
    where_y=linspace(1, rows(A), rows(A));
    M=100*max(abs(custos));
    i=0;
    for c=transpose(custos)
        i++;
        if (c==0)
            j=find(A(:, i)!=0);
            if (sign(A(j, i))*sign(b(j))>=0)
                where_y(find(where_y==j))=[];
                iB=[iB, i];
            else
                iN=[iN, i];
            endif
        else
            iN=[iN, i];
        endif
    endfor
    where_y
    j=columns(A);
    for c=where_y
        j++;
        aux=zeros(rows(A), 1);
        aux(c)=1;
        A=[A, aux];
        custos=[custos; M] ;
        iB=[iB, j];
    endfor
    A
    custos
    iB
    iN 
    
    sol=[];
    B=[];
    Cb=[];
    for i=iB
        B=[B, A(:, i)];
        Cb=[Cb; custos(i)];
    endfor
    while (1)
          Xb=B\b
          transpose(Cb)*Xb
          sol=[sol, [Xb; transpose(iB)]];
          
          lambda_T=transpose(transpose(B)\Cb)
          custos_redux=[];     
          for i=iN
              custos_redux=[custos_redux, custos(i)-lambda_T*A(:, i)];
          endfor
          custos_redux
          [minimo_redux, i]=min(custos_redux)
          if (minimo_redux>=0)
              break;
          endif
          
          y=B\A(:, iN(i))
          epsilons=[];
          j=0;
          for k=iB
              j++;
              if (y(j)>0)
                  epsilons=[epsilons, [Xb(j)/y(j); j]];
              endif
          endfor
          epsilons
          if (!isempty(epsilons))
              [min_epsilon, p]=min(epsilons(1,:))
              p=epsilons(2, p)
              j=iN(i)
              iN(i)=iB(p)
              iB(p)=j
              B(:, p)=A(:, i);
              Cb(p)=custos(i);
          else
              display("Sem solucao limitada.");
              break;
          endif
    end
    iB
    iN
    Xb
    if any(iB>columns(A))
        display("PL infactivel.");
    endif
    
endfunction

sol=Big_M(A, custos, b)

hf=figure;
xlabel("x_1");
ylabel("x_2");
hold on;
print_sol(sol, hf);
for i=1:rows(A)
    if (A(i,2)!=0)
        fplot(@(x) (b(i)-A(i,1)*x)/A(i,2),[0,5.5]);
    end
endfor
line([0.01, 0.01], [-4, 8]);
line([0, 5.5], [0, 0]);
legend("hide");
axis([-0.5, 5.5]);  
print (hf, "c.png");

hold off;