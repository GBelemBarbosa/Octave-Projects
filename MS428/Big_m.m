#Inicializacao e entrada
A=[];
b=[];
custos=[];
m=input("Número de linhas: ");
n=input("Número de colunas: ");
for i=1:m
    aux=[];
    for j=1:n
        v=input(strcat("Elemento a_",num2str(i),num2str(j),":"));
        aux=[aux, v];
    endfor
    A=[A; aux];
endfor
for i=1:m
    v=input(strcat("Elemento b_",num2str(i),":"));
    b=[b; v];
endfor
for i=1:n
    v=input(strcat("Elemento c_",num2str(i),":"));
    custos=[custos; v];
endfor
function sol=Big_M(A, custos, b)
  
    #Forma Big-M e base inicial
    iB=[];
    iN=[];
    m=rows(A);
    n=columns(A);
    where_y=linspace(1, m, m);
    M=100*max(abs(custos));
    i=0;
    for c=transpose(custos)
        i++;
        if (c==0)
            j=find(A(:, i)!=0);
            if ((rows(j)==1) & (sign(A(j, i))*sign(b(j))>=0))
              j=find(where_y==j);
              if (j)
                where_y(j)=[];
                iB=[iB, i];
              else 
                iN=[iN,i];
              endif             
          else
              iN=[iN, i];
          endif
        else
            iN=[iN, i];
        endif
    endfor    
    j=n;
    for i=where_y
        j++;
        aux=zeros(m, 1);
        if (b(i)>=0)
            aux(i)=1;
        else
            aux(i)=-1;
        endif
        A=[A, aux];
        custos=[custos; M] ;
        iB=[iB, j];
    endfor
    
    #Matriz A na forma Big-M
    A
    iB
    sol=[];
    B=[];
    Cb=[];
    limitada=1;
    for i=iB
        B=[B, A(:, i)];
        Cb=[Cb; custos(i)];
    endfor
    
    #Base
    B

    #Metodo simplex	
    while (1)
      
        #Passo 1
        Xb=B\b
        if (any(Xb<0))
            break;
        endif
        
        #f da soluçao basica
        transpose(Cb)*Xb
        sol=[sol, [Xb; transpose(iB)]];
        
        #Passo 2.1
        lambda_T=transpose(transpose(B)\Cb);
        
        #Passo 2.2 e 2.3
        minimo_redux=realmax;
        j=0;     
        for i=iN
            j++;
            custo_redux=custos(i)-lambda_T*A(:, i);
            if (custo_redux<minimo_redux)
                minimo_redux=custo_redux;
                N_index=j;
            endif
        endfor
        in_index=iN(N_index)
        
        #Passo 3
        if (minimo_redux>=0)
            break;
        endif
        
        #Passo 4
        y=B\A(:, in_index);
        
        #Passo 5
        min_epsilon=realmax;
        B_index=0;
        for j=1:m
            if (y(j)>0)
                epsilon=Xb(j)/y(j);
                if (epsilon<=min_epsilon)
                    B_index=j;
                    min_epsilon=epsilon;
                endif
            endif
        endfor
        
        if (B_index)
            out_index=iB(B_index)
            iN(N_index)=out_index;
            iB(B_index)=in_index
            B(:, B_index)=A(:, in_index)
            Cb(B_index)=custos(in_index);
        else
            display("Sem solucao limitada.");
            limitada=0;
            break;
        endif
    endwhile
    
    #Analise da factibilidade do PL original
    if any(iB>n)
        display("PL infactivel.");
    elseif(limitada)
        display("Solucao otima (com indices correspondentes) e valor otimo.")
        Xb
        iB
        transpose(Cb)*Xb
    endif
endfunction

sol=Big_M(A, custos, b);