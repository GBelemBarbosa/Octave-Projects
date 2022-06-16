A=[1,1,1,0;2,3,0,-1];
b=[4;18];
custos=[-3;4;0;0];

function sol=simplex(A, b, custos_original)
    iB=[];
    iN=[];
    where_y=linspace(1, rows(A), rows(A));
    i=0;
    for c=transpose(custos_original)
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
    custos=zeros(columns(A), 1);
    j=columns(A);
    for c=where_y
        j++;
        aux=zeros(rows(A), 1);
        aux(c)=1;
        A=[A, aux];
        custos=[custos; 1] ;
        iB=[iB, j];
    endfor
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
    Xb=B\b
    for k=1:2
        while (1)
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
              Xb=B\b
        end
        Xb
        
        #Fase 2
        custos=custos_original;
        iN(iN>rows(custos))=[]
        A=A(:, 1:rows(custos));
        iB=transpose(sol(end-rows(A)+1:end, end))
        B=[];
        Cb=[];
        for i=iB
            B=[B, A(:, i)];
            Cb=[Cb; custos(i)];
        endfor
    endfor
endfunction

sol=simplex(A, b, custos);
