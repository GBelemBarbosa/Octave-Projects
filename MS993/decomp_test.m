pkg install -forge benchmark

num_iter=1;
n=3;
b=[]
for i=1:n
  b=[b; i]
end

function [Q, R]=qrgivens(A)
  [m, n]=size(A);
  Q=eye(m);
  R=A;

  for j=1:n
    for i=m:-1:(j+1)
      G=givens(R(i-1,j),R(i,j));
      R=G'*R;
      Q=Q*G;
    end
  end
end

function [Q, R]=GS(A)
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
end

function [Q, R]=GSM(A)
  for j=1:n
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j);
        A(:,j)-=R(i,j)*Q(:,i);
    end
    R(j,j)=norm(A(:,j));
    if (R(j,j)==0)
        break
    end
    Q(:,j)=A(:,j)/R(j,j);
  end
end

mse=@(A, B) sum((A.-B).^2)/(size(A, 1)*size(A, 2))
mabse=@(A, B) max(abs.(A.-B))

A_mse_m_chol=0;
A_mse_m_house=0;
A_mse_m_givens=0;
A_mse_m_GS=0;
A_mse_m_GSM=0;
A_mse_m_LU=0;
A_mse_m_LUP=0;
A_mse_m_LUPQ=0;

A_mabse_m_chol=0;
A_mabse_m_house=0;
A_mabse_m_givens=0;
A_mabse_m_GS=0;
A_mabse_m_GSM=0;
A_mabse_m_LU=0;
A_mabse_m_LUP=0;
A_mabse_m_LUPQ=0;

b_mse_m_chol=0;
b_mse_m_house=0;
b_mse_m_givens=0;
b_mse_m_GS=0;
b_mse_m_GSM=0;
b_mse_m_LU=0;
b_mse_m_LUP=0;
b_mse_m_LUPQ=0;

for i=1:num_iter
  A=rand(n, n);
  
  #Cholesky
  R=chol(A'*A);
  Q=R\A;
  x=R\(Q'*b);
  
  Aux=Q*R;
  A_mse_m_chol+=mse(A, Aux);
  A_mabse_m_chol+=mabse(A, Aux);
  b_mse_m_chol+=mse(b, A*x);
  
  #Householder
  [Q, R]=qr(A);
  x=R\(Q'*b);
  
  Aux=Q*R;
  A_mse_m_house+=mse(A, Aux);
  A_mabse_m_house+=mabse(A, Aux);
  b_mse_m_house+=mse(b, A*x);
  
  #Givens
  [Q, R]=qrgivens(A);
  x=R\(Q'*b);
  
  Aux=Q*R;
  A_mse_m_givens+=mse(A, Aux);
  A_mabse_m_givens+=mabse(A, Aux);
  b_mse_m_givens+=mse(b, A*x);
  
  #GS
  [Q, R]=GS(A);
  x=R\(Q'*b);
  
  Aux=Q*R;
  A_mse_m_GS+=mse(A, Aux);
  A_mabse_m_GS+=mabse(A, Aux);
  b_mse_m_GS+=mse(b, A*x);
  
  #GSM
  [Q, R]=GSM(A);
  x=R\(Q'*b);
  
  Aux=Q*R;
  A_mse_m_GSM+=mse(A, Aux);
  A_mabse_m_GSM+=mabse(A, Aux);
  b_mse_m_GSM+=mse(b, A*x);
  
  #LU
  [L, U]=lu(A);
  x=U\(L\b);
  
  Aux=L*U;
  A_mse_m_LU+=mse(A, Aux);
  A_mabse_m_LU+=mabse(A, Aux);
  b_mse_m_LU+=mse(b, A*x);
  
  #LUP
  [L, U, P]=lu(A);
  x=U\(L\(P*b));
  
  Aux=P'*L*U;
  A_mse_m_LUP+=mse(A, Aux);
  A_mabse_m_LUP+=mabse(A, Aux);
  b_mse_m_LUP+=mse(b, A*x);
  
  #LUPQ
  [L, U, P]=lu(A);
  x=Q*(U\(L\(P*b)));
  
  Aux=P'*L*U*Q';
  A_mse_m_LUPQ+=mse(A, Aux);
  A_mabse_m_LUPQ+=mabse(A, Aux);
  b_mse_m_LUPQ+=mse(b, A*x);
end

A_mse_m_chol/=num_iter
A_mse_m_house/=num_iter
A_mse_m_givens/=num_iter
A_mse_m_GS/=num_iter
A_mse_m_GSM/=num_iter
A_mse_m_LU/=num_iter
A_mse_m_LUP/=num_iter
A_mse_m_LUPQ/=num_iter

b_mse_m_chol/=num_iter
b_mse_m_house/=num_iter
b_mse_m_givens/=num_iter
b_mse_m_GS/=num_iter
b_mse_m_GSM/=num_iter
b_mse_m_LU/=num_iter
b_mse_m_LUP/=num_iter
b_mse_m_LUPQ/=num_iter