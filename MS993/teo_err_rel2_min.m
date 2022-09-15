A=rand(3, 3);
for j=1:3
  for p=1:3
    if rand(1)>0.5
      A(j, p)*=-1;
    end
    if rand(1)>0.5
      A_copy=A;
      A_copy(j, p)=0;
      if cond(A_copy)<10^10
        A(j, p)=0;
      end
    end
  endfor
endfor
A
inv_norm_inv=1/norm(inv(A))
cond_A=cond(A)

function weigth=minimize(A, e)
  E=[e(1:3);e(4:6);e(7:9)];
  A+=E;
  v=A(:,1).'*A(:,2)*A(:,2)/norm(A(:,2))+A(:,1).'*A(:,3)*A(:,3)/norm(A(:,3));
  weigth=norm(A(:,1)-v)+norm(E);
end

aux=zeros(1,9);
aux(1)=0.001;

e=fminunc(@(e) minimize(A, e), aux);
E=[e(1:3);e(4:6);e(7:9)]
cond_AE_col=cond(A+E)
A_E_col=A+E
norm_E_col=norm(E)

e=fminunc(@(e) minimize(A.', e), aux);
E=[e(1:3);e(4:6);e(7:9)];
E=E.'
cond_AE_row=cond(A+E)
A_E_row=A+E
norm_E_row=norm(E)
