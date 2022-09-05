b=[1; 2; 3];
norm_b=norm(b)
hf=figure;
zlabel("\\frac{||x-z||}{||x||}, ||A^{-1}||||E||, 2\\kappa(A)(\\epsilon_A+\\epsilon_b)");
xlabel("\\epsilon");
hold on;
for i=1:1
  epsi=0;
  #A=[1, 1; 1, 0.905]
  #A=[0.889102, 0.864793; 0.210131, 0.019724]
  #A=rand(2, 2)
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

  norm_A=norm(A)
  norm_A_inv=norm(inv(A))
  cond_A=cond(A)
  
  E=zeros(3, 3);
  E(3, 3)=epsi;
  f=[0; 0; epsi];
  norm_E=norm(E);
  norm_f=norm(f);

  e_A=norm_E/norm_A;
  e_b=norm_f/norm_b;

  restr_last=norm_E*norm_A_inv;
  limit_last=2*cond_A*(e_A+e_b);

  x=A\b
  z=(A+E)\(b+f)

  err_rel_last=norm(x-z)/norm(x);
  while restr_last<2
    epsi+=0.01;
    E=zeros(3, 3);
    E(3, 3)=epsi;
    f=[0; 0; epsi];
    norm_E=norm(E);
    norm_f=norm(f);

    e_A=norm_E/norm_A;
    e_b=norm_f/norm_b;

    restr=norm_E*norm_A_inv;
    limit=2*cond_A*(e_A+e_b);

    z=(A+E)\(b+f);

    err_rel=norm(x-z)/norm(x);
    
    plot([epsi-0.01, epsi], [restr_last, restr], "r")
    plot([epsi-0.01, epsi], [err_rel_last, err_rel], "g")
    plot([epsi-0.01, epsi], [limit_last, limit], "b")
    
    restr_last=restr;
    err_rel_last=err_rel;
    limit_last=limit;
  end
end
plot([0, epsi], [1/2, 1/2], "k")
hold off;
legend("||A^{-1}||||E||", "||x-z||/||x||", "2\\kappa(A)(\\epsilon_A+\\epsilon_b)")
fig=figure;
quiver3(0, 0, 0, 0, 0, 0)
hold on;
[x, y, z] = sphere(20);
for i=1:21
  for j=1:21
    vect=A*[x(j,i); y(j,i); z(j,i)];
    x(j,i)=vect(1);
    y(j,i)=vect(2);
    z(j,i)=vect(3);
   end
end
surf(x, y, z);
hold off;
##i=0;
##while i<pi
##  j=0;
##  while j<2*pi 
##    if !((i==pi/2 && j==3*pi/2) || (i==0 && j==pi/2))
##      i, j
##      v=i+pi*j/2
##      if i==0 && j==0
##        v=2
##      endif
##      vect=[cos(i)*cos(j); sin(i)*cos(j); sin(j)];
##      h=quiver3(0, 0, 0, vect(1), vect(2), vect(3));
##      set (h, "maxheadsize", 0.1);
##      set (h, "color", [(sin(v)+1)/2, (sin(v+pi/2)+1)/2, (sin(v-pi/2)+1)/2]);
##    end
##    j+=pi/2;
##   endwhile
##   i+=pi/2;
##end
##fig=figure;
##quiver3(0, 0, 0, 0, 0, 0)
##hold on;
##i=0;
##while i<pi
##  j=0;
##  while j<2*pi 
##    if !((i==pi/2 && j==3*pi/2) || (i==0 && j==pi/2))
##      i, j
##      v=i+pi*j/2
##      if i==0 && j==0
##        v=2
##      endif
##      vect=[cos(i)*cos(j); sin(i)*cos(j); sin(j)];
##      vect=A*vect;
##      h=quiver3(0, 0, 0, vect(1), vect(2), vect(3));
##      set (h, "maxheadsize", 0.1/norm(vect));
##      set (h, "color", [(sin(v)+1)/2, (sin(v+pi/2)+1)/2, (sin(v-pi/2)+1)/2]);
##    end
##    j+=pi/2;
##   endwhile
##   i+=pi/2;
##end
aux=int2str(cond_A*1000)
print(hf, strcat(aux, "graph.png"));
print(fig, strcat(aux, "deform.png"));
hold off;
strrep(strrep(mat2str(A, 4),",","\&"),";","\\\\\n")(2:end-1)