hf=figure;
zlabel("\\frac{||x-z||}{||x||}, ||A^{-1}||||E||, 2\\kappa(A)(\\epsilon_A+\\epsilon_b)");
xlabel("\\epsilon");
hold on;
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
alfa=-0.8
beta=0.6
epsi=[0.2; -0.1; 0.1]
A(:,1)=alfa*A(:,2)+beta*A(:,3)+epsi

norm_A=norm(A)
norm_A_inv=1/norm(inv(A))
cond_A=cond(A)

E=zeros(3, 3);

norm_E=norm(E);

vects=eye(3)
epsiAux=0.1*epsi
for i=1:10
  figure
  quiver3(0, 0, 0, 0, 0, 0)
  hold on;
  for i=1:3
    vect=vects(:,i);
    h=quiver3(0, 0, 0, vect(1), vect(2), vect(3), '--');
    set (h, "color", [(sin(vects(1,i))+1)/2, (sin(vects(2,i)+pi)+1)/2, (sin(vects(3,i)-pi/2)+1)/2]);
    Aux=A;
    Aux(:,1)-=epsiAux;
    vect=Aux*vect;
    h=quiver3(0, 0, 0, vect(1), vect(2), vect(3));
    set (h, "maxheadsize", 0.1/norm(vect));
    set (h, "color", [(sin(vects(1,i))+1)/2, (sin(vects(2,i)+pi)+1)/2, (sin(vects(3,i)-pi/2)+1)/2]);
  end
  epsiAux+=0.1*epsi
  hold off;
end

norm(epsi)