clear all;

A_2=[3/4, -1; 1, 3/4];
global  cut=1 e gamma_plus=-A_2(1, 1)/A_2(1, 2) e M=[];
#gamma_plus=-cut*A_2(1, 1)/A_2(1, 2)
#A_1=(M*J)/M
A_1=[13.43538491431723,  -18.14977774439898; 10.01439047867205,  -13.50150142642550];
[M, J]=eig(A_1)

#Se for no sentido horario, colocar t_space decrescente (<=0)
#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar pelo corte multiplas vezes entre dois passos)
global matrices=cell(2, 4) e t_space=linspace(0, 12, 1200) e n_space=-t_space;

matrices{1, 1}=M;
matrices{1, 2}=J;
matrices{1, 4}=[0; 0];
matrices{1, 3}=A_1\matrices{1, 4};

[matrices{2, 1}, matrices{2, 2}]=eig(A_2);
matrices{2, 4}=[0; 0];
matrices{2, 3}=A_2\matrices{2, 4};

function Y=find_Y(i, Y, t_space)
    global cut e matrices;
    x=@(t) real(((matrices{i, 1}*expm(matrices{i, 2}*t)/matrices{i, 1})*([cut; Y]+matrices{i, 3})-matrices{i, 3})(1))-cut;
    xy=@(t) real((matrices{i, 1}*expm(matrices{i, 2}*t)/matrices{i, 1})*([cut; Y]+matrices{i, 3})-matrices{i, 3});
    t_inv=0;
    if (i==1)
          for t=t_space
              if ((x(t)>0) && (t!=0))
                  t_inv=t;
                  break;
              end
          endfor
    else
          for t=t_space
              if (x(t)<0)
                  t_inv=t;
                  break;
              end
          endfor
    end
    if (t_inv)
        aux=fzero(x, [t_inv-t_space(2), t_inv]);
        Y=xy(aux)(2);
    endif
endfunction

y_0_plus=@(y_1) -(y_1+gamma_plus*(1+exp(gamma_plus*pi)))/exp(gamma_plus*pi);
gamma_minus=real(J(1,1))/imag(J(1,1))
yc_minus=-(A_1(1, 1)*cut+matrices{1, 4}(1))/A_1(1, 2)
y_1_minus=@(y_0) -exp(gamma_minus*pi)*y_0+(A_1(2,2)/A_1(1,2))*(1+exp(gamma_minus*pi));
Y=find_Y(1, yc_minus, n_space)
y_space=linspace(Y, 1, 50);
hf=figure();
hold on;
scatter(Y, 0.72, "filled");
axis([Y, Y+0.05, yc_minus-0.5, yc_minus+0.5])
line([y_space(1), y_space(end)], [y_0_plus(y_space(1)), y_0_plus(y_space(end))], "color", "g", "linestyle", "--");
line([y_space(1), y_space(end)], [y_1_minus(y_space(1)), y_1_minus(y_space(end))], "color", "b", "linestyle", "--");
last=[yc_minus, find_Y(2, Y, n_space)];
(y_1_minus(y_space(end))-y_1_minus(y_space(1)))/1.2
h=y_space(2)-y_space(1); 
for y=y_space(2:end)
    Y=find_Y(1, y, t_space);
    line([y, y-h], [Y, last(1)], "color", "b");
    last(1)=Y;
    Y=find_Y(2, y, n_space);
    line([y, y-h], [Y, last(2)], "color", "g");
    last(2)=Y;
endfor
legend ({"C", "Assíntota de (P^+)^{-1}", "Assíntota de P^-", "P^-", "(P^+)^{-1}"}, "location", "northeast");
ylim([0.4,0.8])
xlim("auto")
print(hf, "three.png");    
