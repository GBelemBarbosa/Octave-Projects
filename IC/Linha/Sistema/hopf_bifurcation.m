clear all;

#A_1 e sempre A-
A_1=[1, 4; -1, 1];
A_2=[-1, 2; -1, -1];

#Se for no sentido horario, colocar t_space decrescente (<=0)
global matrices=cell(2, 5) e cut=0;
#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar pelo corte multiplas vezes entre dois passos)
t_space=linspace(0, -18, 1200);
n_space=-t_space;

[matrices{1, 1}, matrices{1, 2}]=eig(A_1);
matrices{1, 4}=[-0.2; 0.05];
matrices{1, 3}=A_1\matrices{1, 4};

[matrices{2, 1}, matrices{2, 2}]=eig(A_2);
matrices{2, 4}=[-0.1; -0.1];
matrices{2, 3}=A_2\matrices{2, 4};

function [Atil_1, Atil_2]=find_elementary(A_1, A_2)
    global matrices e cut;
    Atil_1=-(A_1(1, 1)*cut+matrices{1, 4}(1))/A_1(1, 2);
    Atil_2=-(A_2(1, 1)*cut+matrices{2, 4}(1))/A_2(1, 2);
endfunction

function [Asin_1, Asin_2]=find_singularity(A_1, A_2)
    global matrices;
    Asin_1=transpose(-A_1\matrices{1,4});
    Asin_2=transpose(-A_2\matrices{2,4});
endfunction

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

function cycle_at=step_down(start, finish, t_space)
    h=(finish-start)/20;
    Y=start;
    [start, finish]
    res=[Y, 0];
    cycle_at=[0; 0];
    while (Y<=finish)
        for i=1:2
            Y=find_Y(i, Y, t_space);
        endfor
        si=sign(Y-res(1,1));
        if (res(1,2)*si<0);
            #matricesaso deseja-se uma precisao maior (menor), diminuir (aumentar) o valor abaixo
            if (h>10^(-4))
                cycle_at=step_down(res(1,1)-h, res(1,1), t_space);
            else
                cycle_at=[res(1,1)-h; res(1,1)];
            end
            break;
        end
        res=[res(1,1)+h, si];
        Y=res(1,1);
    endwhile
endfunction

hf=figure;
xlabel ("x_1");
ylabel ("x_2");
hold on;
[matrices{1, 5}, matrices{2, 5}]=find_elementary(A_1, A_2);
scatter(cut, matrices{1, 5}, "filled");
scatter(cut, matrices{2, 5}, "filled");
[Asin_1, Asin_2]=find_singularity(A_1, A_2);
scatter(Asin_1(1), Asin_1(2), "filled");
scatter(Asin_2(1), Asin_2(2), "filled");
#Mudar p/ melhor enquadramento
plot([cut, cut], [-50, 50], "color", "black");

gamma_minus=real(matrices{1, 2}(1,1))/imag(matrices{1, 2}(1,1))
gamma_plus=A_2(1,1)/A_2(2,1)
if (sign(gamma_minus)<0)
    if (gamma_plus<matrices{1, 5})
        Y=find_Y(1, matrices{2, 5}, n_space) 
    else
        Y=find_Y(1, matrices{1, 5}, n_space) 
    end
else
    Y=find_Y(1, matrices{1, 5}, t_space)
    if (Y>gamma_plus)
        Y=find_Y(1, matrices{2, 5}, n_space)        
    else
        Y=matrices{1, 5}
    endif
end
#Redefinir max se necessario(com base no grafico teste)
max=Y+1
#Diminuir o passo se necessario (com base no grafico teste)
h=(max-Y)/30;
res=[Y, 0];
cycles_at=[];
while (Y<=max)
    for i=1:2
        Y=find_Y(i, Y, t_space);
    endfor
    si=sign(Y-res(1,1));
    if (res(1,2)*si<0);
        cycles_at=[cycles_at, step_down(res(1,1)-h, res(1,1), t_space)]
    end
    res=[res(1,1)+h, si];
    Y=res(1,1);
endwhile

#aux=find_Y(2, matrices{1, 5}, t_space)
#cycles_at=[cycles_at, [cycles_at(1,3); cycles_at(1,2)], [cycles_at(1,3); cycles_at(1,3)+0.5]];
#text(1.1, matrices{1, 5}, "y_c^-")
#text(1.1, matrices{2, 5}, "\\gamma^+")
#Y=find_Y(1, matrices{1, 5}, t_space)
#scatter(1, Y, "filled")
#text(1.1, Y, "P^-(y_c^-)")
#cycles_at=[cycles_at, [matrices{1, 5}; matrices{1, 5}]]
#Y=find_Y(1, matrices{2, 5}, n_space)
#scatter(1, Y, "filled")
#text(1.1, Y, "(P^-)^{-1}(\\gamma^+)")
#cycles_at=[cycles_at, [Y; Y]]
cycles_at=[cycles_at, [0.5; 0.5], [0.25; 0.25]]
#cycles_at=[];

for j=1:columns(cycles_at)
    Y=(cycles_at(1,j)+cycles_at(2,j))/2;
    for i=1:2
        x=@(t) real(((matrices{i, 1}*expm(matrices{i, 2}*t)/matrices{i, 1})*([cut; Y]+matrices{i, 3})-matrices{i, 3})(1))-cut;
        xy=@(t) real((matrices{i, 1}*expm(matrices{i, 2}*t)/matrices{i, 1})*([cut; Y]+matrices{i, 3})-matrices{i, 3});
        last=[cut; Y];
        t_inv=0;
        if (i==1)
              for t=t_space
                  aux=xy(t);
                  if ((aux(1)>cut) && (t!=0))
                      t_inv=t;
                      break;
                  end                  
                  line([last(1), aux(1)], [last(2), aux(2)], "color", "black", "linewidth", 1.5)
                  last=aux;
                  if ((aux(1)+0.00001>cut) && (t!=0) && (Y==cycles_at(1,1)))
                      Y=matrices{1, 5}
                      break;
                  end                  
              endfor
        else
              for t=t_space
                  aux=xy(t);
                  if ((aux(1)<cut) && (t!=0))
                      t_inv=t;
                      break;
                  end                  
                  line([last(1), aux(1)], [last(2), aux(2)], "color", "black", "linewidth", 1.5)
                  last=aux;
              endfor
        end
        if (t_inv)
            aux=fzero(x, [t_inv-t_space(2), t_inv]);
            Y=xy(aux)(2);
            line([last(1), cut], [last(2), Y], "color", "black", "linewidth", 1.5);
        endif
    endfor
endfor
ylim([-0.1,0.6])
xlim([-0.5, 0.3])
print (hf, "disturb_mid.png");
hold off;
