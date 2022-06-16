clear all;

A_minus=[13.43538491431723,  -18.14977774439898; 10.01439047867205,  -13.50150142642550];
A_plus=[3/4, -1; 1, 3/4];

#Se for no sentido horario, colocar t_space decrescente (<=0)
data=cell(2, 5);
cut=1;

#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar 
#pelo corte multiplas vezes entre dois passos)
t_space=linspace(0, 18, 1800);
n_space=-t_space;

[data{1, 1}, data{1, 2}]=eig(A_minus);

#b_minus
data{1, 4}=[0; 0];
data{1, 3}=A_minus\data{1, 4};
#b_minus*A_minus^(-1)

[data{2, 1}, data{2, 2}]=eig(A_plus);

#b_plus
data{2, 4}=[0; 0];
#d=b_plus*A_plus^(-1)
data{2, 3}=A_plus\data{2, 4};

function [Atil_1, Atil_2]=find_elementary(A_minus, A_plus, data, cut)
    Atil_1=-(A_minus(1, 1)*cut+data{1, 4}(1))/A_minus(1, 2);
    Atil_2=-(A_plus(1, 1)*cut+data{2, 4}(1))/A_plus(1, 2);
endfunction

function [Asin_1, Asin_2]=find_singularity(A_minus, A_plus, data)
    Asin_1=transpose(-A_minus\data{1,4});
    Asin_2=transpose(-A_plus\data{2,4});
endfunction

function Y=half_poinc(i, Y, t_space, data, A_minus, cut)
    alpha_minus=real(data{i, 2}(1,1));
    beta_minus=imag(data{i, 2}(1,1));
    d=data{i,3};
    Y_p_x=Y+d(2);
    cut_p_x=cut+d(1);
    if (i==2)
            y=@(t) exp(alpha_minus*t)*(sin(beta_minus*t)*cut_p_x+Y_p_x*cos(beta_minus*t))-d(1);
            x=@(t) exp(alpha_minus*t)*(cos(beta_minus*t)*cut_p_x-Y_p_x*sin(beta_minus*t))-cut_p_x;
    else
            gamma_minus=alpha_minus/beta_minus;
            A_m_diag_dif=A_minus(1,1)-A_minus(2,2);
            x=@(t) exp(alpha_minus*t)*(cut_p_x*(cos(beta_minus*t)+A_m_diag_dif*sin(beta_minus*t)/(2*beta_minus))+Y_p_x*A_minus(1,2)*sin(beta_minus*t)/beta_minus)-cut_p_x;
            y=@(t) exp(alpha_minus*t)*(Y_p_x*(cos(beta_minus*t)-A_m_diag_dif*sin(beta_minus*t)/(2*beta_minus))+cut_p_x*A_minus(2,1)*sin(beta_minus*t)/beta_minus)-d(1);
    endif
    
    #Busca da inversao
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
              if (x(t)<0 && (t!=0))
                  t_inv=t;
                  break;
              end
          endfor
    end
    if (t_inv)
      
        #aux e o t onde a inversao ocorre
        aux=fzero(x, [t_inv-t_space(2), t_inv]);
        Y=y(aux);
    endif
endfunction

function cycle_at=step_down(start, finish, t_space, data, A_minus, cut)
    h=(finish-start)/2;
    Y=start;
    [start, finish]
    res=[Y, 0];
    cycle_at=[0; 0];
    k=0;
    do
        for i=1:2
            Y=half_poinc(i, Y, t_space, data, A_minus, cut);
        endfor
        si=sign(Y-res(1,1));
        
        #Se houve inversao de comportamento, buscar no subintervalo
        if (res(1,2)*si<0);
            #Caso deseja-se uma precisao maior (menor), diminuir (aumentar) o valor abaixo
            if (h>10^(-9))
                cycle_at=step_down(res(1,1)-h, res(1,1), t_space, data, A_minus, cut)
            else
                cycle_at=[res(1,1)-h; res(1,1)];
            end
            break;
        end
        res=[res(1,1)+h, si];
        Y=res(1,1);
        k++
    until (k==3)
endfunction

hf=figure;
xlabel ("x_1");
ylabel ("x_2");
hold on;
[data{1, 5}, data{2, 5}]=find_elementary(A_minus, A_plus, data, cut);
scatter(cut, data{1, 5}, "filled");
scatter(cut, data{2, 5}, "filled");
[Asin_1, Asin_2]=find_singularity(A_minus, A_plus, data);
scatter(Asin_1(1), Asin_1(2), "filled");
scatter(Asin_2(1), Asin_2(2), "filled");
#Mudar p/ melhor enquadramento
plot([cut, cut], [-50, 50], "color", "black");

gamma_minus=real(data{1, 2}(1,1))/imag(data{1, 2}(1,1))
gamma_plus=A_plus(1,1)/A_plus(2,1)

#Calculo de y_minus^*
if (sign(gamma_minus)<0)
    if (gamma_plus<data{1, 5})
        Y=half_poinc(1, data{2, 5}, n_space, data, A_minus, cut) 
    else
        Y=half_poinc(1, data{1, 5}, n_space, data, A_minus, cut) 
    end
else
    Y=half_poinc(1, data{1, 5}, t_space, data, A_minus, cut)
    if (Y>gamma_plus)
        Y=half_poinc(1, data{2, 5}, n_space, data, A_minus, cut)        
    else
        Y=data{1, 5}
    endif
end

#Redefinir max se necessario
max=Y+1
#Diminuir o passo se necessario
h=(max-Y)/100;
res=[Y, 0];
cycles_at=[];
while (Y<=max)
    for i=1:2
        Y=half_poinc(i, Y, t_space, data, A_minus, cut);
    endfor
    si=sign(Y-res(1,1));
    
    #Se houve inversao de comportamento, buscar no subintervalo
    if (res(1,2)*si<0);
        cycles_at=[cycles_at, step_down(res(1,1)-h, res(1,1), t_space, data, A_minus, cut)]
    end
    res=[res(1,1)+h, si];
    Y=res(1,1);
endwhile

#aux=half_poinc(2, data{1, 5}, t_space)
#cycles_at=[cycles_at, [cycles_at(1,3); cycles_at(1,2)], [cycles_at(1,3); cycles_at(1,3)+0.5]];
#text(1.1, data{1, 5}, "y_c^-")
#text(1.1, data{2, 5}, "\\gamma^+")
#Y=half_poinc(1, data{1, 5}, t_space)
#scatter(1, Y, "filled")
#text(1.1, Y, "P^-(y_c^-)")
#cycles_at=[cycles_at, [data{1, 5}; data{1, 5}]]
#Y=half_poinc(1, data{2, 5}, n_space)
#scatter(1, Y, "filled")
#text(1.1, Y, "(P^-)^{-1}(\\gamma^+)")
#cycles_at=[cycles_at, [Y; Y]]
cycles_at=[];

for j=1:columns(cycles_at)
    for i=1:2
        Y=(cycles_at(1,j)+cycles_at(2,j))/2
        if (i==2)
            x=@(t) exp(-real(data{i, 2}(1,1))*t)*(cos(-imag(data{i, 2}(1,1))*t)*cut-Y*sin(-imag(data{i, 2}(1,1))*t))-cut;
            y=@(t) exp(-real(data{i, 2}(1,1))*t)*(sin(-imag(data{i, 2}(1,1))*t)*cut+Y*cos(-imag(data{i, 2}(1,1))*t));
        else
            x=@(t) exp(real(data{i, 2}(1,1))*t)*(cut*(cos(imag(data{i, 2}(1,1))*t)+(A_minus(1,1)-A_minus(2,2))*sin(imag(data{i, 2}(1,1))*t)/(2*imag(data{i, 2}(1,1))))+Y*A_minus(1,2)*sin(imag(data{i, 2}(1,1))*t)/imag(data{i, 2}(1,1)))-cut;
            y=@(t) exp(real(data{i, 2}(1,1))*t)*(Y*(cos(imag(data{i, 2}(1,1))*t)-(A_minus(1,1)-A_minus(2,2))*sin(imag(data{i, 2}(1,1))*t)/(2*imag(data{i, 2}(1,1))))+cut*A_minus(2,1)*sin(imag(data{i, 2}(1,1))*t)/imag(data{i, 2}(1,1)));
        endif
        last=[cut; Y];
        t_inv=0;
        if (i==1)
              for t=t_space
                  X=x(t);
                  Y=y(t);
                  if ((X>0) && (t!=0))
                      t_inv=t
                      break;
                  end                  
                  line([last(1), X+cut], [last(2), Y], "color", "black", "linewidth", 1.5)
                  last=[X+cut, Y];
                  if ((X+0.00001>0) && (t!=0) && (Y==cycles_at(1,1)))
                      Y=data{1, 5}
                      break;
                  end                  
              endfor
        else
              for t=t_space
                  X=x(t);
                  Y=y(t);
                  if ((X<0) && (t!=0))
                      t_inv=t
                      break;
                  end                  
                  line([last(1), X+cut], [last(2), Y], "color", "black", "linewidth", 1.5)
                  last=[X+cut, Y];
              endfor
        end
        if (t_inv)
            aux=fzero(x, [t_inv-t_space(2), t_inv]);
            x(aux)
            Y=y(aux)
            line([last(1), cut], [last(2), Y], "color", "black", "linewidth", 1.5);
        endif
    endfor
endfor
ylim([0.5,1])
xlim([0.99,1.015])
print(hf, "triz.png");
hold off;
