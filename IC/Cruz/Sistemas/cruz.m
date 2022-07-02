clear all;

A_1=[0,  1; -1,  0];
A_2=[0,  1; -1,  0];
A_3=[0,  1; -1,  0];
A_4=[0,  1; -1,  0];

#Se for no sentido horario, colocar t_space decrescente (<=0)
data=cell(4, 4);

#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar 
#pelo corte multiplas vezes entre dois passos)
t_space=linspace(0, 1, 500);
n_space=-t_space;

[data{1, 1}, data{1, 2}]=eig(A_1);

#b_1
data{1, 4}=[2; 1];
#d=b_1*A_2^(-1)
data{1, 3}=A_1\data{1, 4};

[data{2, 1}, data{2, 2}]=eig(A_2);

#b_2
data{2, 4}=[3; 4];
#d=b_2*A_2^(-1)
data{2, 3}=A_2\data{2, 4};

[data{3, 1}, data{3, 2}]=eig(A_3);

#b_3
data{3, 4}=[1; -1];
#d=b_3*A_3^(-1)
data{3, 3}=A_3\data{3, 4};

[data{4, 1}, data{4, 2}]=eig(A_4);

#b_4
data{4, 4}=[-2; 2];
#d=b_4*A_4^(-1)
data{4, 3}=A_4\data{4, 4};

function tangents=find_elementary(A_1, A_2, A_3, A_4, data)
    Atil_1_y=-data{1, 4}(1)/A_1(1, 2);
    Atil_1_x=-data{1, 4}(2)/A_1(2, 1);
    Atil_2_y=-data{2, 4}(1)/A_2(1, 2);
    Atil_2_x=-data{2, 4}(2)/A_2(2, 1);
    Atil_3_y=-data{3, 4}(1)/A_3(1, 2);
    Atil_3_x=-data{3, 4}(2)/A_3(2, 1);
    Atil_4_y=-data{4, 4}(1)/A_4(1, 2);
    Atil_4_x=-data{4, 4}(2)/A_4(2, 1);
    tangents=[Atil_1_y, Atil_1_x, Atil_2_y, Atil_2_x, Atil_3_y, Atil_3_x, Atil_4_y, Atil_4_x]
endfunction

function singularities=find_singularity(A_1, A_2, A_3, A_4, data)
    Asin_1=transpose(-A_1\data{1,4});
    Asin_2=transpose(-A_2\data{2,4});
    Asin_3=transpose(-A_3\data{3,4});
    Asin_4=transpose(-A_4\data{4,4});
    singularities=[Asin_1; Asin_2; Asin_3; Asin_4]
endfunction

function point_res=quarter_poinc(i, point, t_space, data, A)
    alpha=real(data{i, 2}(1,1));
    beta=imag(data{i, 2}(1,1));
    d=data{i,3};
    gamma=alpha/beta;
    A_m_diag_dif=A(1,1)-A(2,2);
    if !mod(i,2)
            point_p_b=point+d(2);
            x=@(t) exp(alpha*t)*(d(1)*(cos(beta*t)+A_m_diag_dif*sin(beta*t)/(2*beta))+point_p_b*A(1,2)*sin(beta*t)/beta)-d(1);
            y=@(t) exp(alpha*t)*(point_p_b*(cos(beta*t)-A_m_diag_dif*sin(beta*t)/(2*beta))+d(1)*A(2,1)*sin(beta*t)/beta)-d(2);
    else
            point_p_b=point+d(1);
            x=@(t) exp(alpha*t)*(point_p_b*(cos(beta*t)+A_m_diag_dif*sin(beta*t)/(2*beta))+d(2)*A(1,2)*sin(beta*t)/beta)-d(1);
            y=@(t) exp(alpha*t)*(d(2)*(cos(beta*t)-A_m_diag_dif*sin(beta*t)/(2*beta))+point_p_b*A(2,1)*sin(beta*t)/beta)-d(2);
    endif
    
    #Busca da inversao
    t_inv=0;
    if !mod(i, 2) && (point*(i-3)>0)
          last=[0, point]
          for t=t_space
              Y=y(t);
              X=x(t);
              line([last(1), X], [last(2), Y], "color", "black", "linewidth", 1.5)
              last=[X, Y];
              if ((Y*(i-3)<0) && (t!=0))
                  t_inv=t
                  break;
              end
          endfor
    elseif (point*(i-2)<0)
          last=[point, 0]
          for t=t_space
              Y=y(t);
              X=x(t);
              line([last(1), X], [last(2), Y], "color", "black", "linewidth", 1.5)
              last=[X, Y];
              if (X*(i-2)>-0.00001 && (t!=0))
                  t_inv=t
                  break;
              end             
          endfor
    end
    if (t_inv)
      
        #aux e o t onde a inversao ocorre
        if !mod(i,2)
            aux=fzero(y, [t_inv-t_space(2), t_inv]);
            point_res=x(aux);
        else
            aux=fzero(x, [t_inv-t_space(2), t_inv]);
            point_res=y(aux);
        end
    endif
endfunction

hf=figure;
xlabel ("x_1");
ylabel ("x_2");
hold on;
plot([0, 0], [-10, 10], "color", "black");
plot([-10, 10], [0, 0], "color", "black");

tangents=find_elementary(A_1, A_2, A_3, A_4, data);
for i=1:length(tangents)
    if mod(i, 2)
        scatter(0, tangents(i), "filled");
        text(0, tangents(i), strcat('Ã_', num2str(round(i/2)), '^y'));
        if i==3
            quarter_poinc(2, tangents(i), n_space, data, A_2)
        elseif i==7
            quarter_poinc(4, tangents(i), n_space, data, A_4)
        endif
    else
        scatter(tangents(i), 0, "filled");
        text(tangents(i), 0, strcat('Ã_', num2str(round(i/2)), '^x'));
        if i==2
            quarter_poinc(1, tangents(i), n_space, data, A_1)
        elseif i==6
            quarter_poinc(3, tangents(i), n_space, data, A_3)
        endif
    endif
endfor

singularities=find_singularity(A_1, A_2, A_3, A_4, data);
for i=1:length(singularities)
    scatter(singularities(i,1), singularities(i,2), "filled");
    text(singularities(i,1), singularities(i,2), strcat('A_', num2str(i)));
endfor