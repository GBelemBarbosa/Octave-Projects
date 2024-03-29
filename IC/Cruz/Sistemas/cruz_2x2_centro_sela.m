a=-1;

A_1=[0,  -1/a; a,  0];
M=[1,  1; sqrt(3),  1/sqrt(3)];
J=[1,  0; 0,  -1];
A_2=M*J/M;

#Se for no sentido horario, colocar t_space decrescente (<=0)
data=cell(2, 4);

#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar 
#pelo corte multiplas vezes entre dois passos)
t_space=linspace(0, 7, 700);
n_space=-t_space;

[data{1, 1}, data{1, 2}]=eig(A_1);

x=1;
y=1;

#b_1
data{1, 4}=[0.3; 0.3];
data{1, 3}=A_1\data{1, 4};
#b_1*A_1^(-1)

data{2, 1}=M;
data{2, 2}=J;

#b_2
data{2, 4}=[0.5; -0.5];
#d=b_2*A_2^(-1)
data{2, 3}=A_2\data{2, 4};

function tangents=find_elementary(A_1, A_2, data)
    Atil_1_y=-data{1, 4}(1)/A_1(1, 2);
    Atil_1_x=-data{1, 4}(2)/A_1(2, 1);
    tangents=[Atil_1_y, Atil_1_x]
endfunction

function singularities=find_singularity(A_1, A_2, data)
    Asin_1=transpose(-A_1\data{1,4});
    Asin_2=transpose(-A_2\data{2,4});
    singularities=[Asin_1; Asin_2]
endfunction

function point_res=quarter_poinc2(i, point, t_space, data, A)    
    d=data{2-mod(i,2),3};
    A_m_diag_dif=A(1,1)-A(2,2);
    if !mod(i,2)
            lambda_1=data{2-mod(i,2), 2}(1,1);
            lambda_2=data{2-mod(i,2), 2}(2,2);
            C=data{2-mod(i,2), 1}\[d(1); point+d(2)];
            x=@(t) [e^(lambda_1*t), e^(lambda_2*t)]*C-d(1);
            y=@(t) [data{2-mod(i,2), 1}(2,1)*e^(lambda_1*t), data{2-mod(i,2), 1}(2,2)*e^(lambda_2*t)]*C-d(2);
    else
            point_p_b=point+d(1);
            x=@(t) point_p_b*(cos(t)+A_m_diag_dif*sin(t)/2)+d(2)*A(1,2)*sin(t)-d(1);
            y=@(t) d(2)*(cos(t)-A_m_diag_dif*sin(t)/2)+point_p_b*A(2,1)*sin(t)-d(2);
    endif
    
    #Busca da inversao
    t_inv=0;
    if !mod(i, 2) && (point*(i-3)>0)
          last=[0, point]
          for t=t_space
              Y=y(t);
              X=x(t);
              if ((Y*(i-3)<0) && (t!=0))
                  t_inv=t
                  break;
              else
                  line([last(1), X], [last(2), Y], "color", "black", "linewidth", 1.5)
                  last=[X, Y];
              end
          endfor
    elseif (point*(i-2)<0)
          last=[point, 0]
          for t=t_space
              Y=y(t);
              X=x(t);              
              if (X*(i-2)>-0.00001 && (t!=0))
                  t_inv=t
                  break;
              else 
                  line([last(1), X], [last(2), Y], "color", "black", "linewidth", 1.5)
                  last=[X, Y];
              end             
          endfor
    end
    if (t_inv)
      
        #aux e o t onde a inversao ocorre
        try
            if !mod(i,2)
                aux=fzero(y, [t_inv-t_space(2), t_inv]);
                point_res=x(aux);
            else
                aux=fzero(x, [t_inv-t_space(2), t_inv]);
                point_res=y(aux);
            end           
            line([last(1), x(aux)], [last(2), y(aux)], "color", "black", "linewidth", 1.5)
        catch
            point_res=0;
        end_try_catch
    endif
endfunction

function point_res=quarter_poinc2_solver(i, point, t_space, data, A)    
    d=data{2-mod(i,2),3};
    A_m_diag_dif=A(1,1)-A(2,2);
    if !mod(i,2)
            lambda_1=data{2-mod(i,2), 2}(1,1);
            lambda_2=data{2-mod(i,2), 2}(2,2);
            C=data{2-mod(i,2), 1}\[d(1); point+d(2)];
            x=@(t) [e^(lambda_1*t), e^(lambda_2*t)]*C-d(1);
            y=@(t) [data{2-mod(i,2), 1}(2,1)*e^(lambda_1*t), data{2-mod(i,2), 1}(2,2)*e^(lambda_2*t)]*C-d(2);
    else
            point_p_b=point+d(1);
            x=@(t) point_p_b*(cos(t)+A_m_diag_dif*sin(t)/2)+d(2)*A(1,2)*sin(t)-d(1);
            y=@(t) d(2)*(cos(t)-A_m_diag_dif*sin(t)/2)+point_p_b*A(2,1)*sin(t)-d(2);
    endif
    
    #Busca da inversao
    t_inv=0;
    if !mod(i, 2) && (point*(i-3)>0)
          for t=t_space
              Y=y(t);
              if ((Y*(i-3)<0) && (t!=0))
                  t_inv=t
                  break;              
              end
          endfor
    elseif (point*(i-2)<0)
          for t=t_space
              X=x(t);              
              if (X*(i-2)>-0.00001 && (t!=0))
                  t_inv=t
                  break;
              end             
          endfor
    end
    if (t_inv)
      
        #aux e o t onde a inversao ocorre
        try
            if !mod(i,2)
                aux=fzero(y, [t_inv-t_space(2), t_inv]);
                point_res=x(aux);
            else
                aux=fzero(x, [t_inv-t_space(2), t_inv]);
                point_res=y(aux);
            end           
        catch
            point_res=0;
        end_try_catch
    endif
endfunction

function cycle_res=cycle2_solver(y_1, t_space, data, A_1, A_2)
    x_1=quarter_poinc2_solver(4, y_1, t_space, data, A_2);
    y_2=quarter_poinc2_solver(3, x_1, t_space, data, A_1);
    x_2=quarter_poinc2_solver(2, y_2, t_space, data, A_2);
    cycle_res=quarter_poinc2_solver(1, x_2, t_space, data, A_1);
endfunction

function cycle_res=cycle2(y_1, t_space, data, A_1, A_2)
    x_1=quarter_poinc2(4, y_1, t_space, data, A_2);
    y_2=quarter_poinc2(3, x_1, t_space, data, A_1);
    x_2=quarter_poinc2(2, y_2, t_space, data, A_2);
    y_1=quarter_poinc2(1, x_2, t_space, data, A_1);
endfunction

hf=figure;
xlabel("x");
ylabel("y");
hold on;
plot([0, 0], [-10, 10], "color", "black");
plot([-10, 10], [0, 0], "color", "black");
plot([-10-data{2, 3}(1), 10-data{2, 3}(1)], [-10*M(2,1)-data{2, 3}(2), 10*M(2,1)-data{2, 3}(2)], "color", "blue");
plot([-10-data{2, 3}(1), 10-data{2, 3}(1)], [-10*M(2,2)-data{2, 3}(2), 10*M(2,2)-data{2, 3}(2)], "color", "blue");
c=0.2;
type=1;

tangents=find_elementary(A_1, A_2, data);
singularities=find_singularity(A_1, A_2, data);

scatter(0, tangents(1), "filled", "markerfacecolor", "k");
text(0-2*c, tangents(1)-c, strcat('�_1^y'));
scatter(tangents(2), 0, "filled", "markerfacecolor", "k");
text(tangents(2)-2*c, 0-c, strcat('�_1^x'));

for i=1:length(singularities)
    scatter(singularities(i,1), singularities(i,2), "filled", "markerfacecolor", "k");
    text(singularities(i,1)-2*c, singularities(i,2)-c, strcat('A_', num2str(i)));
endfor

#cycle2(8, n_space, data, A_1, A_2)
#cycle2(0.45, n_space, data, A_1, A_2)
#cycle2(0.01, n_space, data, A_1, A_2)
#cycle2(1.8, n_space, data, A_1, A_2)
#cycle2(1.3, n_space, data, A_1, A_2)
#cycle2(1.5, n_space, data, A_1, A_2)

#solver=@(y_1) cycle2_solver(y_1, n_space, data, A_1, A_2)-y_1;
#cycle2(fzero(solver, [1.3, 1.8]), n_space, data, A_1, A_2);

axis([-3,3,-2.5,2.5])

print(hf, "prov.png");
hold off;