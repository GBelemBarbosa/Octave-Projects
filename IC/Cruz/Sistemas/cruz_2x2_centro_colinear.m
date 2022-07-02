a=-10/11;

A_1=[0,  -4/(4*a); a,  0];
A_2=[0,  1; -1,  0];

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

beta=2;
#b_1
data{1, 4}=[beta*x*4/(4*a); -beta*y*a];
data{1, 3}=A_1\data{1, 4};
#b_1*A_1^(-1)

[data{2, 1}, data{2, 2}]=eig(A_2);

alfa=3;
#b_2
data{2, 4}=[-alfa*x; alfa*y];
#d=b_2*A_2^(-1)
data{2, 3}=A_2\data{2, 4};

function tangents=find_elementary(A_1, A_2, data)
    Atil_1_y=-data{1, 4}(1)/A_1(1, 2);
    Atil_1_x=-data{1, 4}(2)/A_1(2, 1);
    Atil_2_y=-data{2, 4}(1)/A_2(1, 2);
    Atil_2_x=-data{2, 4}(2)/A_2(2, 1);
    tangents=[Atil_1_y, Atil_1_x, Atil_2_y, Atil_2_x]
endfunction

function singularities=find_singularity(A_1, A_2, data)
    Asin_1=transpose(-A_1\data{1,4});
    Asin_2=transpose(-A_2\data{2,4});
    singularities=[Asin_1; Asin_2]
endfunction

function point_res=quarter_poinc1(i, point, t_space, data, A)
    d=data{2-mod(i,2),3};
    A_m_diag_dif=A(1,1)-A(2,2);
    if i>2
            point_p_b=point+d(2);
            x=@(t) d(1)*(cos(t)+A_m_diag_dif*sin(t)/2)+point_p_b*A(1,2)*sin(t)-d(1);
            y=@(t) point_p_b*(cos(t)-A_m_diag_dif*sin(t)/2)+d(1)*A(2,1)*sin(t)-d(2);
    else
            point_p_b=point+d(1);
            x=@(t) point_p_b*(cos(t)+A_m_diag_dif*sin(t)/2)+d(2)*A(1,2)*sin(t)-d(1);
            y=@(t) d(2)*(cos(t)-A_m_diag_dif*sin(t)/2)+point_p_b*A(2,1)*sin(t)-d(2);
    endif
    
    #Busca da inversao
    t_inv=0;
    if i>2
        last=[0, point]
    else
        last=[point,0]
    endif
    if mod(i,4)<2 && (point>0)
          for t=t_space
              Y=y(t);
              X=x(t);
              line([last(1), X], [last(2), Y], "color", "black", "linewidth", 1.5)
              last=[X, Y];
              if ((X*(i-2)>-0.00001) && (t!=0))
                  t_inv=t
                  break;
              end
          endfor
    elseif (point>0)
          for t=t_space
              Y=y(t);
              X=x(t);
              line([last(1), X], [last(2), Y], "color", "black", "linewidth", 1.5)
              last=[X, Y];
              if (Y*(i-2.5)<0 && (t!=0))
                  t_inv=t
                  break;
              end
          endfor
    end
    if (t_inv)
      
        #aux e o t onde a inversao ocorre
        try
            if mod(i,4)>1
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

function y_3=cycle1(y_1, t_space, data, A_1, A_2)
    y_2=quarter_poinc1(4, y_1, t_space, data, A_2);
    x_1=quarter_poinc1(3, y_2, t_space, data, A_1);
    x_2=quarter_poinc1(2, x_1, t_space, data, A_2);
    y_3=quarter_poinc1(1, x_2, t_space, data, A_1)
endfunction

function point_res=quarter_poinc3(i, point, t_space, data, A)
    d=data{i,3};
    A_m_diag_dif=A(1,1)-A(2,2);
    point_p_b=point+d(2);
    x=@(t) d(1)*(cos(t)+A_m_diag_dif*sin(t)/2)+point_p_b*A(1,2)*sin(t)-d(1);
    y=@(t) point_p_b*(cos(t)-A_m_diag_dif*sin(t)/2)+d(1)*A(2,1)*sin(t)-d(2);
    
    #Busca da inversao
    t_inv=0;
    last=[0, point]
    for t=t_space
        Y=y(t);
        X=x(t);
        line([last(1), X], [last(2), Y], "color", "black", "linewidth", 1.5)
        last=[X, Y];
        if ((X*(1.5-i)<0) && (t!=0))
            t_inv=t
            break;
        end
    endfor
    if (t_inv)
      
        #aux e o t onde a inversao ocorre
        try
            aux=fzero(x, [t_inv-t_space(2), t_inv]);
            point_res=y(aux);
        catch
            point_res=0;
        end_try_catch
    endif
endfunction

function cycle_res=cycle3(y_1, t_space, data, A_1, A_2)
    y_2=quarter_poinc3(2, y_1, t_space, data, A_2);
    y_3=quarter_poinc3(1, y_2, t_space, data, A_1);
endfunction

function point_res=quarter_poinc2(i, point, t_space, data, A)
    d=data{2-mod(i,2),3};
    A_m_diag_dif=A(1,1)-A(2,2);
    if !mod(i,2)
            point_p_b=point+d(2);
            x=@(t) d(1)*(cos(t)+A_m_diag_dif*sin(t)/2)+point_p_b*A(1,2)*sin(t)-d(1);
            y=@(t) point_p_b*(cos(t)-A_m_diag_dif*sin(t)/2)+d(1)*A(2,1)*sin(t)-d(2);
    else
            point_p_b=point+d(1);
            x=@(t) exp(alpha*t)*(point_p_b*(cos(t)+A_m_diag_dif*sin(t)/2)+d(2)*A(1,2)*sin(t))-d(1);
            y=@(t) exp(alpha*t)*(d(2)*(cos(t)-A_m_diag_dif*sin(t)/2)+point_p_b*A(2,1)*sin(t))-d(2);
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
c=0.2;
type=1;

tangents=find_elementary(A_1, A_2, data);
singularities=find_singularity(A_1, A_2, data);

if type
    for i=1:length(tangents)
        if mod(i, 2)
            scatter(0, tangents(i), "filled", "markerfacecolor", "k");
            text(0+c, tangents(i)-c, strcat('Ã_', num2str(round(i/2)), '^y'));
        else
            scatter(tangents(i), 0, "filled", "markerfacecolor", "k");
            text(tangents(i)+c, 0-c, strcat('Ã_', num2str(round(i/2)), '^x'));
        endif
    endfor

    for i=1:length(singularities)
        scatter(singularities(i,1), singularities(i,2), "filled", "markerfacecolor", "k");
        text(singularities(i,1)+c, singularities(i,2)-c, strcat('A_', num2str(i)));
    endfor
else
    scatter(0, tangents(1), "filled", "markerfacecolor", "k");
    text(0+c, tangents(1)-c, "Ã_1^y=Ã_2^y");
    scatter(tangents(3), 0, "filled", "markerfacecolor", "k");
    text(tangents(3)+c, 0-c, "Ã_1^x=Ã_2^x");

    scatter(singularities(1,1), singularities(1,2), "filled", "markerfacecolor", "k");
    text(singularities(1,1)+c, singularities(1,2)-c, "A_1=A_2");
end

#y=beta*(1+sqrt(1-a^2))
#text(-2*c, y, "y_{min}^{1*}");
#cycle1(beta+sqrt(beta^2*(1-a^2)), n_space, data, A_1, A_2)

#y=2*alfa-beta + sqrt(beta^2 + a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2));
#cycle1(y, n_space, data, A_1, A_2)
#text(c, y, "y_{min}^1");

#y=2*alfa+0.001;
#y=2*alfa-(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2);
#cycle2(y, n_space, data, A_1, A_2)
#text(c, y, "y_{lim}^{1-2}");

#cycle2(2*alfa+1, n_space, data, A_1, A_2)
    
#cycle2(2*alfa-(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2), n_space, data, A_1, A_2)

cycle3(2*alfa-beta+sqrt(-beta^2*(a^2-1)), n_space, data, A_1, A_2)
cycle3(4.01, n_space, data, A_1, A_2)



axis([-1,5.5,-1,5], "equal")

print(hf, "prov.png");
hold off;