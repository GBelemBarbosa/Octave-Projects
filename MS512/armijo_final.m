
   
clear all;

A_2=[3/4, -1; 1, 3/4];
global  cut=1 e gamma_plus=-A_2(1, 1)/A_2(1, 2) e C=0.72 e beta_minus=0.6 e folga=0.09;
#gamma_plus=-cut*A_2(1, 1)/A_2(1, 2)

#A_1=(M*J)/M

#Se for no sentido horario, colocar t_space decrescente (<=0)
#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar pelo corte multiplas vezes entre dois passos)
global matrices=cell(2, 4) e t_space=linspace(0, 12, 1200) e n_space=-t_space;

h=-0.018
J=[h+beta_minus*i, 0; 0, h-beta_minus*i];
matrices{1, 2}=J;
A_1=[h-beta_minus*C/h, beta_minus/h; -beta_minus*(h+C^2/h), h+beta_minus*C/h]
matrices{1, 4}=[0; 0];
matrices{1, 3}=A_1\matrices{1, 4};

[matrices{2, 1}, matrices{2, 2}]=eig(A_2);
matrices{2, 4}=[0; 0];
matrices{2, 3}=A_2\matrices{2, 4};

function aux=find_Y(i, Y, t_space, mode)
    global cut e matrices;
    x=@(t) real(((matrices{i, 1}*expm(matrices{i, 2}*t)/matrices{i, 1})*([cut; Y]+matrices{i, 3})-matrices{i, 3})(1))-cut;
    xy=@(t) real((matrices{i, 1}*expm(matrices{i, 2}*t)/matrices{i, 1})*([cut; Y]+matrices{i, 3})-matrices{i, 3});
    t_inv=0;
    aux=0;
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
    if (mode)
        aux=Y;
    endif
endfunction

function phi=calc_phi(tau, gamma)
    phi=1-exp(gamma*tau)*(cos(tau)-gamma*sin(tau));
endfunction

function resto=res(gamma_minus, x, beta_minus)
    global matrices e n_space e t_space e C;
    M=[1, 1; x+gamma_minus*beta_minus*i, x-gamma_minus*beta_minus*i];
    matrices{1, 1}=M;
    resto=(C+exp(gamma_minus*pi)*find_Y(1, -gamma_minus*beta_minus*gamma_minus+x, n_space, 1))/(1+exp(gamma_minus*pi))-gamma_minus*beta_minus*gamma_minus-x;
endfunction
 
function f12=calc_f12(g)
    global gamma_plus e matrices e n_space e t_space e beta_minus e gamma_minus e folga;
    J=[g+beta_minus*i, 0; 0, g-beta_minus*i];
    gamma_minus=g/beta_minus
    matrices{1, 2}=J;
    x=fzero(@(x) res(gamma_minus, x, beta_minus), [-100, 100])
    A_1=[g-beta_minus*x/g, beta_minus/g; -beta_minus*(g+x^2/g), g+beta_minus*x/g];
    
    yc_minus=-A_1(1, 1)/A_1(1, 2)
    #yc_minus=-cut*A_1(1, 1)/A_1(1, 2);
    tau_plus=find_Y(2, yc_minus, t_space, 0);
    #tau_minus=find_Y(2, yc_minus, t_space, 0)*beta_plus;
    tau_minus=-find_Y(1, yc_minus, n_space, 0)*beta_minus;
    phi_plus_pos=calc_phi(tau_plus, gamma_plus);
    phi_plus_neg=calc_phi(tau_plus, -gamma_plus);
    phi_minus=calc_phi(tau_minus, gamma_minus);
    f12=log(sin(tau_minus)*A_1(1,2)*(exp(-gamma_plus*tau_plus)*phi_plus_pos+exp(gamma_plus*tau_plus)*phi_plus_neg)/(sin(tau_plus)*beta_minus*phi_minus))/tau_minus+gamma_minus+folga
    scatter(g, f12);
    if (sign(f12)>=0)
        display("gamma_minus=<f12<0");
    endif
endfunction

hf=figure;
xlabel ("\gamma_-");
ylabel ("f_{12}-\gamma_-");
hold on;

k=0;
if (sign(calc_f12(h))>0)
    finish=h
    do  
        k++
        h*=1.07
        aux=calc_f12(h);
    until (sign(aux)<0 || k>10)
    start=h;
else
    start=h
    do
        k++
        h/=2
        aux=calc_f12(h);
    until (sign(aux)>0 || k>10)
    finish=h;
endif
g=fzero(@(g) calc_f12(g), [start, finish])
J=[g+beta_minus*i, 0; 0, g-beta_minus*i];
gamma_minus=g/beta_minus
matrices{1, 2}=J;
x=fzero(@(x) res(gamma_minus, x, beta_minus), [0, 100])
A_1=[g-beta_minus*x/g, beta_minus/g; -beta_minus*(g+x^2/g), g+beta_minus*x/g]
text(g+0.002, 0.005, "f_{12}=\gamma_-");
hold off;

new=figure;
y_0_plus=@(y_1) -(y_1+gamma_plus*(1+exp(gamma_plus*pi)))/exp(gamma_plus*pi);
yc_minus=-(A_1(1, 1)+matrices{1, 4}(1))/A_1(1, 2)
y_1_minus=@(y_0) -exp(gamma_minus*pi)*y_0+(A_1(2,2)/A_1(1,2))*(1+exp(gamma_minus*pi));
Y=find_Y(1, yc_minus, n_space, 1)
y_space=linspace(Y, Y+1.2, 50);
hold on;
axis([Y, Y+0.05, yc_minus-0.5, yc_minus+0.5])
line([y_space(1), y_space(end)], [y_0_plus(y_space(1)), y_0_plus(y_space(end))], "color", "g", "linestyle", "--");
line([y_space(1), y_space(end)], [y_1_minus(y_space(1)), y_1_minus(y_space(end))], "color", "b", "linestyle", "--");
last=[yc_minus, find_Y(2, Y, n_space, 1)];
(y_1_minus(y_space(end))-y_1_minus(y_space(1)))/1.2
scatter(Y, yc_minus)
scatter(Y, last(2))
h=y_space(2)-y_space(1); 
for y=y_space(2:end)
    Y=find_Y(1, y, t_space, 1);
    line([y, y-h], [Y, last(1)], "color", "b");
    last(1)=Y;
    Y=find_Y(2, y, n_space, 1);
    line([y, y-h], [Y, last(2)], "color", "g");
    last(2)=Y;
endfor        
