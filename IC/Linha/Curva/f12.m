clear all;

A_2=[3/4, -1; 1, 3/4];
global  cut=1 e gamma_plus=-A_2(1, 1)/A_2(1, 2) e M=[1, 1; gamma_plus-0.1-i, gamma_plus-0.1+i];
#gamma_plus=-cut*A_2(1, 1)/A_2(1, 2)
J=[i, 0; 0, -i]
#A_1=(M*J)/M
A_1=[gamma_plus-0.1, -1; 1+(gamma_plus-0.1)^2, -gamma_plus+0.1]

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

function aux=find_Y(i, Y, t_space)
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
endfunction

function phi=calc_phi(tau, gamma)
      phi=1-exp(gamma*tau)*(cos(tau)-gamma*sin(tau));
endfunction

beta_minus=imag(J(1,1))
beta_plus=A_2(2,1)
gamma_minus=real(J(1,1))
#gamma_minus=real(J(1,1))/beta_minus
yc_minus=-A_1(1, 1)/A_1(1, 2)
#yc_minus=-cut*A_1(1, 1)/A_1(1, 2)
tau_plus=find_Y(2, yc_minus, t_space)
#tau_plus=find_Y(2, yc_minus, t_space)/beta_plus
phi_plus=calc_phi(tau_plus, gamma_plus)
#tau_minus=-find_Y(1, yc_minus, n_space)
#tau_minus=-find_Y(1, yc_minus, n_space)/beta_minus

function f12=calc_f12(g)
    global gamma_plus e matrices e n_space e t_space e M;
    J=[g+i, 0; 0, g-i];
    A_1=[g+gamma_plus-0.1, -1; 1+(gamma_plus-0.1)^2, g-gamma_plus+0.1];
    #matrices{1, 1}=M;
    matrices{1, 2}=J;
    yc_minus=A_1(1, 1);
    #yc_minus=-cut*A_1(1, 1)/A_1(1, 2);
    tau_plus=find_Y(2, yc_minus, t_space);
    #tau_minus=find_Y(2, yc_minus, t_space, 0)*beta_plus;
    tau_minus=-find_Y(1, yc_minus, n_space);
    #tau_minus=-find_Y(1, yc_minus, n_space)*beta_minus;
    phi_plus_pos=calc_phi(tau_plus, gamma_plus);
    phi_plus_neg=calc_phi(tau_plus, -gamma_plus);
    phi_minus=calc_phi(tau_minus, g);
    f12=log(sin(tau_minus)*A_1(1,2)*(exp(-gamma_plus*tau_plus)*phi_plus_pos+exp(gamma_plus*tau_plus)*phi_plus_neg)/(sin(tau_plus)*phi_minus))/tau_minus+g
    #f12=-log(sin(tau_minus)*A_1(1,2)*(exp(-gamma_plus*tau_plus)*phi_plus_pos+exp(gamma_plus*tau_plus)*phi_plus_neg)/(sin(tau_plus)*beta_minus*phi_minus))/tau_minus-gamma_minus
    scatter(g, f12);
    if (sign(f12)>=0)
        display("gamma_minus=<f12<0");
    endif
endfunction

hf=figure;
xlabel ("\gamma_-");
ylabel ("f_{12}-\gamma_-");
hold on;

h=-0.1;
if (sign(calc_f12(h))>0)
    finish=h;
    do
        h*=2;
        aux=calc_f12(h);
    until (sign(aux)<0)
    start=h;
else
    start=h;
    do
        h/=2;
        aux=calc_f12(h);
    until (sign(aux)>0)
    finish=h;
endif    
g=fzero(@(g) calc_f12(g), [start, finish]);
text(g+0.002, 0.005, "f_{12}=\gamma_-");
hold off;    
A_1=[g+gamma_plus-0.1, -1; 1+(gamma_plus-0.1)^2, g-gamma_plus+0.1]