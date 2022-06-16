clear all;

A_2=[3/4, -1; 1, 3/4];
global  cut=1 e gamma_plus=-A_2(1, 1)/A_2(1, 2) e C=0.72 e beta_minus=0.6 e folga=0;
#gamma_plus=-cut*A_2(1, 1)/A_2(1, 2)

#A_1=(M*J)/M

#Se for no sentido horario, colocar t_space decrescente (<=0)
#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar pelo corte multiplas vezes entre dois passos)
global matrices=cell(2, 4) e t_space=linspace(0, 12, 1200) e n_space=-t_space;

#A_1=[h-beta_minus*C/h, beta_minus/h; -beta_minus*(h+C^2/h), h+beta_minus*C/h]
matrices{1, 4}=[0; 0];
matrices{1, 3}=[0; 0];
#matrices{1, 3}=A_1\matrices{1, 4};

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

function resto=res(gamma_minus, alpha_minus, m, beta_minus)
    global matrices e n_space e t_space e C;
    M=[1, 1; m+alpha_minus*i, m-alpha_minus*i];
    matrices{1, 1}=M;
    resto=(C+exp(gamma_minus*pi)*find_Y(1, m-alpha_minus*gamma_minus, n_space, 1))/(1+exp(gamma_minus*pi))-alpha_minus*gamma_minus-m;
endfunction
 
function f12=calc_f12(alpha_minus)
    global last e gamma_plus e matrices e n_space e t_space e beta_minus e gamma_minus e folga;
    J=[alpha_minus+beta_minus*i, 0; 0, alpha_minus-beta_minus*i];
    gamma_minus=alpha_minus/beta_minus
    matrices{1, 2}=J;
    m=fzero(@(m) res(gamma_minus, alpha_minus, m, beta_minus), [-100, 100])
    plot([alpha_minus, alpha_minus+0.02*(alpha_minus!=-0.0005)], [m, last(1)], "b");
    
    A_1=[alpha_minus-m/gamma_minus, 1/gamma_minus; -beta_minus*(alpha_minus+m^2/alpha_minus), alpha_minus+m/gamma_minus];
    #matrices{1, 3}=A_1\matrices{1, 4};

    yc_minus=-A_1(1, 1)/A_1(1, 2)
    plot([alpha_minus, alpha_minus+0.02*(alpha_minus!=-0.0005)], [yc_minus, last(2)], "g");
    
    #yc_minus=-cut*A_1(1, 1)/A_1(1, 2);
    tau_plus=find_Y(2, yc_minus, t_space, 0)
    #tau_minus=find_Y(2, yc_minus, t_space, 0)*beta_plus;
    tau_minus=-find_Y(1, yc_minus, n_space, 0)*beta_minus;
    phi_plus_pos=calc_phi(tau_plus, gamma_plus);
    phi_plus_neg=calc_phi(tau_plus, -gamma_plus);
    phi_minus=calc_phi(tau_minus, gamma_minus);
    
    f12=log(sin(tau_minus)*A_1(1,2)*(exp(-gamma_plus*tau_plus)*phi_plus_pos+exp(gamma_plus*tau_plus)*phi_plus_neg)/(sin(tau_plus)*beta_minus*phi_minus))/tau_minus+gamma_minus+folga
    plot([alpha_minus, alpha_minus+0.02*(alpha_minus!=-0.0005)], [f12, last(3)], "r");
    last=[m, yc_minus, f12];
    if (sign(f12)>=0)
        display("gamma_minus=<f12<0");
    endif
endfunction

hf=figure;
xlabel("\\alpha^-");
hold on;

global last=[0.300038883047789, 0.300038633047789, 1.645341196444514];
h=-0.0005;
k=0;
if (sign(calc_f12(h))>0)
    do  
        k++
        h-=0.02
        aux=calc_f12(h);
    until (k>10)
    start=h;
    finish=h/1.8;
else
    do
        k++
        h/=2
        aux=calc_f12(h);
    until (sign(aux)>0 || k>10)
    finish=h;
    start=h*=2;
endif
#alpha_minus=fzero(@(alpha_minus) calc_f12(alpha_minus), [start, finish])
#plot([-0.484815333006630, -0.484815333006630], [-1.5, 2], "k");
#text(-0.47, 1.2, "\\alpha^-=\\alpha^-_*")
plot([0, 0], [-1.5, 2], "k");
legend({"m^-", "y_c^-", "f(\\alpha^-)"}, "location", "southeast");
ylim("auto")
xlim("auto")
#print(hf, "curve.png");
hold off;
J=[alpha_minus+beta_minus*i, 0; 0, alpha_minus-beta_minus*i];
gamma_minus=alpha_minus/beta_minus
matrices{1, 2}=J;
m=fzero(@(m) res(gamma_minus, m, beta_minus), [0, 100])
A_1=[alpha_minus-m/gamma_minus, 1/gamma_minus; -beta_minus*(alpha_minus+m^2/alpha_minus), alpha_minus+m/gamma_minus]
text(alpha_minus+0.002, 0.005, "f_{12}=\\gamma_-");