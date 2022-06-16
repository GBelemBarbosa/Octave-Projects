A_plus=[3/4, -1; 1, 3/4];
cut=1;
gamma_plus=-A_plus(1, 1)/A_plus(1, 2);
c=0.72;
beta_minus=0.6;
p=0.09;
#gamma_plus=-cut*A_plus(1, 1)/A_plus(1, 2)

#Se for no sentido horario, colocar t_space decrescente (<=0)
#Mudar o modulo do max t_space e o numero de intervalos se necessario 
#(este ultimo no caso de um sistema passar pelo corte multiplas vezes entre dois passos)
data=cell(2, 2);
t_space=linspace(0, 12, 1200);
n_space=-t_space;

[data{2, 1}, data{2, 2}]=eig(A_plus);

function aux=half_poinc(i, Y, t_space, mode, cut, data)
    alpha_minus=real(data{i, 2}(1,1));
    beta_minus=imag(data{i, 2}(1,1));
    if (i==2)
            x=@(t) exp(alpha_minus*t)*(cos(beta_minus*t)*cut-Y*sin(beta_minus*t))-cut;
            y=@(t) exp(alpha_minus*t)*(sin(beta_minus*t)*cut+Y*cos(beta_minus*t));
    else
            m=real(data{i, 1}(2,1));           
            gamma_minus=alpha_minus/beta_minus;
            A_m_diag_dif=-2*m/gamma_minus;
            A_minus_12=1/gamma_minus;
            A_minus_21=-beta_minus*(alpha_minus+m^2/alpha_minus);
            x=@(t) exp(alpha_minus*t)*(cut*(cos(beta_minus*t)+A_m_diag_dif*sin(beta_minus*t)/(2*beta_minus))+Y*A_minus_12*sin(beta_minus*t)/beta_minus)-cut;
            y=@(t) exp(alpha_minus*t)*(Y*(cos(beta_minus*t)-A_m_diag_dif*sin(beta_minus*t)/(2*beta_minus))+cut*A_minus_21*sin(beta_minus*t)/beta_minus);
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
    
    #Retorna Y ao inves de t da inversao
    if (mode)
        aux=Y;
    endif
endfunction

function phi=calc_phi(tau, gamma)
    phi=1-exp(gamma*tau)*(cos(tau)-gamma*sin(tau));
endfunction

function g=calc_g(gamma_minus, alpha_minus, m, beta_minus, cut, data, n_space, t_space, c)
    M=[1, 1; m+alpha_minus*i, m-alpha_minus*i];
    data{1, 1}=M;
    g=(c+exp(gamma_minus*pi)*half_poinc(1, m-alpha_minus*gamma_minus, n_space, 1, cut, data))/(1+exp(gamma_minus*pi))-alpha_minus*gamma_minus-m;
endfunction
 
function f=calc_f(alpha_minus, gamma_plus, data, n_space, t_space, beta_minus, p, cut, c)
    J=[alpha_minus+beta_minus*i, 0; 0, alpha_minus-beta_minus*i];
    gamma_minus=alpha_minus/beta_minus
    data{1, 2}=J;
    m=fzero(@(m) calc_g(gamma_minus, alpha_minus, m, beta_minus, cut, data, n_space, t_space, c), [-100, 100])
    data{1, 1}=[1, 1; m+alpha_minus*i, m-alpha_minus*i];
    scatter(alpha_minus, m, "b");
    
    A_minus=[alpha_minus-m/gamma_minus, 1/gamma_minus; -beta_minus*(alpha_minus+m^2/alpha_minus), alpha_minus+m/gamma_minus];
    
    yc_minus=-A_minus(1, 1)/A_minus(1, 2)
    scatter(alpha_minus, yc_minus, "g");
    
    #yc_minus=-cut*A_minus(1, 1)/A_minus(1, 2);
    tau_plus=half_poinc(2, yc_minus, t_space, 0, cut, data);
    #tau_minus=half_poinc(2, yc_minus, t_space, 0, cut, data)*beta_plus;
    tau_minus=-half_poinc(1, yc_minus, n_space, 0, cut, data)*beta_minus;
    phi_plus_pos=calc_phi(tau_plus, gamma_plus);
    phi_plus_neg=calc_phi(tau_plus, -gamma_plus);
    phi_minus=calc_phi(tau_minus, gamma_minus);
    
    f=log(sin(tau_minus)*A_minus(1,2)*(exp(-gamma_plus*tau_plus)*phi_plus_pos+exp(gamma_plus*tau_plus)*phi_plus_neg)/(sin(tau_plus)*beta_minus*phi_minus))/tau_minus+gamma_minus+p
    scatter(alpha_minus, f, "r");
endfunction

hf=figure;
xlabel("\\alpha^-");
hold on;

alpha_minus=-0.018;
k=0;

#Busca de inversão de f
if (sign(calc_f(alpha_minus, gamma_plus, data, n_space, t_space, beta_minus, p, cut, c))>0)
    do  
        k++
        
        #Passo (alterar se necessario)
        alpha_minus*=1.1
        aux=calc_f(alpha_minus, gamma_plus, data, n_space, t_space, beta_minus, p, cut, c);
        
        #k>n condição de parada secundaria (alterar se necessario)
    until (sign(aux)<0 || k>10)
    start=alpha_minus;
    finish=alpha_minus/1.1;
else
    do
        k++
        
        #Passo (alterar se necessario)
        alpha_minus/=2
        aux=calc_f(alpha_minus, gamma_plus, data, n_space, t_space, beta_minus, p, cut, c);
        
        #k>n condição de parada secundaria (alterar se necessario)
    until (sign(aux)>0 || k>10)
    finish=alpha_minus;
    start=alpha_minus*2;
endif
alpha_minus=fzero(@(alpha_minus) calc_f(alpha_minus, gamma_plus, data, n_space, t_space, beta_minus, p, cut, c), [start, finish])
legend({"m^-", "y_c^-", "f(\\alpha^-)"}, "location", "southeast");
J=[alpha_minus+beta_minus*i, 0; 0, alpha_minus-beta_minus*i];
gamma_minus=alpha_minus/beta_minus
data{1, 2}=J;
m=fzero(@(m) calc_g(gamma_minus, alpha_minus, m, beta_minus, cut, data, n_space, t_space, c), [0, 100])
A_minus=[alpha_minus-m/gamma_minus, 1/gamma_minus; -beta_minus*(alpha_minus+m^2/alpha_minus), alpha_minus+m/gamma_minus]
text(alpha_minus+0.002, 0.005, "f(\\alpha^-)=0");
hold off;  
