a=2;

function res=cycle_2(y_1, alfa, beta, a)
  
  x_1 = 2*alfa-y_1;
  
  y_2 = beta - sqrt(a^2*x_1*(x_1 - 2*beta) + beta^2);
  
  x_2 = 2*alfa-y_2;
  
  res = x_2 - beta - sqrt(a^2*(a^2*beta^2 + y_1*(y_1 - 2*beta)))/a^2;
    
endfunction

function res=cycle_1(y_1, alfa, beta, a)
  
  y_2 = 2*alfa-y_1;
  
  x_1 =(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*y_2*(y_2 - 2*beta)))/(2*a^2);
  
  x_2 = 2*alfa-x_1;
    
  res = x_2 - beta - sqrt(a^2*(a^2*beta^2 + y_1*(y_1 - 2*beta)))/a^2;
    
endfunction

function res=cycle_3(y_1, alfa, beta, a)
  
  y_2 = 2*alfa-y_1;
  
  res = y_2 - 2*beta+y_1;
    
endfunction

hf=figure;
xlabel("y_1");
ylabel("r_i");
hold on;
line([0, 0], [0, 0], "color", "green", "linewidth", 1.5)
line([0, 0], [0, 0], "color", "red", "linewidth", 1.5)
line([0, 0], [0, 0], "color", "blue", "linewidth", 1.5)
legend({"Tipo 3", "Tipo 1", "Tipo 2"}, "location", "east");
legend("autoupdate", "off");

last=[0,0];
h=1/100;

alfa=0.7;
beta=2;
a=-10/11;

y_3_min=0;
y_1_min=0;
y_2_min=0;

if alfa>beta/2
    if alfa>=beta
        #Tipo 3
        #y_1 da menor órbita
        y_1=2*alfa-beta
        cycle_3(y_1, alfa, beta, a)
        last=[y_1, cycle_3(y_1, alfa, beta, a)];
        
        y_3_min=y_1;

        #y_1 da órbita limítrofe 3-1
        y_max=2*alfa-beta+sqrt(-beta^2*(a^2-1))
        
        y_3_1_lim=y_max;

        while y_1<y_max
            y_1+=h;
            y_f=cycle_3(y_1, alfa, beta, a);
            
            line([last(1), y_1], [last(2), y_f], "color", "green", "linewidth", 1.5)
            last=[y_1, y_f];
        endwhile        
    else
        #y_1 da menor órbita
        y_1_min_inv=beta*(1+sqrt(1-a^2))
        y_1_min=2*alfa-beta + sqrt(beta^2 + a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2))
        
        y_max=max(y_1_min_inv, y_1_min)          
    endif

    #Tipo 1
    y_1=y_max
    cycle_3(y_1, alfa, beta, a)
    cycle_1(y_1, alfa, beta, a)
    last=[y_1, cycle_1(y_1, alfa, beta, a)];

    #y_1 da órbita limítrofe 1-2
    y_max=2*alfa
    y_1_2_lim=y_max;

    while y_1<y_max
        y_1+=h;
        y_f=cycle_1(y_1, alfa, beta, a);
        
        line([last(1), y_1], [last(2), y_f], "color", "red", "linewidth", 1.5)
        last=[y_1, y_f];
    endwhile
else
    #y_1 da menor órbita
    y_2_min=2*alfa-(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2)
    y_2_min_inv=y_2_min

    y_max=max(y_2_min, y_2_min_inv)    
endif

#Tipo 2
y_1=max(y_max, y_1_min_inv)
last=[y_1, cycle_2(y_1, alfa, beta, a)];

#y_1 máximo
y_max=y_1+2

while y_1<y_max
    y_1+=h;
    y_f=cycle_2(y_1, alfa, beta, a);
    line([last(1), y_1], [last(2), y_f], "color", "blue", "linewidth", 1.5)
    last=[y_1, y_f];
endwhile

c=0.05;

if y_3_min || y_1_min
    y_max_graph=max([cycle_1(y_1_2_lim, alfa, beta, a), cycle_2(y_1_2_lim, alfa, beta, a), cycle_2(y_max, alfa, beta, a)]);
    y_min=min([cycle_1(y_1_2_lim, alfa, beta, a), cycle_2(y_1_2_lim, alfa, beta, a), cycle_2(y_max, alfa, beta, a)]);

    if y_3_min
        y_max_graph=max([y_max_graph, cycle_3(y_3_min, alfa, beta, a), cycle_1(y_3_1_lim, alfa, beta, a), cycle_3(y_3_1_lim, alfa, beta, a)]);
        y_min=min([y_max_graph, cycle_3(y_3_min, alfa, beta, a), cycle_1(y_3_1_lim, alfa, beta, a), cycle_3(y_3_1_lim, alfa, beta, a)]);
        
        line([y_3_min, y_3_min], [y_min-1, y_max_graph+1], "linestyle", "--", "color", "black", "linewidth", 1.5)
        text(y_3_min+c, y_max_graph+c, "y_{min}^3");
        
        line([y_3_1_lim, y_3_1_lim], [y_min-1, y_max_graph+1], "linestyle", "--", "color", "black", "linewidth", 1.5)
        text(y_3_1_lim+c, y_max_graph+c, "y_{lim}^{3-1}");
    endif
        
    if y_1_min
        y_max_graph=max(y_max_graph, real(cycle_1(y_1_min, alfa, beta, a)));
        if y_1_min_inv<y_1_2_lim
            y_max_graph=max(y_max_graph, cycle_1(y_1_min_inv, alfa, beta, a));
        endif
        y_min=min([y_max_graph, cycle_1(y_1_min, alfa, beta, a), cycle_2(y_max, alfa, beta, a)]);
        
        line([y_1_min, y_1_min], [y_min-1, y_max_graph+1], "linestyle", "--", "color", "black", "linewidth", 1.5)
        text(y_1_min-6*c, y_max_graph+c, "y_{min}^1");
        
        line([y_1_min_inv, y_1_min_inv], [y_min-1, y_max_graph+1], "linestyle", "--", "color", "black", "linewidth", 1.5)
        text(y_1_min_inv-5*c, y_max_graph+c, "y_{min}^{1*}");
    endif
    
    line([y_1_2_lim, y_1_2_lim], [y_min-1, y_max_graph+1], "linestyle", "--", "color", "black", "linewidth", 1.5)
    text(y_1_2_lim+c, y_max_graph+c, "y_{lim}^{1-2}");
endif

if y_2_min
    y_max_graph=max(cycle_2(y_2_min, alfa, beta, a), cycle_2(y_max, alfa, beta, a));
    y_min=min(cycle_2(y_2_min, alfa, beta, a), cycle_2(y_max, alfa, beta, a));
    
    line([y_2_min, y_2_min], [y_min-1, y_max_graph+1], "linestyle", "--", "color", "black", "linewidth", 1.5)
    text(y_2_min+c, y_max_graph+c, "y_{lim}^{1-2}");
endif
  
x_min=max([y_1_min, y_2_min, y_3_min])-1;
  
plot([0, 0], [y_min-1, y_max_graph+1], "color", "black");
plot([x_min, y_max], [0, 0], "color", "black");

axis([x_min, y_max, y_min-1, y_max_graph+1], "equal")

print(hf, "prov.png");
hold off;

#y_1=2*alfa -(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2)

#cycle_2(y_1, alfa, beta, a)







