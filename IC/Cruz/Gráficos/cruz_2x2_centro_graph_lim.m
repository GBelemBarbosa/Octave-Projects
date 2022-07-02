hf=figure;
xlabel("\alpha");
ylabel("y_{min}^{2*}/y_{min}^2");
hold on;

alfa=1;
beta=2;
a=-0.95;

y_2_min_inv=beta*(1+sqrt(1-a^2));
y_2_min_prev=2*alfa-(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2);

while alfa>-2*beta
    alfa_prev=alfa;
    alfa=alfa-0.1;
    y_2_min=2*alfa-(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2);
    
    line([alfa_prev, alfa], [y_2_min_prev, y_2_min], "color", "cyan", "linewidth", 1.5)
    line([alfa_prev, alfa], [y_2_min_inv, y_2_min_inv], "color", "magenta", "linewidth", 1.5)
    y_2_min_prev=y_2_min;
endwhile

f = @(alfa) 2*alfa-(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2)-y_2_min_inv;
alfa_null_prev=fzero(@(alfa) f(alfa), [-100,1]);  

while a<-0.1
    a_prev=a;
    a=a+0.05;
    f = @(alfa) 2*alfa-(2*a^2*beta - sqrt(4*a^4*beta^2 + 4*a^2*(4*alfa^2 - 8*alfa*beta + 3*beta^2)))/(2*a^2)-y_2_min_inv;
    alfa_null=fzero(@(alfa) f(alfa), [-100,1]);    
    
    line([a_prev, a], [alfa_null_prev, alfa_null], "color", "yellow", "linewidth", 1.5)
    alfa_null_prev=alfa_null;
endwhile
hold off;