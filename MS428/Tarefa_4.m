restri=[-3,4,1,0,0;1,-1,0,1,0;1,1,0,0,1];
res=[12;4;6];
hf=figure;
xlabel ("x1");
ylabel ("x2");
solsol=[];
k=0;
hold on;
for i=1:size(restri,2)
    for j=i+1:size(restri,2)
        k++;
        aux=restri;
        aux(:,[i,j])=[];
        sol=aux\res;
        if (aux*sol-res<=0.001)
            solsol=[solsol,[sol; i; j]];
            if (i==1)
                if (j==2)
                    text(0.1, 0.1, strcat("X", num2str(k)));
                    scatter(0, 0, "filled");
                else
                    text(0.1, sol(1)+0.1, strcat("X", num2str(k)));
                    scatter(0, sol(1), "filled");
                end
            elseif (i==2)
                text(sol(1)+0.1, 0.1, strcat("X", num2str(k)));
                scatter(sol(1), 0, "filled");
            else
                text(sol(1)+0.1, sol(2)+0.1, strcat("X", num2str(k)));
                scatter(sol(1), sol(2), "filled");
            end
        end
    endfor
endfor
for i=1:size(restri,1)
    if (restri(i,2)!=0)
        fplot(@(x) (res(i)-restri(i,1)*x)/restri(i,2),[0,6.5]);
    end
endfor
line([0.01,0.01], [-4,6]);
line([0,6.5], [0,0]);
axis([-0.5, 6.5]);  
legend("hide");
print (hf, "c.png");
hold off;