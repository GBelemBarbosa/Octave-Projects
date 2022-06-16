n=13;
a=sqrt(1/2);
DIAG=[a,1,1,0,0,1,1,1,a,1,1,a,1];
Pb=[0,0,10,0,0,0,0,0,15,0,20,0,0];
ENVsup=[-1,0,0,-a,0,0,a,-1,0,0,0,-1,0,0,-a,0,a,-1,0,0,0,-a,0,0,0,-1,0,0];
ENVcolsup=[1,1,1,1,4,8,12,12,15,18,22,22,26,29];
ENVlinsup=[1,2,3,1,2,3,4,2,3,4,5,5,6,7,6,7,8,6,7,8,9,8,9,10,11,10,11,12];

ENVinf=[a,0,1,1,a,a,0,1,0,a,0,1,-a];
ENVlininf=[1,1,1,1,4,5,6,6,6,10,10,10,13,14];
ENVcolinf=[1,2,3,4,5,5,6,7,8,9,10,11,12];

for j=1:n-1
    for p=1:length(ENVcolinf)
        if ENVcolinf(p)==j && ENVinf(p)!=0
            c=-ENVinf(p)/DIAG(j);
            ENVinf(p)=-c;
            linha=2;
            while ENVlininf(linha+1)<=p
                linha++;
            end
            for x=1:length(ENVlinsup)
                if ENVlinsup(x)==j && ENVsup(x)!=0;
                    coluna=2;
                    while ENVcolsup(coluna+1)<=x
                        coluna++;
                    end
                    if linha>coluna
                        ENVinf(p+coluna-j)+=c*ENVsup(x);
                    elseif linha<coluna
                        ENVsup(x+linha-j)+=c*ENVsup(x);
                    else
                        DIAG(coluna)+=c*ENVsup(x);
                    end
                end
            end
        end
    end                                    
end
y=[];
x=[];
for j=1:n
    yj=Pb(j);
    y=[y;yj];
    for p=1:length(ENVcolinf)
        if ENVcolinf(p)==j && ENVinf(p)!=0
            linha=2;
            while ENVlininf(linha+1)<=p
                linha++;
            end
            Pb(linha)-=ENVinf(p)*yj;
        end
    end
end
for j=n:-1:1
    xj=y(j)/DIAG(j);
    x=[xj;x];
    for p=ENVcolsup(j):ENVcolsup(j+1)-1
        y(ENVlinsup(p))-=ENVsup(p)*xj;
    end
end
