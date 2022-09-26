cIni=200000;
mFin=307600;

mc=mFin/cIni

jc=@(j, n) (1+j)^n

sol=@(j1, j2, n, mFin) log(mFin/jc(j2, n))/log(jc(j1, 1)/jc(j2, 1))

sol(4/100, 5/100, 10, mc)

G=@(m, j1, j2, x) (x<m)*jc(j1, x)+(x>=m)*jc(j1, m)*jc(j2, x-m)

G(6, 4/100, 5/100, 10)*cIni

T=[]
for t=1:20
  T=[T; t, G(6, 4/100, 5/100, t)];
end

plot(T(:,1), T(:,2))
hold on;

T=[]
for t=1:20
  T=[T; t, G(6, 5/100, 4/100, t)];
end

plot(T(:,1), T(:,2))
hold off;

real=@(m, inf) m/inf

T=[-4; -3; -2; -1; 0]
Inf=[1.02; 1.03; 1.029; 1.032; 1.025]
L=length(T)

a=(L*sum(T.*Inf)-sum(T)*sum(Inf))/(L*sum(T.^2)-sum(T)^2)

b=mean(Inf)-a*mean(T)

fit=@(x) a*x+b

fig=figure
line([T(1), T(end)], [fit(T(1)), fit(T(end))])
hold on;
scatter(T, Inf)
T=[ones(size(T), 1), T]
[p, e_var, r, p_var, fit_var] = LinearRegression (T, Inf)