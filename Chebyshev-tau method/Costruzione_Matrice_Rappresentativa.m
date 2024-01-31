function[F]=Costruzione_Matrice_Rappresentativa(N,vettcoeff,order)
%%%% N è il numero di polinomi di Chebyshev nella base per approx theta;
%%%% vettcoeff è il vettore dei coefficienti nello 
%%%% sviluppo della funzione f davanti ai polinomi di Chebyshev;
%%%% order è il grado del polinomio di Taylor di f;

theta=zeros(N,order+N);
primo=zeros(1,order+N);
secondo=zeros(1,order+N);
contributo=zeros(1,order+N);
for k=1:N
for i=1:order
    primo(1,i+k-1)=vettcoeff(i)/2;
    secondo(1,abs(i-k)+1)=vettcoeff(i)/2;
    contributo=primo+secondo+contributo;
    primo=zeros(1,order+N);
    secondo=zeros(1,order+N);
end
theta(k,:)=contributo;
contributo=zeros(1,order+N);
end
F=transpose(theta);
F=F(1:N,:);
end