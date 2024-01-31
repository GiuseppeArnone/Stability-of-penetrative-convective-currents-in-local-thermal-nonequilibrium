function[F,err_approx_Taylor,err_approx_Cheby]=Matrice_Rappresentativa_Esponenziale(fsym,c,N)
%%%% fsym è la funzione esponenziale che vogliamo sviluppare in serie di
%%%% Chebyshev e di cui vogliamo la matrice rappresentativa;
%%%% c è il coefficiente che compare nell'esponenziale;
%%%% N è la dimensione della base di polinomi di Chebyshev;
c=round(c);
order=c+2;
if c==0
    order=1;
end
if order >=8
    order=8;
end
f=fsym;
t=taylor(f,'ExpansionPoint', 0,'Order',order);
t=expand(t);
ff=matlabFunction(f);
tt=matlabFunction(t);
zz=0:0.01:1;
if c~=0
    err_approx_Taylor=max(abs(ff(zz)-tt(zz)));
else
    err_approx_Taylor=0;
end

coeff=CoeffExtraction(t,order); % Estrazione coefficienti del polinomio

%% Potenze di z in funzione dei polinomi di Chebyshev

cheb=zeros(order,order); % è la matrice per scrivere
% le potenze come combinazione di polinomi di Chebyshev.
% la seconda riga è x^1, la terza è x^2, la quarta è x^3

T=zeros(1,order); % vettore per i coefficienti davanti a T_0, T_1, T_2 ecc

for i = 2:order
    for j=1:i       % formula ricorsiva per z^(i-1)
        if (i-j)/2==fix((i-j)/2)
            T(1,j)=nchoosek(i-1,(i-j)/2);
            if j==1
                T(1,:)=T(1,:)/2;
            end
        end
    end
    T(1,:)= T(1,:)*2^(2-i);
    cheb(i,:) = T(1,:);
    T=zeros(1,order);
end
cheb(1,1)=1;

%% Costruzione vettore dei coefficienti davanti ai polinomi di Chebyshev nello sviluppo di f

vettcoeff=zeros(1,order); % Vettore dei coefficienti
% nello sviluppo di f davanti ai polinomi di Chebyshev
for i = 1:order
    vettcoeff=cheb(i,:)*coeff(i)+vettcoeff;
end

%% Stima della bontà dello sviluppo di f in polinomi di Chebyshev
syms x
chebpoly=chebyshevT(0:order-1,x);
g=0;
for i = 1:order
    g=chebpoly(i)*vettcoeff(i)+g;
end
gg=matlabFunction(g); % Sviluppo di f nella base dei polinomi di Chebyshev
ff=matlabFunction(f);
zz=0:0.01:1;
if c~=0
    err_approx_Cheby=max(abs(ff(zz)-gg(zz)));% errore commesso nello sviluppo
else
    err_approx_Cheby=0;
end

F=Costruzione_Matrice_Rappresentativa(N,vettcoeff,order);
end
