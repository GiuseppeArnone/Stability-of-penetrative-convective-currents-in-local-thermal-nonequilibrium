function[coeff]=CoeffExtraction(t,ORDER)
% function per ottenere un vettore di coefficienti di un polinomio
% t � una funzione simbolica polinomiale
% ORDER � l'ordine del polinomio
syms z
for i = 1:ORDER
    coeff(i)=subs(t,z,0);
    syms z
    t(z)=simplify((t-coeff(i))/z);
end
end