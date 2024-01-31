%costruzione della matrice D^2
%D^2 è matrice delle derivate seconde
function[matrixout]=chder2(rows,columns)
matrix=zeros(rows,columns);

% coefficiente dovuto alla discretizzazione dell'operatore derivata
ck = 2;
for i= 1:rows
for  j = i + 2: 2: columns
matrix(i, j) = ((j-1) .* ((j-1).^2 - (i-1).^2)) / ck;
end
ck=1;
end
%j = fix(columns+ 2);
matrixout=matrix;
end

