function[A,B]=fullmatrix_5eq(n,x,H,gamma,F,Da,r,bc,D2)
N=n;
R=r;

MM=zeros(4*N,4*N);
M=zeros(4*N,4*N);
I=eye(N,N);

A11=D2-x*I;
A12=-I;
A13=0*I;
A14=0*I;
A21=0*I;
A22=-I+Da*(D2-x*I);
A23=R*F*x;
A24=0*I;
A31=-R*I;
A32=0*I;
A33=D2-x*I-H*I;
A34=H*I;
A41=0*I;
A42=0*I;
A43=H*gamma*I;
A44=D2-x*I-H*gamma*I;


A=1;
B11=0*I;
B22=0*I;
B33=I;
B44=A*I;


%Costruzione matrici a blocchi
for i=1:N-2
    for j=1:N
        MM(i,j)= A11(i,j);
    end
    for j=N+1:2*N
        MM(i,j)= A12(i,j-N);
    end
    for j=2*N+1:3*N
        MM(i,j)=A13(i,j-2*N);
    end
    for j=3*N+1:4*N
        MM(i,j)=A14(i,j-3*N);
    end
end
for i=N+1:2*N-2             
    for j=1:N             
        MM(i,j)=A21(i-N,j); 
    end
    for j=N+1:2*N
        MM(i,j)=A22(i-N,j-N);
    end
    for j=2*N+1:3*N
        MM(i,j)=A23(i-N,j-2*N);
    end
    for j=3*N+1:4*N
        MM(i,j)=A24(i-N,j-3*N);
    end
end
for i=2*N+1:3*N-2
    for j=1:N
        MM(i,j)=A31(i-2*N,j);
    end
    for j=N+1:2*N
        MM(i,j)=A32(i-2*N,j-N);
    end
    for j=2*N+1:3*N
        MM(i,j)=A33(i-2*N,j-2*N);
    end
    for j=3*N+1:4*N
        MM(i,j)=A34(i-2*N,j-3*N);
    end
end
for i=3*N+1:4*N-2
    for j=1:N
        MM(i,j)=A41(i-3*N,j);
    end
    for j=N+1:2*N
        MM(i,j)=A42(i-3*N,j-N);
    end
    for j=2*N+1:3*N
        MM(i,j)=A43(i-3*N,j-2*N);
    end
    for j=3*N+1:4*N
        MM(i,j)=A44(i-3*N,j-3*N);
    end
end


for i=1:N-2
    for j=1:N
        M(i,j)= B11(i,j);
    end
end
for i=N+1:2*N-2
    for j=N+1:2*N
        M(i,j)=B22(i-N,j-N);
    end
end
for i=2*N+1:3*N-2
    for j=2*N+1:3*N
        M(i,j)=B33(i-2*N,j-2*N);
    end
end
for i=3*N+1:4*N-2
    for j=3*N+1:4*N
        M(i,j)=B44(i-3*N,j-3*N);
    end
end


% Righe condizioni al bordo
MM(N-1,1:N)=bc(1,:);
MM(N,1:N)=bc(2,:); 
MM(2*N-1,N+1:2*N)=bc(1,:);
MM(2*N,N+1:2*N)=bc(2,:);
MM(3*N-1,2*N+1:3*N)=bc(1,:);
MM(3*N,2*N+1:3*N)=bc(2,:);
MM(4*N-1,3*N+1:4*N)=bc(1,:);
MM(4*N,3*N+1:4*N)=bc(2,:);


% Sostituzione condizione al bordo
index=[N-1 N 2*N-1 2*N 3*N-1 3*N 4*N-1 4*N  ];
for h=1:length(index)
    if h==1 || h==2
        j=1;
        i=2;
            m=MM(index(i),index(j))/MM(index(j),index(j));
            MM(index(i),:)=MM(index(i),:)-m*MM(index(j),:);
    end
    if h==3 || h==4
        j=3;
        i=4;
            m=MM(index(i),index(j))/MM(index(j),index(j));
            MM(index(i),:)=MM(index(i),:)-m*MM(index(j),:);
    end
    if h==5 || h==6
        j=5;
        i=6;
            m=MM(index(i),index(j))/MM(index(j),index(j));
            MM(index(i),:)=MM(index(i),:)-m*MM(index(j),:);
    end
    if h==7 || h==8
        j=7;
        i=8;
            m=MM(index(i),index(j))/MM(index(j),index(j));
            MM(index(i),:)=MM(index(i),:)-m*MM(index(j),:);
    end
%     if h==9 || h==10
%         j=9;
%         i=10;
%             m=MM(index(i),index(j))/MM(index(j),index(j));
%             MM(index(i),:)=MM(index(i),:)-m*MM(index(j),:);
%     end
%     if h==11 || h==12
%         j=11;
%         i=12;
%             m=MM(index(i),index(j))/MM(index(j),index(j));
%             MM(index(i),:)=MM(index(i),:)-m*MM(index(j),:);
%     end
for k=1:length(MM)
    if k~=index(j) && k~=index(j)+1
    MM(k,index(i)-N+1:index(j)-1)=MM(k,index(i)-N+1:index(j)-1)-MM(index(i),index(i)-N+1:index(j)-1)/MM(index(i),index(i))*MM(k,index(i));
    M(k,index(i)-N+1:index(j)-1)=M(k,index(i)-N+1:index(j)-1)-MM(index(i),index(i)-N+1:index(j)-1)/MM(index(i),index(i))*M(k,index(i));
    end
end

MM(N-1,1:N)=bc(1,:);
MM(N,1:N)=bc(2,:); 
MM(2*N-1,N+1:2*N)=bc(1,:);
MM(2*N,N+1:2*N)=bc(2,:);
MM(3*N-1,2*N+1:3*N)=bc(1,:);
MM(3*N,2*N+1:3*N)=bc(2,:);
MM(4*N-1,3*N+1:4*N)=bc(1,:);
MM(4*N,3*N+1:4*N)=bc(2,:);


MM(:,[index(h),index(i)])=MM(:,[index(i),index(h)]);
MM([index(h),index(i)],:)=MM([index(i),index(h)],:);
M(:,[index(h),index(i)])=M(:,[index(i),index(h)]);
M([index(h),index(i)],:)=M([index(i),index(h)],:);
end

%index_new=[index,2*N-1,2*N];
index_new=index;
MM(index_new,:)=[];
MM(:,index_new)=[];
M(index_new,:)=[];
M(:,index_new)=[];


A=MM;
B=M;

end




