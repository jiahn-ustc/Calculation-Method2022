n=160;

errors = 1e-10;
h=1/n;
A = zeros(n-1,n-1);
D = zeros(n-1,n-1);
I = eye(n-1,n-1);
b = zeros(n-1,1);

ue = zeros(n-1,1);
%get b,ue
for i=1:n-1
    b(i,1)= f(i*h);
    ue(i,1)= u(i*h);
end
%get A,D
for i=1:n-1
    
    A(i,i)=2/(h*h);
    D(i,i)=2/(h*h);
    if(i>1)
        A(i,i-1)=-1/(h*h);
    end
    if(i<n-1)
        A(i,i+1)=-1/(h*h);
    end
    
end
InvD = inv(D);
R = I-InvD*A;
g = InvD*b;
x1 = zeros(n-1,1);
x2 = x1+1;
num = 0;

while norm(x1-x2,inf)>errors
    x1 = x2;
    x2 = R*x1+g;
    num = num + 1;
end

fprintf("n=\n");
n
fprintf("x2=\n");
x2
fprintf("ue=\n");
ue
fprintf("eh=\n");
eh=norm(x2-ue,2);
eh
fprintf("迭代次数为：\n");
num



function result = f(x)
    result = pi*pi*sin(pi*x);
end

function result = u(x)
    result = sin(pi*x);
end