n=160;
errors = 1e-8;
h=1/n;
b = zeros(n-1,1);
ue = zeros(n-1,1);
delta_u = zeros(n-1,1);

%get b,ue
for i=1:n-1
    b(i,1)= f(i*h);
    ue(i,1)= u(i*h);
end

u1 = zeros(n-1,1);
u2 = u1+1;
B = zeros(n-1,1);
A = zeros(n-1,n-1);
num=0;
while norm(u1-u2,inf)>errors
    u1 = u2;
    for i=1:n-1
        A(i,i) = (2+3*h^2*u2(i,1)^2);
        
        if(i==1)
            A(i,i+1) = -1;
            B(i,1) = -(2*u2(1,1)-u2(2,1)+h^2*(u2(1,1)^3-b(1,1)));
        end
        if( i>=2 && i<= n-2)
            A(i,i+1) = -1;
            A(i,i-1) = -1;
            B(i,1) = -(2*u2(i,1)-u2(i-1,1)-u2(i+1,1)+h^2*(u2(i,1)^3-b(i,1)));
        end
        if(i==n-1)
            A(i,i-1) = -1;
            B(i,1) = -(2*u2(n-1,1)-u2(n-2,1)+h^2*(u2(n-1,1)^3-b(n-1,1)));
        end
    end
    delta_u = A\B;
    u2 = u1+delta_u;
    num = num + 1;
end
fprintf("n=\n");
n
fprintf("u2=\n");
u2
fprintf("ue=\n");
ue
fprintf("eh=\n");
eh=norm(u2-ue,2);
eh
fprintf("迭代次数为：\n");
num




function result = f(x)
    result = pi*pi*sin(pi*x)+sin(pi*x).^3;
end

function result = u(x)
    result = sin(pi*x);
end