#include<iostream>
#include<vector>
#include<map>
#include<cmath>
#include<iomanip>
using namespace std;

class solution{
public:
    solution(double(*f)(double), double a, double b, double e){
        this->f = f;
        this->a = a;
        this->b = b;
        this->e = e;
    };
    //得到S(n)
    double get_S(int n){
        int m=n/2;
        double h=(b-a)/n;
        vector<double> x(n+1);
        for(int i=0;i<=n;i++)
        {
            x[i]=a+i*h;
        }
        double sum_odd=0,sum_even=0;
        for(int i=0;i<=m-1;i++){
            sum_odd += f(x[2*i+1]);
        }
        sum_odd *= 4;
        for(int i=1;i<=m-1;i++)
        {
            sum_even += f(x[2*i]);
        }
        sum_even *= 2;
        return h/3*(f(a)+f(b)+sum_odd+sum_even);
    }
    //得到H(n)的值
    double get_H(int n){
        double h_n=(b-a)/n;
        vector<double> x(n);
        for(int i=0;i<n;i++){
            x[i]=a+(i+0.5)*h_n;
        }
        double result=0;
        for(int i=0;i<n;i++)
        {
            result += f(x[i]);
        }
        result *= h_n;
        return result;
    }
    //得到积分的近似值
    double get_solution()
    {
        double S1=0,S2=0;
        int n=10,m=5;
        S2=get_S(n);
        S1 = S2 +100;
        while(fabs(S1-S2)>e)
        {
            S1=S2;
            double H_n=get_H(n);
            double H_2n=get_H(2*n);
            S2 = 0.5*S1 +(4*H_2n-H_n)/6;
            n*=2;
        }
        return S2;
    }

private:
    double (*f)(double);
    double a,b,e;
};

//test function
double f1(double x){
    return x*x;
}
double f2(double x){
    return x*x*x+f1(x);
}

int main()
{
    //test log(x)
    solution s1(log,1,2,1e-6);
    cout<<"test log(x) from 1 to 2:"<<endl;
    cout<<setprecision(6)<<s1.get_solution()<<endl;
    

    //test exp(x)
    solution s2(exp,1,2,1e-8);
    cout<<"test exp(x) from 1 to 2:"<<endl;
    cout<< setprecision(8)<<s2.get_solution()<<endl;

    //test x^2
    solution s3(f1,1,2,1e-8);
    cout<<"test x^2 from 1 to 2:"<<endl;
    cout<< setprecision(8)<<s3.get_solution()<<endl;

    //test x^3+x^2
    solution s4(f2,1,2,1e-8);
    cout<<"test x^3+x^2 from 1 to 2:"<<endl;
    cout<< setprecision(8)<<s4.get_solution()<<endl;

    return 0;
}