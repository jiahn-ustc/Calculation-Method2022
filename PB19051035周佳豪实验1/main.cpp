#include <iostream>
#include <map>
#include <algorithm>
using namespace std;

class Solution
{
public:
    Solution()
    {
        f[1920] = 105711;
        f[1930] = 123203;
        f[1940] = 131669;
        f[1950] = 150697;
        f[1960] = 179323;
        f[1970] = 203212;
        n = 5;
        x = new double[n + 1];
        x[0] = 1920;
        x[1] = 1930;
        x[2] = 1940;
        x[3] = 1950;
        x[4] = 1960;
        x[5] = 1970;
        g = new double[n + 1];
        get_g();
    }
    double L(int year)
    {
        double result = 0;

        for (int i = 0; i <= n; i++)
        {
            double temp = 1;
            for (int j = 0; j <= i - 1; j++)
            {
                // cout<<"year-x[j]="<<year-x[j]<<endl;
                //  cout<<"x[i]-x[j]="<<x[i]-x[j]<<endl;
                temp *= (year - x[j]) / (x[i] - x[j]);
            }

            for (int j = i + 1; j <= n; j++)
            {
                // cout<<"year-x[j]="<<year-x[j]<<endl;
                // cout<<"x[i]-x[j]="<<x[i]-x[j]<<endl;
                temp *= (year - x[j]) / (x[i] - x[j]);
            }

            /*
            cout<<"temp = "<<temp<<endl;
            cout<<"f[x[i]] = "<<f[x[i]]<<endl;
            cout<<"temp*f[x[i]] = "<<temp*f[x[i]]<<endl;*/
            result += temp * f[x[i]];
        }
        return result;
    }
    double _L(int year)
    {
        x[0] = 1910;
        f[1910] = 91772;
        double result = L(year);
        x[0] = 1920;
        return result;
    }
    double get_error(int year)
    {
        return (year - 1920) / (1920 - 1910) * (L(year) - _L(year));
    }
    //得到函数g
    void get_g()
    {
        map<string, double> temp;
        for (int i = 0; i <= n; i++)
        {
            temp[to_string(i)] = f[x[i]];
            // cout<<"temp["<<i<<"] = "<<temp[to_string(i)]<<endl;
        }
        for (int i = 2; i <= n + 1; i++)
        {
            for (int j = 0; j <= n + 1 - i; j++)
            {
                string s, s1, s2;
                for (int k = j; k < j + i; k++)
                {
                    s += to_string(k);
                }
                for (int k = j + 1; k < j + i; k++)
                    s1 += to_string(k);
                for (int k = j; k < j + i - 1; k++)
                    s2 += to_string(k);
                temp[s] = (temp[s1] - temp[s2]) / (x[j + i - 1] - x[j]);

                // cout<<"temp["<<s<<"] = "<<temp[s]<<"=(temp["<<s1<<"]-temp["<<s2<<"])/("<<x[j+i-1]<<"-"<<x[j]<<")"<<endl;
            }
        }
        string str;
        for (int i = 0; i <= n; i++)
        {
            str += to_string(i);
            g[i] = temp[str];
        }
    }
    double N(int year)
    {
        get_g();
        double t = 1, newton = g[0];
        for (int k = 1; k <= n; k++)
        {
            t *= (year - x[k - 1]);
            newton += t * g[k];
        }
        return newton;
    }

private:
    double *x;
    double *g;
    map<int, double> f;
    int n;
};

int main()
{
    Solution s;
    // Lagrange interpolation
    cout << "----test Lagrange interpolation:----" << endl;
    cout << "1910: " << s.L(1910) << endl;
    cout << "1965: " << s.L(1965) << endl;
    cout << "2002: " << s.L(2002) << endl;
    cout << "---test accuracy:---" << endl;
    cout << "1965:" << s.get_error(1965) << endl;
    cout << "2002:" << s.get_error(2002) << endl;
    // Newton interpolation
    cout << "----test Newton interpolation:----" << endl;
    cout << "1965: " << s.N(1965) << endl;
    cout << "2012: " << s.N(2012) << endl;

    return 0;
}