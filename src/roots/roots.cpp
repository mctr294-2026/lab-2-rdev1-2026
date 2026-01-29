#include "roots.hpp"

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) 
{
    if (!root) return false;

    const int MAX = 1000;
    double fa = f(a);
    double fb = f(b);

    
    if (fa == 0)
    {
        *root=a; 
        return true;
    }

    if (fb == 0)
    {
        *root=b;
        return true;
    }
    

    if (fa * fb > 0) return false;

      for (int i = 0; i < MAX; ++i)
      {
        double c = (a + b) / 2.0;
        double fc = f(c);

        //check if midpoint is the root
        if (fc == 0)
        {
            *root = c;
            return true;
        }

        //decide the side to repeat the steps
        if (fc * fa < 0)
        {
            b = c;
            fb = fc;
        }
        
        else
        {
            a = c;
            fa = fc;
        }

      }

      *root = (a + b) / 2.0;

      return true;
    }

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) 
{
    if (!root) return false;
    
    const int MAX = 10000000;
    double fa = f(a);
    double fb = f(b);

    if (fa == 0)
    {
        *root=a; 
        return true;
    }

    if (fb == 0)
    {
        *root=b;
        return true;
    }

    if (fa * fb > 0) return false;

    double c = a;
    double last_c = c;
    
    for (int i = 0; i < MAX; ++i)
    {
        double denom = (fb - fa);
        if (denom == 0) return false;

        c = ((a * fb)-(b * fa)) / denom;

        double fc = f(c);

        if (fc == 0)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }

        else
        {
            a = c;
            fa = fc;
        }

    last_c = c;
    }

    *root = c;

    return true;


}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) 

    {
    if (!root) return false;
    
    const int MAX = 1000;
    double fa = f(a);
    double fb = f(b);

    if (c < a || c > b) return false; 
    
    double x = c;

    for (int i = 0; i < MAX; ++i)
    {
        double fx = f(x);
        if (fx == 0)
        {
            *root = x;
            return true;
        }

        double gx = g(x);
        if (gx == 0) return false;

        double x2 = x - fx / gx;

        if (x2 < a || x2 > b) return false;

        if (x2 == 0)
        {
            *root = x2;
            return true;
        }

        x = x2;

    }

    return false;

}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) 
{

if (!root) return false;
    
    const int MAX = 1000;
    double fa = f(a);
    double fb = f(b);

    if ((c < a) || c > b) return false; 
    
    double x0 = c;
    double x1 = (std::abs(c - a) > std::abs(c - b)) ? a : b;
    
    double f0 = f(x0);
    double f1 = f(x1);

    if (f0 == 0) {*root = x0; return true;}
    if (f1 == 0) {*root = x1; return true;}

    for (int i = 0; i < MAX; ++i)
    {
        double denom = (f1 - f0);
        
        if (denom == 0) return false;

        double x2 = x1 - f1 * (x1 - x0) / denom;

        if (x2 < a || x2 > b) return false;
      
        double f2 = f(x2);

        if (f2 == 0 || (x2 - x1) == 0)
        {
            *root = x2;
            return true;

        }
        
    x0 = x1;
    f0 = f1;
    x1 = x2;
    f1 = f2;

    }

    return false;
}