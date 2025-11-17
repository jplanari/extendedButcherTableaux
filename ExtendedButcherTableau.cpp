#include "ExtendedButcherTableau.h"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <map>
#include <functional>
#include <complex>

ExtendedButcherTableau::ExtendedButcherTableau(unsigned int s,
    const std::vector<double>& b,
    const std::vector<double>& A,
    const std::string& name,
    unsigned int n_steps,
    const std::vector<double>& beta,
    const std::vector<double>& alpha)
{
  fill(s, b, A, name, n_steps, beta, alpha);
}

void ExtendedButcherTableau::fill(unsigned int s, const std::vector<double>& b,
    const std::vector<double>& A, const std::string& name, unsigned int n_steps,
    const std::vector<double>& beta, const std::vector<double>& alpha)
{
  this->s = s;
  this->name = name;
  this->n_steps = n_steps;

  // resize vectors
  this->b = b;
  this->A = A;
  this->A.resize(s * (s - 1) / 2); // ensure correct size
  
  this->beta = beta;
  if (!beta.empty())
    this->beta.resize(n_steps); // ensure correct size
  this->alpha = alpha;
}

void ExtendedButcherTableau::print(void) const
{
  std::cout << "Butcher Tableau: " << name << std::endl;
  std::cout << "Number of stages: " << s << std::endl;
  std::cout << "b: ";
  for (const auto& bi : b)
    std::cout << bi << " ";
  std::cout << std::endl;

  std::cout << "A: " << std::endl;
  for (const auto&ai : A)
    std::cout << ai << " ";
  std::cout << std::endl;

  std::cout << "Number of steps: " << n_steps << std::endl;
  std::cout << "beta: ";
  for (const auto& betai : beta)
    std::cout << betai << " ";
  std::cout << std::endl;

  std::cout << "alpha: ";
  for (const auto& alphai : alpha)
    std::cout << alphai << " ";
  std::cout << std::endl;
}

double bisection(std::function<double(double)> f,
               double r_min = 0.01,
               double r_max = 5.0,
               double tol = 1e-8)
{
  double a = r_min, b = r_max;
  double fa = f(a), fb = f(b);
  if (fa * fb > 0)
  {
    throw std::runtime_error("Root not bracketed in bisection method.");
  }

  while ((b - a) > tol)
  {
    double c = 0.5 * (a + b);
    double fc = f(c);
    if (fa * fc < 0){
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }

  return 0.5 * (a + b);
}

double factorial(unsigned int x)
{
  if (x == 0 || x == 1)
    return 1.0;
  else
    return x * factorial(x - 1);
}

double ExtendedButcherTableau::stabilityFunction(double r, double phi) const
{
  std::complex<double> z = std::polar(r, M_PI - phi);
  if (n_steps > 1)
  {
    if (beta.size() != 2)
    {
      throw std::runtime_error("Stability function for multi-step methods only implemented for AB2.");
    }
    std::complex<double> A = 1.0 + 1.5 * z;
    std::complex<double> B = -0.5*z;
    std::complex<double> sqrt_part = std::sqrt(A*A + 4.0*B);
    return std::max(std::abs((A + sqrt_part)/2.0), std::abs((A - sqrt_part)/2.0));
  }
  else
  {
    if(b.size() > 4)
    {
      throw std::runtime_error("Stability function for multi-stage methods only implemented for up to 4 stages.");
    }
    std::complex<double> sum = 1.0;
    double coef, ii;
    for (long unsigned int i = 1; i <= b.size(); ++i)
    {
      coef = 1/factorial(i);
      ii = static_cast<double>(i);
      sum += coef * std::pow(z, ii);
    }
    return std::abs(sum);
  }
}

double ExtendedButcherTableau::stabilityRadius(double phi) const
{
  auto stabilityFunc = [this, phi](double r) {
    return this->stabilityFunction(r, phi) - 1.0;
  };

  return bisection(stabilityFunc);
}

// ------------------------------------------------------

std::map<std::string, ExtendedButcherTableau::FunctionTypeVoid>
    ExtendedButcherTableau::registryVoid = {
    {"AB2", ExtendedButcherTableau::AB2},
    {"Forward Euler", ExtendedButcherTableau::ForwardEuler},
    {"Heun RK2", ExtendedButcherTableau::heunRK2},
    {"Midpoint RK2", ExtendedButcherTableau::midpointRK2},
    {"Standard RK3", ExtendedButcherTableau::stdRK3},
    {"Heun RK3", ExtendedButcherTableau::heunRK3},
    {"Wray RK3", ExtendedButcherTableau::wrayRK3},
    {"Nystrom RK3", ExtendedButcherTableau::nystromRK3},
    {"Ralston RK3", ExtendedButcherTableau::ralstonRK3},
    {"SSPRK3", ExtendedButcherTableau::SSPRK3},
    {"Standard RK4", ExtendedButcherTableau::stdRK4},
    {"Kutta RK4", ExtendedButcherTableau::varRK4}
    };

std::map<std::string, ExtendedButcherTableau::FunctionTypeDouble>
    ExtendedButcherTableau::registryDouble = {
    {"k1L2", ExtendedButcherTableau::k1L2},
    {"paramRK2", ExtendedButcherTableau::paramRK2},
    {"3p5q", ExtendedButcherTableau::PS3p5q}
    };

std::map<std::string, ExtendedButcherTableau::FunctionTypeDoubleDouble>
    ExtendedButcherTableau::registryDoubleDouble = {};

ExtendedButcherTableau ExtendedButcherTableau::create(
    const std::string& name,
    double param1,
    double param2)
{
    // --- 0 parameters ---
    auto it0 = registryVoid.find(name);
    if (it0 != registryVoid.end())
        return (it0->second)();

    // --- 1 parameter ---
    auto it1 = registryDouble.find(name);
    if (it1 != registryDouble.end())
        return (it1->second)(param1);

    // --- 2 parameters ---
    auto it2 = registryDoubleDouble.find(name);
    if (it2 != registryDoubleDouble.end())
        return (it2->second)(param1, param2);

    throw std::runtime_error("Unknown integration scheme: " + name);
}

// ------------------------------------------------------
// Factory methods
// ------------------------------------------------------

//MULTI-STEP METHODS
// ------------------------------------------------------

ExtendedButcherTableau ExtendedButcherTableau::AB2(void)
{
  std::vector<double> b = {1.0};
  std::vector<double> A = {}; // flattened lower triangular
  std::vector<double> beta = {1.5, -0.5};
  std::string name = "AB2";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();

  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::k1L2(double k)
{
  std::vector<double> b = {1.0};
  std::vector<double> A = {}; // flattened lower triangular
  std::vector<double> beta = {(1.0 + k)/(k+0.5), -k(k+0.5)};
  std::vector<double> alpha = {2*k/(k+0.5), -(k-0.5)/(k+0.5)};
  std::string name = "k1L2(" + std::to_string(k) + ")";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();

  return ExtendedButcherTableau(s, b, A, name, n_steps, beta, alpha);
}

//MULTI-STAGE METHODS
// ------------------------------------------------------

ExtendedButcherTableau ExtendedButcherTableau::ForwardEuler(void)
{
  std::vector<double> b = {1.0};
  std::vector<double> A = {}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "Forward Euler";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();

  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::paramRK2(double theta)
{
  std::vector<double> b = {1 - theta, theta};
  std::vector<double> A = {1.0 / (2.0 * theta)}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "paramRK2";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();

  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::heunRK2(void)
{
  ExtendedButcherTableau tableau = paramRK2(0.5);
  tableau.name = "Heun RK2";
  return tableau;
}

ExtendedButcherTableau ExtendedButcherTableau::midpointRK2(void)
{
  ExtendedButcherTableau tableau = paramRK2(1.0);
  tableau.name = "Midpoint RK2";
  return tableau;
}

ExtendedButcherTableau ExtendedButcherTableau::stdRK3(void)
{
  std::vector<double> b = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
  std::vector<double> A = {0.5,  -1.0, 2.0}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "Standard RK3";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();

  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::heunRK3(void)
{
  std::vector<double> b = {0.25, 0.0, 0.75};
  std::vector<double> A = {1.0/3.0, 0.0, 2.0/3.0}; // flattened lower triangular
  std::vector<double> beta = {1.0};                                                   
  std::string name = "Heun RK3";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();

  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::wrayRK3(void)
{
  std::vector<double> b = {0.25, 0.0, 0.75};
  std::vector<double> A = {8.0/15.0, 0.25, 5.0/12.0}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "Wray RK3";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();
  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::nystromRK3(void)
{
  std::vector<double> b = {0.25, 3.0/8.0, 3.0/8.0};
  std::vector<double> A = {2.0/3.0, 0.0, 2.0/3.0}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "Nystrom RK3";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();
  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::ralstonRK3(void)
{
  std::vector<double> b = {2.0/9.0, 1.0/3.0, 4.0/9.0};
  std::vector<double> A = {0.5, 0.0, 2.0/3.0}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "Ralston RK3";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();
  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::SSPRK3(void)
{
  std::vector<double> b = {1.0/6.0, 1.0/6.0, 2.0/3.0};
  std::vector<double> A = {1.0, 0.25, 0.25}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "SSPRK3";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();
  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::PS3p5q(double c3)
{
  double a21 = (c3-1)/(4*c3-3);
  double a32 = (2*c3-1)/(2*a21);
  double a31 = c3-a32;
  double a42 = 1./3.0*((4*c3-3)/(a21*(2*c3-1))+0.5*pow(c3,2.0)/((c3-1)*(2*c3-1)));
  double a43 = a21/(2*c3-1);
  double a41 = 1-a42-a43;
  std::vector<double> A = {a21,a31,a32,a41,a42,a43};
  double b1 = 1/(12*(c3-1));
  double b2 = b1*pow(4*c3-3,2)/(2*c3-1);
  double b3 = -b1/(2*c3-1);
  double b4 = b1*(4*c3-3);
  std::vector<double> b = {b1,b2,b3,b4};
  std::vector<double> beta = {1.0};
  std::string name = "3p5q(" + std::to_string(c3) + ")";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();
  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::stdRK4(void)
{
  std::vector<double> b = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
  std::vector<double> A = {0.5,
                           0.0, 0.5,
                           0.0, 0.0, 1.0}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "Standard RK4";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();

  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}

ExtendedButcherTableau ExtendedButcherTableau::varRK4(void)
{
  std::vector<double> b = {1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0};
  std::vector<double> A = {1.0/3.0,
                           -1.0/3.0, 1.0,
                           1.0, -1.0, 1.0}; // flattened lower triangular
  std::vector<double> beta = {1.0};
  std::string name = "Kutta RK4";
  unsigned int s = b.size();
  unsigned int n_steps = beta.size();
  return ExtendedButcherTableau(s, b, A, name, n_steps, beta);
}



