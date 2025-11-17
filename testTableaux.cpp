#include "ExtendedButcherTableau.h"
#include <iostream>
#include <cmath>

int main()
{
  auto rk2 = ExtendedButcherTableau::create("paramRK2", 0.5);

  std::cout << "RK2 Butcher Tableaux:" << std::endl;
  rk2.print();

  auto ab2 = ExtendedButcherTableau::create("AB2");
  std::cout << "AB2 Butcher Tableaux:" << std::endl;
  ab2.print();
  std::cout << "Radius at phi=1.0 rad: " << ab2.stabilityRadius(1.0) << std::endl;
  std::cout << "Radius at phi=0.0: " << ab2.stabilityRadius(0.0) << std::endl;

  auto k1l2 = ExtendedButcherTableau::create("k1L2",0.5); //Should be equivalent to AB2
  std::cout << "k1L2 Butcher Tableaux:" << std::endl;
  k1l2.print();

 
  auto rk3 = ExtendedButcherTableau::create("Heun RK3");
  std::cout << "Heun RK3 Butcher Tableaux:" << std::endl;
  rk3.print();

  std::cout << "Radius at phi=0.0: " << rk3.stabilityRadius(0.0) << std::endl;
  std::cout << "Radius at phi=pi/2: " << rk3.stabilityRadius(M_PI/2) << std::endl;

  return 0;
}
