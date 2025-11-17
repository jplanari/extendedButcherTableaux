#include <string>
#include <vector>
#include <cstdint>
#include <map>
#include <functional>

/// Index mapping for strictly lower triangular storage of A

#define ID(eps, sig) ((eps - 1) * (eps - 2) / 2 + (sig - 1))

class ExtendedButcherTableau 
{
public:
  // --- Core data ---
  unsigned int s = 1;               // Number of stages
  unsigned int n_steps = 1;         // Number of steps
  std::string name;

  std::vector<double> b;            // Weights (size s)
  std::vector<double> A;            // Coefficients (size s*(s-1)/2), flattened
  std::vector<double> beta;         // Step weights (size n_steps)
  std::vector<double> alpha;        // 
                                    //
public:
  // --- Constructors ---
  ExtendedButcherTableau() = default;  

  ExtendedButcherTableau(unsigned int s,
      const std::vector<double>& b,
      const std::vector<double>& A,
      const std::string& name,
      unsigned int n_steps = 1,
      const std::vector<double>& beta = {},
      const std::vector<double>& alpha = {1.0,0.0});

  // --- Fill function ---

  void fill(unsigned int s,
      const std::vector<double>& b,
      const std::vector<double>& A,
      const std::string& name,
      unsigned int n_steps = 1,
      const std::vector<double>& beta = {},
      const std::vector<double>& alpha = {1.0,0.0});

  // --- Print function ---
  void print(void) const;

  // --- Create method ---

  static ExtendedButcherTableau create(
      const std::string& name,
      double param1 = 0.0,
      double param2 = 0.0);

  // --- Self-adaptive timestep methods ---
  
  double stabilityFunction(double r, double phi) const;
  double stabilityRadius(double phi) const;

private:

  // -- Function types ---
  using FunctionTypeVoid          = std::function<ExtendedButcherTableau(void)>;
  using FunctionTypeDouble        = std::function<ExtendedButcherTableau(double)>;
  using FunctionTypeDoubleDouble  = std::function<ExtendedButcherTableau(double, double)>;

  // --- Registries ---
  // Maps from string names to factory functions
  
  static std::map<std::string, FunctionTypeVoid> registryVoid;
  static std::map<std::string, FunctionTypeDouble> registryDouble;
  static std::map<std::string, FunctionTypeDoubleDouble> registryDoubleDouble;

    // --- Factory methods ---
  // Multi-step methods
  static ExtendedButcherTableau AB2(void);
  static ExtendedButcherTableau k1L2(double k);
  // Multi-stage methods
  static ExtendedButcherTableau ForwardEuler(void);
  static ExtendedButcherTableau paramRK2(double theta);
  static ExtendedButcherTableau heunRK2(void);
  static ExtendedButcherTableau midpointRK2(void);
  static ExtendedButcherTableau stdRK3(void);
  static ExtendedButcherTableau heunRK3(void);
  static ExtendedButcherTableau wrayRK3(void);
  static ExtendedButcherTableau nystromRK3(void);
  static ExtendedButcherTableau ralstonRK3(void);
  static ExtendedButcherTableau SSPRK3(void);
  static ExtendedButcherTableau PS3p5q(double c3);
  static ExtendedButcherTableau stdRK4(void);
  static ExtendedButcherTableau varRK4(void);

};
