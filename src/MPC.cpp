#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <iostream>

using CppAD::AD;


// TODO: Set the timestep length and duration
//size_t N = 0;
//double dt = 0;
const size_t N = 10; // long — in types.h:94
const double dt = 0.1; //0.1 ;
const double psiDesired = 0;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
double referenceVelocity = 10;


// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
const int x_start = 0;
const int y_start = x_start + N;
const int psi_start = y_start + N;
const int v_start = psi_start + N;
const int cte_start = v_start + N;
const int epsi_start = cte_start + N;
const int delta_start = epsi_start + N;
const int a_start = delta_start + N - 1;

const CppAD::AD<double> fOfX(const Eigen::VectorXd thePolynomialCoefficients, const CppAD::AD<double> theValueX) {
  //const int numberOfCoefficients=sizeof(thePolynomialCoefficients)/sizeof(thePolynomialCoefficients[0]);
  const int numberOfCoefficients = thePolynomialCoefficients.size();
  CppAD::AD<double> valueAtX=thePolynomialCoefficients[0];
  for (int c=1;c<numberOfCoefficients;c++) {
    valueAtX+=thePolynomialCoefficients[c]*CppAD::pow(theValueX,c);
  }
  //cout << "CppADfOfX-numberOfCoefficients:" << numberOfCoefficients << ", valueAtX:" << valueAtX << std::endl;
  return valueAtX;
}

const CppAD::AD<double> fPrimeOfX(const Eigen::VectorXd thePolynomialCoefficients, const CppAD::AD<double> theValueX) {
  //const int numberOfCoefficients=sizeof(thePolynomialCoefficients)/sizeof(thePolynomialCoefficients[0]);
  const int numberOfCoefficients = thePolynomialCoefficients.size();
  if (numberOfCoefficients > 0) {
    CppAD::AD<double> valueAtX=thePolynomialCoefficients[1];
    for (int c=2;c<numberOfCoefficients;c++) {
      CppAD::AD<double> power=c;
      CppAD::AD<double> coefficient=thePolynomialCoefficients[c];
      valueAtX+=power*coefficient*CppAD::pow(theValueX,c-1);
    }
    //cout << "CppADfPrimeOfX-numberOfCoefficients:" << numberOfCoefficients << ", valueAtX:" << valueAtX << std::endl;
    return valueAtX;
  } else {
    return 0.;
  }
}

const CppAD::AD<double> tangentialAngle(const Eigen::VectorXd thePolynomialCoefficients, const CppAD::AD<double> theValueX) {
  return CppAD::atan(fPrimeOfX(thePolynomialCoefficients, theValueX));
}

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    
    CppAD::AD<double> cteCost=0;
    for (int t = 0; t < N; t++) {
      //const CppAD::AD<double> yState = vars[y_start+t];
      //const CppAD::AD<double> yPath = polyeval(coeffs, t);
      //cout << "yState[" << t << "]:" << yState << ", yPath[" << t << "]:" << yPath << std::endl;
      //cteCost += pow(yState-yPath, 2); // minimize path errors in y (cte predicted)
      const CppAD::AD<double> cteError = vars[cte_start+t];
      cteCost += CppAD::pow(cteError, 2); // minimize path errors in y (cte predicted)
    }
    
    CppAD::AD<double> vCost=0;
    for (int t = 0; t < N; t++) {
      const CppAD::AD<double> vState = vars[v_start+t];
      //cout << "vState[" << t << "]:" << vState << std::endl;
      vCost += CppAD::pow(referenceVelocity-vState, 2); // minimize path errors in y (cte predicted)
    }
    
    CppAD::AD<double> psiCost=0;
    for (int t = 0; t < N; t++) {
      const CppAD::AD<double> psiState = vars[psi_start+t];
      //cout << "psiState[" << t << "]:" << psiState << std::endl;
      const CppAD::AD<double> psiError=vars[epsi_start+t];
      psiCost += CppAD::pow(psiError, 2); // minimize path errors in y (cte predicted)
    }
    
    CppAD::AD<double> daCost=0.;
    for (int t = 0; t < N-1; t++) {// control signals have one less entry
      //cost += pow(delta[t], 2) // minimize large control changes
      daCost += CppAD::pow(vars[delta_start+t], 2); // minimize control changes
      daCost += CppAD::pow(vars[a_start+t], 2); // minimize control changes
    }
    
    CppAD::AD<double> daDeltaCost=0.;
    //for (int t = 0; t < N-1; t++) {
    for (int t = 0; t < N-2; t++) {// indexes N+1
      //cost += pow(delta[t+1] - delta[t], 2) // where control change between control steps
      daDeltaCost += CppAD::pow(vars[delta_start+t+1] - vars[delta_start+t], 2); // where control change between control steps
      //cost += pow(a[t+1] - a[t], 2)
      daDeltaCost += CppAD::pow(vars[a_start+t+1] - vars[a_start+t], 2);
    }
    
    cout << "yCost:" << cteCost << ", psiCost:" << psiCost << ", dCost:" << daCost << ", daCost:" << daDeltaCost << ", vCost:" << vCost << std::endl;
    //fg[0] = 0;
    fg[0] = cteCost+psiCost+daCost+daDeltaCost+vCost;
    cout << "total cost:" << fg[0] << std::endl;
    
    // Reference State Cost
    // TODO: Define the cost related the reference state??? and
    // any anything you think may be beneficial.
    
    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.
    
    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];
    
    // The rest of the constraints
    for (int t = 1; t < N; t++) {
      AD<double> x1 = vars[x_start + t];
      
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      
      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.
      
      // TODO: Setup the rest of the model constraints
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      
      // y(t+1)=y(t)+v(t)*sin(phi(t))*dt
      const AD<double> y1 = vars[y_start + t];
      const AD<double> y0 = vars[y_start + t - 1];
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      //cout << "fg.y[" << t << "]:" << fg[1 + y_start + t] << std::endl;
      
      
      // ψ(t+1)=ψ(t)+v(t)/Lf*δ(t)*dt
      const AD<double> psi1 = vars[psi_start + t];
      //const AD<double> psi0 = vars[psi_start + t - 1];
      const AD<double> d0 = vars[delta_start + t -1];
      //fg[1 + psi_start + t] = psi1-(psi0+v0/Lf*d0*dt);
      fg[1 + psi_start + t] = psi1-(psi0-v0/Lf*d0*dt);// project notes: δ is positive rotating counter clockwise

      
      // v(t+1) = v(t)+a(t)*dt
      const AD<double> v1 = vars[v_start + t];
      const AD<double> a0 = vars[a_start + t - 1];
      fg[1 + v_start + t] = v1-(v0+a0*dt);
      
      // cte(t+1) = f(x(t)) - y(t) + v(t)*sin(eψ(t))*dt
      // cte(t+1) = cte(t)         + v(t)*sin(eψ(t))*dt
      // eψ(t)=ψ(t)−ψdes(t)
      const AD<double> cte1 = vars[cte_start + t];
      const AD<double> epsi0 = vars[epsi_start + t - 1];
      //const AD<double> cte0 = vars[cte_start + t - 1];
      const AD<double> cte0 = fOfX(coeffs, x0)-y0;
      fg[1 + cte_start + t] = cte1-(cte0 + v0*CppAD::sin(epsi0)*dt);
      
      // eψ(t+1)=eψ(t)+(v(t)/Lf)*δ(t)*dt
      // eψ(t)=ψ(t)−ψdes(t)
      const AD<double> epsi1 = vars[epsi_start + t];
      const AD<double> psiDesiredAtX0 = tangentialAngle(coeffs, x0);
      const AD<double> psiError = psi0-psiDesiredAtX0;
      fg[1 + epsi_start + t] = epsi1-(psiError + (v0/Lf)*d0*dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

//vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
vector<double> MPC::Solve(Eigen::VectorXd state/*x,y,psi,v,cte,epsi*/, Eigen::VectorXd coeffs/*=polyfit(ptsx,ptsy,1)*/) {

  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  //size_t n_vars = 0;
  // TODO: Set the number of constraints
  //size_t n_constraints = 0;

  const double x = state[0];
  const double y = state[1];
  const double psi = state[2];
  const double v = state[3];
  const double cte = state[4];
  const double epsi = state[5];
  
  // number of independent variables
  // N timesteps == N - 1 actuations
  // number of time steps * number of states + (number of time steps-1) * number of actuators
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Number of constraints: number of time steps * number of states
  size_t n_constraints = N * 6;
  
  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;
  
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  
  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;
  
  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;
  
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);// constructor -> saves off coeffs

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {solution.x[x_start + 1],   solution.x[y_start + 1],
    solution.x[psi_start + 1], solution.x[v_start + 1],
    solution.x[cte_start + 1], solution.x[epsi_start + 1],
    solution.x[delta_start],   solution.x[a_start]};
}
