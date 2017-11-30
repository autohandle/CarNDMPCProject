# Model Predictive Control (MPC)

#### Compiling
##### Code must compile without errors with cmake and make.

``` shell
Softwares-MacBook-Pro:tmp david$ git clone http://github.com/autohandle/CarNDMPCProject.git
Cloning into 'CarNDMPCProject'...
warning: redirecting to https://github.com/autohandle/CarNDMPCProject.git/
remote: Counting objects: 1982, done.
remote: Compressing objects: 100% (1469/1469), done.
remote: Total 1982 (delta 457), reused 1982 (delta 457), pack-reused 0
Receiving objects: 100% (1982/1982), 2.52 MiB | 2.69 MiB/s, done.
Resolving deltas: 100% (457/457), done.
Softwares-MacBook-Pro:tmp david$ cd CarNDMPCProject/
Softwares-MacBook-Pro:CarNDMPCProject david$ mkdir build
Softwares-MacBook-Pro:CarNDMPCProject david$ cd build
Softwares-MacBook-Pro:build david$ cmake ..
-- The C compiler identification is AppleClang 9.0.0.9000037
-- The CXX compiler identification is AppleClang 9.0.0.9000037
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Configuring done
-- Generating done
-- Build files have been written to: /tmp/CarNDMPCProject/build
Softwares-MacBook-Pro:build david$ make
Scanning dependencies of target mpc
[ 33%] Building CXX object CMakeFiles/mpc.dir/src/MPC.cpp.o
/tmp/CarNDMPCProject/src/MPC.cpp:229:10: warning: unused variable 'i' [-Wunused-variable]
  size_t i;
         ^
/tmp/CarNDMPCProject/src/MPC.cpp:29:14: warning: unused variable 'psiDesired' [-Wunused-const-variable]
const double psiDesired = 0.;
             ^
2 warnings generated.
[ 66%] Building CXX object CMakeFiles/mpc.dir/src/main.cpp.o
/tmp/CarNDMPCProject/src/main.cpp:786:24: warning: unused variable 'sV' [-Wunused-variable]
          const double sV = solutionState[3];
                       ^
/tmp/CarNDMPCProject/src/main.cpp:738:24: warning: unused variable 'localX' [-Wunused-variable]
          const double localX=0.;
                       ^
/tmp/CarNDMPCProject/src/main.cpp:739:24: warning: unused variable 'localY' [-Wunused-variable]
          const double localY=0.;
                       ^
/tmp/CarNDMPCProject/src/main.cpp:785:24: warning: unused variable 'sPsi' [-Wunused-variable]
          const double sPsi = solutionState[2];
                       ^
/tmp/CarNDMPCProject/src/main.cpp:784:24: warning: unused variable 'sY' [-Wunused-variable]
          const double sY = solutionState[1];
                       ^
/tmp/CarNDMPCProject/src/main.cpp:787:24: warning: unused variable 'sCte' [-Wunused-variable]
          const double sCte = solutionState[4];
                       ^
/tmp/CarNDMPCProject/src/main.cpp:783:24: warning: unused variable 'sX' [-Wunused-variable]
          const double sX = solutionState[0];
                       ^
/tmp/CarNDMPCProject/src/main.cpp:788:24: warning: unused variable 'sEpsi' [-Wunused-variable]
          const double sEpsi = solutionState[5];
                       ^
8 warnings generated.
[100%] Linking CXX executable mpc
ld: warning: directory not found for option '-L/usr/local/Cellar/libuv/1*/lib'
[100%] Built target mpc
Softwares-MacBook-Pro:build david$ ls
CMakeCache.txt    CMakeFiles    Makefile    cmake_install.cmake mpc
Softwares-MacBook-Pro:build david$ ./mpc
Listening to port 4567
^C
Softwares-MacBook-Pro:build david$
```

#### Implementation
##### Polynomial Fitting and MPC Preprocessing

The state for each iteration is obtained from the sensor telemetry data transmitted by the simulator, e.g.

```
42["telemetry",{
"ptsx":[-93.05002,-107.7717,-123.3917,-134.97,-145.1165,-158.3417],
"ptsy":[65.34102,50.57938,33.37102,18.404,4.339378,-17.42898],
"psi_unity":3.957497,"psi":3.896485,
"x":-93.00126,"y":65.01852,
"steering_angle":-0.00553279,"throttle":0.1,"speed":10.52398}]
```

The sensor telemetry data is [pulled from the json packet](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L692-L697)

```  C++
vector<double> globalPtsX = j[1]["ptsx"];
vector<double> globalPtsY = j[1]["ptsy"];
double globalPx = j[1]["x"];
double globalPy = j[1]["y"];
double globalPsi = j[1]["psi"];
double v = j[1]["speed"];
```

The waypoints `globalPtsX` and the `globalPtsY` from the simulator sensor data were [transformed to car coordinates](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L742-L743)
``` C++
const vector<Eigen::Vector2d> xyCarSensorPath=transformMapToCar(globalPx, globalPy, -globalPsi, globalPtsX, globalPtsY);
```

by using an [affine transformation](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L270-L283) consisting of a rotation followed by a translation:

``` C++
const Eigen::Affine2d translateToOrigin(Eigen::Translation<double,2>(-theCarMapPositionX,-theCarMapPositionY));
const Eigen::Affine2d rotation((Eigen::Rotation2D<double>(theCarMapRotation)));
```
The transformed waypoints were then [fitted to a polynomial](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L742-L743)

``` 
const Eigen::VectorXd localCoeffs = polyfit(pullXCoordinate(xyCarSensorPath), pullYCoordinate(xyCarSensorPath), POLYFIT);// polynomial fit for midline of road
```

##### Model Predictive Control with Latency

The velocity from the car sensor data, `v`, is [collected in a fixed length buffer](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L414-L431)

```
addVelocityToBuffer(v);
const double a0=averageAcceleration();
```

and the velocity deltas between iterations are use to [get the average car acceleration](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L414-L431).
In the same way, the average amount of time for the last few iterations is [tracked](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L400-L412)
and used to determine how much time will pass before the car receives a control signal:

```
const double dt=averageLoopTime();
```

In my case and on my computer, this value was around .12 seconds: .1 seconds of imposed delay and an additional .02 seconds of calculation delay.

With the accleration, `a0`, the iteration time, `dt`, and the polynomial coefficients fitted to the waypoints, `localCoeffs`; the same state model used in the optimizer could be use to [calculate the initial state](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L738-L775) of the car at time `dt` in the future

``` C++
const double localPsi=0.;
const double localPsi0 = localPsi;
// fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
const double localX0 = 0;
const double localX1 = localX0+v0*cos(localPsi0)*dt;// car will continue in the local x direction, cos(0)=1
// fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
const double localY0 = 0;
const double localY1 = localY0+v0*sin(localPsi0)*dt;// car will continue in the local x direction, sin(0)=0
// fg[1 + psi_start + t] = psi1-(psi0-(v0/Lf)*delta0*dt);
const double Lf = 2.67;
const double delta0 = j[1]["steering_angle"];
const double localPsi1 = localPsi0-(v0/Lf)*delta0*dt;// car will continue in the localPsi=0;
// fg[1 + v_start + t] = v1-(v0+a0*dt);
const double v1 = v0+a0*dt;
// fg[1 + cte_start + t] = cte1-(cte0 + v0*CppAD::sin(epsi0)*dt);
const double localCte0 = polyeval(localCoeffs, localX0);
const double localEpsi0 = tangentialAngle(localCoeffs, localX0);
const double localCte1 = localCte0+v0*sin(localEpsi0)*dt;
// fg[1 + epsi_start + t] = epsi1-(psiError + (v0/Lf)*delta0*dt);
const double localEps1 = localEpsi0+(v0/Lf)*delta0*dt;
```

The state values at `dt` are stored in a state vector, `localState` and [passed to MPC](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L738-L775) to solve the nonlinear optimization problem:

``` C++
Eigen::VectorXd localState(6);
localState << localX1,localY1,localPsi1,v1,localCte1,localEps1;
const vector<vector<double>> mpcSolution = mpc.Solve(localState, localCoeffs);
```

##### The Model

The nonlinear optimization model has 2 parts: [the constraints](https://github.com/autohandle/CarNDMPCProject/blob/def705e76b0e3c370ba75dac1220d152f887586b/src/MPC.cpp#L144-L214) and [the costs](https://github.com/autohandle/CarNDMPCProject/blob/def705e76b0e3c370ba75dac1220d152f887586b/src/MPC.cpp#L88-L138) (to be minimized)

###### The Constraints

[The constraints](https://github.com/autohandle/CarNDMPCProject/blob/def705e76b0e3c370ba75dac1220d152f887586b/src/MPC.cpp#L144-L214) consist of 3 parts: the equations of motion,
``` C++
AD<double> x1 = vars[x_start + t];

AD<double> x0 = vars[x_start + t - 1];
AD<double> psi0 = vars[psi_start + t - 1];
AD<double> v0 = vars[v_start + t - 1];

fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);

// y(t+1)=y(t)+v(t)*sin(phi(t))*dt
const AD<double> y1 = vars[y_start + t];
const AD<double> y0 = vars[y_start + t - 1];
fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);

// ψ(t+1)=ψ(t)+v(t)/Lf*δ(t)*dt
const AD<double> psi1 = vars[psi_start + t];
//const AD<double> psi0 = vars[psi_start + t - 1];
const AD<double> delta0 = vars[delta_start + t -1];
//fg[1 + psi_start + t] = psi1-(psi0+v0/Lf*d0*dt);
fg[1 + psi_start + t] = psi1-(psi0-(v0/Lf)*delta0*dt);// project notes tips & tricks: δ is positive rotating counter clockwise

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
fg[1 + epsi_start + t] = epsi1-(psiError + (v0/Lf)*delta0*dt);
```

the [bounds on the state variables](https://github.com/autohandle/CarNDMPCProject/blob/def705e76b0e3c370ba75dac1220d152f887586b/src/MPC.cpp#L288-L302),
``` C++
// Acceleration/decceleration upper and lower limits.
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
```

and the [initial values](https://github.com/autohandle/CarNDMPCProject/blob/def705e76b0e3c370ba75dac1220d152f887586b/src/MPC.cpp#L304-L316) of the [state variables](https://github.com/autohandle/CarNDMPCProject/blob/4e954c12caf3172bb16529301f37b68dd96015e6/src/main.cpp#L738-L775) passed into [MPC::Solve](https://github.com/autohandle/CarNDMPCProject/blob/def705e76b0e3c370ba75dac1220d152f887586b/src/MPC.cpp#L225-L377)

``` C++
const double x = state[0];
const double y = state[1];
const double psi = state[2];
const double v = state[3];
const double cte = state[4];
const double epsi = state[5];

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
```

###### The Costs

[The Costs](https://github.com/autohandle/CarNDMPCProject/blob/def705e76b0e3c370ba75dac1220d152f887586b/src/MPC.cpp#L88-L138) to be minimized consist of:

the desired cross track error,

``` C++
const CppAD::AD<double> cteError = vars[cte_start+t]-cteDesired;
cteCost += CppAD::pow(cteError, 2); // minimize path errors in y (cte predicted)
```

the desired velocity

```C++
const CppAD::AD<double> vState = vars[v_start+t];
vCost += CppAD::pow(vState-velocityDesired, 2);
```

the desired orientation

```C++
const CppAD::AD<double> epsiState = vars[psi_start+t];
const CppAD::AD<double> epsiError=vars[epsi_start+t]-epsiDesired;
epsiCost += CppAD::pow(epsiError, 2);
```

the use of steering and acceleration

``` C++
controlSignalCost += CppAD::pow(vars[delta_start+t], 2); // minimize control changes
controlSignalCost += CppAD::pow(vars[a_start+t], 2); // minimize control changes
```

the size of the control change

``` C++
controlSignalDeltaCost += CppAD::pow(vars[delta_start+t+1] - vars[delta_start+t], 2);
controlSignalDeltaCost += CppAD::pow(vars[a_start+t+1] - vars[a_start+t], 2);
```

All of the cost components were weighted and combined

``` C++
const double CTECOST=500.;
const double VCOST=1.;
const double EPSICOST=7000.;
const double CONTROLSIGNALCOST=5.;
const double CONTROLSIGNALDELTACOST=500.;

fg[0] = (CTECOST*cteCost)+(EPSICOST*epsiCost)+(CONTROLSIGNALCOST*controlSignalCost)+(CONTROLSIGNALDELTACOST*controlSignalDeltaCost)+(VCOST*vCost);
```

The cost weights were determined by trial and error.

###### The Actuators

###### The Update Equations.

##### Timestep Length and Elapsed Duration (N & dt)


#### Simulation - successfully drive a lap around the track

The video of the car with the parameters in the checked-in code:
[third order polynomial](https://s3.amazonaws.com/autohandle.com/video/CarNDMPCProject_Polyfit3.mp4)
[second order polnomial](https://s3.amazonaws.com/autohandle.com/video/CarNDMPCProject_Polyfit2.mp4)

Too late in the project, I discovered that a negative throttle would apply the brakes, so I changed the clamp for the throttle to vary from -1 to +1
<br>
![-1 to +1 Logistic Function](./images/1Minus2xAbsLogFuncAt095.png)
<br>
Thr car is now a bit more aggresive and skids around the first turn after the bridge with the brakes on:
[Braking PID Controller](https://s3.amazonaws.com/autohandle.com/video/CarNDPIDControlProjectBraking.mp4)


The video was created by using a [screen recording tool](http://www.idownloadblog.com/2016/02/26/how-to-record-part-of-mac-screen-quicktime/).

