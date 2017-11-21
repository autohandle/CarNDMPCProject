#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include <Eigen/Dense>
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

const Eigen::Affine3d rotate(double ax, double ay, double az) {
  const Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(ax, Eigen::Vector3d(1, 0, 0)));
  const Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(ay, Eigen::Vector3d(0, 1, 0)));
  const Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
  return rz * ry * rx;
}

const Eigen::Rotation2D<double> rotateZ(double az) {
  //const Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
  Eigen::Rotation2D<double> rotate(az);
  return rotate;
}

const Eigen::Affine3d rotate(double az) {
  return rotate(0., 0., az);
}
/*
const Eigen::Affine3d mapToCarTransformation(const double theCarMapPositionX, const double theCarMapPositionY, const double theCarMapRotation) {
  cout << "theCarMapPositionX:" << theCarMapPositionX << ", theCarMapPositionY:" << theCarMapPositionY << ", theCarMapRotation:" << theCarMapRotation << std::endl;
  const Eigen::Affine3d translate(Eigen::Translation<double,3>(-theCarMapPositionX,-theCarMapPositionY, 1.));
  const Eigen::Affine3d rotation=rotate(theCarMapRotation);
  cout << "translate:" << std::endl << translate.matrix() << " ," << std::endl << "rotation:" << std::endl << rotation.matrix() << std::endl << " ," << std::endl << "rotate*translate:" << std::endl << (rotation*translate).matrix() << std::endl;
  throw 221;
  return rotation*translate;
}
*/
const Eigen::Affine2d mapToCarTransformation(const double theCarMapPositionX, const double theCarMapPositionY, const double theCarMapRotation) {
  cout << "theCarMapPositionX:" << theCarMapPositionX << ", theCarMapPositionY:" << theCarMapPositionY << ", theCarMapRotation:" << theCarMapRotation << std::endl;
  const Eigen::Affine2d translation(Eigen::Translation<double,2>(-theCarMapPositionX,-theCarMapPositionY));
  const Eigen::Affine2d  rotation((Eigen::Rotation2D<double>(theCarMapRotation)));
  cout << "translation:" << std::endl << translation.matrix() << " ," << std::endl << "rotation:" << std::endl << rotation.matrix() << std::endl << " ," << std::endl << "rotation*translation:" << std::endl << (rotation*translation).matrix() << std::endl;
  return rotation*translation;
}

const Eigen::Affine2d mapToCarTransformation(const Eigen::Vector2d theCarMapPosition, const double theCarMapRotation) {
  return mapToCarTransformation(theCarMapPosition[0], theCarMapPosition[1], theCarMapRotation);
}

/*
const Eigen::Affine3d mapToCarTransformation(const Eigen::Vector2d theCarMapPosition, const double theCarMapRotation) {
  return mapToCarTransformation(theCarMapPosition[0], theCarMapPosition[1], theCarMapRotation);
}
*/
const Eigen::Vector3d transformMapToCar(const Eigen::Affine2d theTransformationFromMapToCar, const Eigen::Vector3d theMapPosition) {
  return theTransformationFromMapToCar*theMapPosition;
}

const Eigen::Vector2d transformMapToCar(const Eigen::Affine2d theTransformationFromMapToCar, const Eigen::Vector2d theMapPosition) {
  const Eigen::Vector3d homogeneousPosition(theMapPosition[0], theMapPosition[1], 1.);
  const Eigen::Vector3d transformedPosition=transformMapToCar(theTransformationFromMapToCar,homogeneousPosition);
  const Eigen::Vector2d xyPosition(transformedPosition[0], transformedPosition[1]);
  return xyPosition;
}

const vector<Eigen::Vector2d> transformMapToCar(const Eigen::Affine2d theTransformationFromMapToCar, const vector<Eigen::Vector2d> theMapPositions) {
  vector<Eigen::Vector2d> transformedPositions;
  for (int p=0; p<theMapPositions.size(); p++) {
    transformedPositions.push_back(transformMapToCar(theTransformationFromMapToCar, theMapPositions[p]));
  }
  return transformedPositions;
}

bool areSame(const double a, const double b, const double epsilon) {
  return fabs(a - b) < epsilon;
}

bool isZero(const double a, const double epsilon) {
  return fabs(a) < epsilon;
}

bool isNotZero(const double a, const double epsilon) {
  return !isZero(a, epsilon);
}

bool areSame(const  Eigen::VectorXd a, const  Eigen::VectorXd b, const double epsilon) {
  for (int r=0; r<a.rows(); r++) {
    if (!areSame(a(r),b(r), epsilon)) {
      std::cout << std::endl
      << "a(" << r << "):" << a(r) << " != "
      << "b(" << r << "):" << b(r) << " == "
      << "fabs(a-b):" << fabs(a(r) - b(r))
      << std::endl;
      return false;
    }
  }
  return true;
}

bool areSame(const  Eigen::MatrixXd a, const  Eigen::MatrixXd b, const double epsilon) {
  for (int r=0; r<a.rows(); r++) {
    for (int c=0; c<a.cols(); c++) {
      if (!areSame(a(r,c),b(r,c), epsilon)) {
        std::cout << std::endl
        << "a(" << r << "," << c << "):" << a(r,c) << " != "
        << "b(" << r << "," << c << "):" << b(r,c) << " == "
        << "fabs(a-b):" << fabs(a(r,c) - b(r,c))
        << std::endl;
        return false;
      }
    }
  }
  return true;
}

static const double EPSILON=1.e-5;

bool areSame(const  Eigen::VectorXd a, const  Eigen::VectorXd b) {
  return areSame(a, b, EPSILON);
}

bool areSame( Eigen::MatrixXd a,  Eigen::MatrixXd b) {
  return areSame(a, b, EPSILON);
}

bool areSame(double a, double b) {
  return areSame(a, b, EPSILON);
}

bool isZero(double a) {
  return isZero(a, EPSILON);
  
}

bool isNotZero(double a) {
  return isNotZero(a, EPSILON);
}

const int testCarSameAsMap() {
  const Eigen::Vector2d carLocationOnMap(0.,0.);
  const double carRotation=0.;
  cout << "testCarCoordinatesSameAsMap-carLocationOnMap:" << carLocationOnMap.matrix() << std::endl << "mapToCarTransformation:" << mapToCarTransformation(carLocationOnMap, carRotation).matrix() << std::endl;
  const Eigen::Vector2d mapLocation0(0.,0.);
  const Eigen::Vector2d mapLocationFromCar0=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation0);
  assert(areSame(mapLocationFromCar0[0],0.) && areSame(mapLocationFromCar0[1],0.));
  const Eigen::Vector2d mapLocation1(1.,1.);
  const Eigen::Vector2d mapLocationFromCar1=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation1);
  cout << "testCarCoordinatesSameAsMap-mapLocationFromCar1:" << mapLocationFromCar1.matrix() << std::endl;
  assert(areSame(mapLocationFromCar1[0],1.) && areSame(mapLocationFromCar1[1],1.));
  const Eigen::Vector2d mapLocation2(-1.,1.);
  const Eigen::Vector2d mapLocationFromCar2=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation2);
  assert(areSame(mapLocationFromCar2[0],-1.) && areSame(mapLocationFromCar2[1],1.));
  const Eigen::Vector2d mapLocation3(-1.,-1.);
  const Eigen::Vector2d mapLocationFromCar3=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation3);
  assert(areSame(mapLocationFromCar3[0],-1.) && areSame(mapLocationFromCar3[1],-1.));
  const Eigen::Vector2d mapLocation4(1.,-1.);
  const Eigen::Vector2d mapLocationFromCar4=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation4);
  assert(areSame(mapLocationFromCar4[0],1.) && areSame(mapLocationFromCar4[1],-1.));
  return 0;
}

const int rotateCarAtOriginByPiOver2() {
  const Eigen::Vector2d carLocationOnMap(0.,0.);
  const double carRotation=-pi()/2.;
  cout << "rotateCarCoordinatesAtOriginByPiOver2-carLocationOnMap:" << carLocationOnMap.matrix() << std::endl << "mapToCarTransformation:" << mapToCarTransformation(carLocationOnMap, carRotation).matrix() << std::endl;
  const Eigen::Vector2d mapLocation0(0.,0.);// should stay same
  const Eigen::Vector2d mapLocationFromCar0=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation0);
  assert(areSame(mapLocationFromCar0[0],0.) && areSame(mapLocationFromCar0[1],0.));
  const Eigen::Vector2d mapLocation1(1.,1.);
  const Eigen::Vector2d mapLocationFromCar1=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation1);
  cout << "rotateCarCoordinatesAtOriginByPiOver2-mapLocationFromCar1:" << mapLocationFromCar1.matrix() << std::endl;
  assert(areSame(mapLocationFromCar1[0],1.) && areSame(mapLocationFromCar1[1],-1.));
  const Eigen::Vector2d mapLocation2(-1.,1.);
  const Eigen::Vector2d mapLocationFromCar2=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation2);
  assert(areSame(mapLocationFromCar2[0],1.) && areSame(mapLocationFromCar2[1],1.));
  const Eigen::Vector2d mapLocation3(-1.,-1.);
  const Eigen::Vector2d mapLocationFromCar3=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation3);
  assert(areSame(mapLocationFromCar3[0],-1.) && (mapLocationFromCar3[1],1.));
  const Eigen::Vector2d mapLocation4(1.,-1.);
  const Eigen::Vector2d mapLocationFromCar4=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation4);
  assert(areSame(mapLocationFromCar4[0],-1.) && areSame(mapLocationFromCar4[1],-1.));
  return 0;
}

const int rotateCarAtOriginByPi() {
  const Eigen::Vector2d carLocationOnMap(0.,0.);
  const double carRotation=-pi();
  cout << "rotateCarCoordinatesAtOriginByPiOver2-carLocationOnMap:" << carLocationOnMap.matrix() << std::endl << "transformMapToCar:" << mapToCarTransformation(carLocationOnMap, carRotation).matrix() << std::endl;
  const Eigen::Vector2d mapLocation0(0.,0.);// should stay same
  const Eigen::Vector2d mapLocationFromCar0=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation0);
  assert(areSame(mapLocationFromCar0[0],0.) && areSame(mapLocationFromCar0[1],0.));
  const Eigen::Vector2d mapLocation1(1.,1.);
  const Eigen::Vector2d mapLocationFromCar1=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation1);
  cout << "rotateCarCoordinatesAtOriginByPiOver2-mapLocationFromCar1:" << mapLocationFromCar1.matrix() << std::endl;
  assert(areSame(mapLocationFromCar1[0],-1.) && areSame(mapLocationFromCar1[1],-1.));
  const Eigen::Vector2d mapLocation2(-1.,1.);
  const Eigen::Vector2d mapLocationFromCar2=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation2);
  assert(areSame(mapLocationFromCar2[0],1.) && areSame(mapLocationFromCar2[1],-1.));
  const Eigen::Vector2d mapLocation3(-1.,-1.);
  const Eigen::Vector2d mapLocationFromCar3=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation3);
  assert(areSame(mapLocationFromCar3[0],1.) && (mapLocationFromCar3[1],1.));
  const Eigen::Vector2d mapLocation4(1.,-1.);
  const Eigen::Vector2d mapLocationFromCar4=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation4);
  assert(areSame(mapLocationFromCar4[0],-1.) && areSame(mapLocationFromCar4[1],1.));
  return 0;
}

const int translateCarToOneOne() {
  const Eigen::Vector2d carLocationOnMap(1.,1.);
  const double carRotation=0.;
  cout << "translateCarCoordinatesToOneOne-carLocationOnMap:" << carLocationOnMap.matrix() << std::endl << "transformMapToCar:" << mapToCarTransformation(carLocationOnMap, carRotation).matrix() << std::endl;
  const Eigen::Vector2d mapLocation0(0.,0.);// should stay same
  const Eigen::Vector2d mapLocationFromCar0=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation0);
  cout << "translateCarCoordinatesToOneOne-mapLocationFromCar0:" << mapLocationFromCar0.matrix() << std::endl;
  assert(areSame(mapLocationFromCar0[0],-1.) && areSame(mapLocationFromCar0[1],-1.));
  const Eigen::Vector2d mapLocation1(1.,1.);
  const Eigen::Vector2d mapLocationFromCar1=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation1);
  cout << "translateCarCoordinatesToOneOne-mapLocationFromCar1:" << mapLocationFromCar1.matrix() << std::endl;
  assert(areSame(mapLocationFromCar1[0],0.) && areSame(mapLocationFromCar1[1],0.));
  const Eigen::Vector2d mapLocation2(-1.,1.);
  const Eigen::Vector2d mapLocationFromCar2=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation2);
  assert(areSame(mapLocationFromCar2[0],-2.) && areSame(mapLocationFromCar2[1],0.));
  const Eigen::Vector2d mapLocation3(-1.,-1.);
  const Eigen::Vector2d mapLocationFromCar3=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation3);
  assert(areSame(mapLocationFromCar3[0],-2.) && (mapLocationFromCar3[1],-2.));
  const Eigen::Vector2d mapLocation4(1.,-1.);
  const Eigen::Vector2d mapLocationFromCar4=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation4);
  assert(areSame(mapLocationFromCar4[0],0.) && areSame(mapLocationFromCar4[1],-2.));
  return 0;
}

const int translateCarToOneOneRotatePi() {
  const Eigen::Vector2d carLocationOnMap(1.,1.);
  const double carRotation=pi();
  cout << "translateCarToOneOneRotatePi-carLocationOnMap:" << carLocationOnMap.matrix() << std::endl << "transformMapToCar:" << mapToCarTransformation(carLocationOnMap, carRotation).matrix() << std::endl;
  const Eigen::Vector2d mapLocation0(0.,0.);
  const Eigen::Vector2d mapLocationFromCar0=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation0);
  cout << "translateCarToOneOneRotatePi-mapLocationFromCar0:" << mapLocationFromCar0.matrix() << std::endl;
  assert(areSame(mapLocationFromCar0[0],1.) && areSame(mapLocationFromCar0[1],1.));
  const Eigen::Vector2d mapLocation1(1.,1.);
  const Eigen::Vector2d mapLocationFromCar1=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation1);
  cout << "translateCarToOneOneRotatePi-mapLocationFromCar1:" << mapLocationFromCar1.matrix() << std::endl;
  assert(areSame(mapLocationFromCar1[0],0.) && areSame(mapLocationFromCar1[1],0.));
  const Eigen::Vector2d mapLocation2(-1.,1.);
  const Eigen::Vector2d mapLocationFromCar2=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation2);
  assert(areSame(mapLocationFromCar2[0],2.) && areSame(mapLocationFromCar2[1],0.));
  const Eigen::Vector2d mapLocation3(-1.,-1.);
  const Eigen::Vector2d mapLocationFromCar3=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation3);
  assert(areSame(mapLocationFromCar3[0],2.) && (mapLocationFromCar3[1],2.));
  const Eigen::Vector2d mapLocation4(1.,-1.);
  const Eigen::Vector2d mapLocationFromCar4=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation4);
  assert(areSame(mapLocationFromCar4[0],0.) && areSame(mapLocationFromCar4[1],2.));
  return 0;
}

const int testMappingToCar() {
  testCarSameAsMap();// simplest car coordinates same as map
  rotateCarAtOriginByPiOver2();// rotate to the north
  rotateCarAtOriginByPi();// rotate to the west
  translateCarToOneOne();
  translateCarToOneOneRotatePi();
  return 0;
}

const bool RUNTESTS = true;

int main() {
  
  if (RUNTESTS) {
    testMappingToCar();
    return 0;
  }
  std::cout << "mapToCarTransform:" << mapToCarTransformation(-40., 108., 3.733).matrix() << std::endl;
  
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value;
          double throttle_value;
          
          throttle_value=0.1;

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
          
          mpc_x_vals.push_back(10.);mpc_x_vals.push_back(20.);mpc_x_vals.push_back(30.);mpc_x_vals.push_back(30.);
          mpc_y_vals.push_back(0.1);mpc_y_vals.push_back(0.2);mpc_y_vals.push_back(0.3);mpc_y_vals.push_back(0.4);
          
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          next_x_vals.push_back(0.);next_x_vals.push_back(10.);next_x_vals.push_back(20.);next_x_vals.push_back(30.);
          next_y_vals.push_back(0.1);next_y_vals.push_back(0.4);next_y_vals.push_back(0.9);next_y_vals.push_back(0.16);

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}


