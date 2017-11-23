#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <string>
#include <algorithm>    // std::min

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include <Eigen/Dense>
#include "MPC.h"
#include "json.hpp"
#include "assert.h"


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

vector<double> polyeval(Eigen::VectorXd theCoefficients, vector<double> theXValues) {
  vector<double> yValues;
  for (int x=0; x<theXValues.size(); x++) {
    double yValue=polyeval(theCoefficients, theXValues[x]);
    yValues.push_back(yValue);
  }
  return yValues;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
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

const string toString(vector<double> theVector) {
  std::ostringstream oss;
  
  if (!theVector.empty())
  {
    // Convert all but the last element to avoid a trailing ","
    std::copy(theVector.begin(), theVector.end()-1,
              std::ostream_iterator<int>(oss, ", "));
    
    // Now add the last element with no delimiter
    oss << theVector.back();
  }
  return oss.str();
}

// Eigen::Vector3d v2(v1.data())
// VectorXcd v3 = VectorXcd::Map(v2.data(), v2.size());
Eigen::VectorXd polyfit(const vector<double> theXValues, const vector<double> theYValues, const int thePolynomialOrder) {
  Eigen::VectorXd eigenX = Eigen::VectorXd::Map(theXValues.data(), theXValues.size());
  Eigen::VectorXd eigenY = Eigen::VectorXd::Map(theYValues.data(), theYValues.size());
  cout  << "eigenX:" << std::endl << eigenX.matrix() << ", " << std::endl
        << "eigenY:" << std::endl << eigenY.matrix() << ", " << std::endl
        << "thePolynomialOrder:" << thePolynomialOrder << std::endl
        << "polyfit:" << std::endl << polyfit(eigenX, eigenY, thePolynomialOrder) << std::endl;
  return polyfit(eigenX, eigenY, thePolynomialOrder);
}

const vector<double> selectValuesFollowingGradient(const vector<double> theValues) {
  const double EPSILON = 0.5;
  enum SIGNOFGRADIENT{POSITIVE,NEGATIVE};
  const SIGNOFGRADIENT theGradient = (theValues[1]-theValues[0])>0.?POSITIVE:NEGATIVE;
  vector<double> valuesFollowingTheGradient={theValues[0]};
  if (theGradient==POSITIVE) {// POSITIVE gradient
    for (int v=1; v<theValues.size(); v++) {
      if ((theValues[v]>theValues[v-1]) || abs(theValues[v]-theValues[v-1])>EPSILON) {// still increasing or close by
        valuesFollowingTheGradient.push_back(theValues[v]);
      } else {
        break;
      }
    }
  } else {// NEGATIVE gradient
    for (int v=1; v<theValues.size(); v++) {
      if ((theValues[v]<theValues[v-1]) || abs(theValues[v]-theValues[v-1])>EPSILON) {// still decreasing or close by
        valuesFollowingTheGradient.push_back(theValues[v]);
      } else {
        break;
      }
    }
  }
  return valuesFollowingTheGradient;
}

const vector<double> keepValues(const vector<double> theValues, const int theNumberOfValues) {
  vector<double> values;
  for (int v=0; v<theNumberOfValues; v++) {
      values.push_back(theValues[v]);
  }
  return values;
}

Eigen::VectorXd trimmedPolyfit(const vector<double> theXValues, const vector<double> theYValues, const int thePolynomialOrder) {
  const vector<double> xTrackingGradient=selectValuesFollowingGradient(theXValues);
  const vector<double> matchingY=keepValues(theYValues, xTrackingGradient.size());
  cout  << "xTrackingGradient:" << std::endl << toString(xTrackingGradient) << ", " << std::endl
  << "matchingY:" << std::endl << toString(matchingY) << ", " << std::endl;
  Eigen::VectorXd coefficients=polyfit(xTrackingGradient, matchingY, thePolynomialOrder<(xTrackingGradient.size()-1)?thePolynomialOrder:(xTrackingGradient.size()-1));
  coefficients.resize(thePolynomialOrder+1, 1);
  for (int c=xTrackingGradient.size(); c<coefficients.size(); c++) {
    coefficients[c]=0.;
  }
  return coefficients;
}

const double fPrimeOfX(const Eigen::VectorXd thePolynomialCoefficients, const double theValueX) {
  //const int numberOfCoefficients=sizeof(thePolynomialCoefficients)/sizeof(thePolynomialCoefficients[0]);
  const int numberOfCoefficients = thePolynomialCoefficients.size();
  if (numberOfCoefficients > 0) {
    double valueAtX=thePolynomialCoefficients[1];
    for (int c=2;c<numberOfCoefficients;c++) {
      valueAtX+=c*thePolynomialCoefficients[c]*pow(theValueX,c-1);
    }
    cout << "fPrimeOfX-numberOfCoefficients:" << numberOfCoefficients << ", valueAtX:" << valueAtX << std::endl;
    return valueAtX;
  } else {
    return 0.;
  }
}

const double tangentialAngle(const Eigen::VectorXd thePolynomialCoefficients, const double theValueX) {
  return atan(fPrimeOfX(thePolynomialCoefficients, theValueX));
}

const Eigen::Affine2d affineMapToCarTransformation(const double theCarMapPositionX, const double theCarMapPositionY, const double theCarMapRotation) {
  cout << "theCarMapPositionX:" << theCarMapPositionX << ", theCarMapPositionY:" << theCarMapPositionY << ", theCarMapRotation:" << theCarMapRotation << std::endl;
  const Eigen::Affine2d translateToOrigin(Eigen::Translation<double,2>(-theCarMapPositionX,-theCarMapPositionY));
  const Eigen::Affine2d rotation((Eigen::Rotation2D<double>(theCarMapRotation)));
  //const Eigen::Affine2d translateToMap(Eigen::Translation<double,2>(theCarMapPositionX,theCarMapPositionY));
  cout  << "translateToOrigin:" << std::endl << translateToOrigin.matrix() << " ," << std::endl
        << "rotation:" << std::endl << rotation.matrix() << std::endl << " ," << std::endl
        //<< "translateToMap:" << std::endl << translateToMap.matrix() << " ," << std::endl
        << "rotation*translateToOrigin:" << std::endl << (rotation*translateToOrigin).matrix()
        << std::endl;
  return rotation*translateToOrigin;
}

const Eigen::Matrix3d matrixMapToCarTransformation(const double theMapPositionX, const double theMapPositionY, const double theMapRotation) {
  cout << "theMapPositionX:" << theMapPositionX << ", theMapPositionY:" << theMapPositionY << ", theMapRotation:" << theMapRotation << std::endl;
  const double sinMapRotation=sin(theMapRotation);
  const double cosMapRotation=cos(theMapRotation);
  Eigen::Matrix3d rotationAndTranslation;
  rotationAndTranslation <<
    cosMapRotation, -sinMapRotation,  -theMapPositionX,
    sinMapRotation, cosMapRotation,   -theMapPositionY,
    0.,             0.,               1;
  cout << "rotationAndTranslation:" << std::endl << rotationAndTranslation.matrix() << std::endl;
  return rotationAndTranslation;
}

const Eigen::Affine2d mapToCarTransformation(const Eigen::Vector2d theCarMapPosition, const double theCarMapRotation) {
  return affineMapToCarTransformation(theCarMapPosition[0], theCarMapPosition[1], theCarMapRotation);
}

const Eigen::Matrix3d mapToCarTransformationX(const Eigen::Vector2d theCarMapPosition, const double theCarMapRotation) {
  return matrixMapToCarTransformation(theCarMapPosition[0], theCarMapPosition[1], theCarMapRotation);
}

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

const Eigen::Vector3d transformMapToCar(const Eigen::Matrix3d theTransformationFromMapToCar, const Eigen::Vector3d theMapPosition) {
  return theTransformationFromMapToCar*theMapPosition;
}

const Eigen::Vector2d transformMapToCar(const Eigen::Matrix3d theTransformationFromMapToCar, const Eigen::Vector2d theMapPosition) {
  const Eigen::Vector3d homogeneousPosition(theMapPosition[0], theMapPosition[1], 1.);
  const Eigen::Vector3d transformedPosition=transformMapToCar(theTransformationFromMapToCar,homogeneousPosition);
  const Eigen::Vector2d xyPosition(transformedPosition[0], transformedPosition[1]);
  return xyPosition;
}

const vector<Eigen::Vector2d> transformMapToCar(const Eigen::Matrix3d theTransformationFromMapToCar, const vector<Eigen::Vector2d> theMapPositions) {
  vector<Eigen::Vector2d> transformedPositions;
  for (int p=0; p<theMapPositions.size(); p++) {
    transformedPositions.push_back(transformMapToCar(theTransformationFromMapToCar, theMapPositions[p]));
  }
  return transformedPositions;
}

const vector<Eigen::Vector2d> transformMapToCar(const double theCarXPosition, const double theCarYPosition, const double theCarRotation, const vector<double> theMapXPositions, const vector<double> theMapYPositions) {
  assert(theMapXPositions.size() == theMapYPositions.size());
  vector<Eigen::Vector2d> mapPositions;
  for (int p=0; p<theMapXPositions.size(); p++) {
    mapPositions.push_back(Eigen::Vector2d(theMapXPositions[p], theMapYPositions[p]));
  }
  //const Eigen::Matrix3d transformation=matrixMapToCarTransformation(theCarXPosition, theCarYPosition, theCarRotation);
  const Eigen::Affine2d transformation=affineMapToCarTransformation(theCarXPosition, theCarYPosition, theCarRotation);
  return transformMapToCar(transformation, mapPositions);
}

const vector<double> pullCoordinate(vector<Eigen::Vector2d> thePoints, const int theCoordinate) {
  vector<double> coordinates;
  for (int p=0; p<thePoints.size(); p++) {
    coordinates.push_back(thePoints[p][theCoordinate]);
  }
  return coordinates;
}

const vector<double> pullXCoordinate(vector<Eigen::Vector2d> thePoints) {
  return pullCoordinate(thePoints,0);
}

const vector<double> pullYCoordinate(vector<Eigen::Vector2d> thePoints) {
  return pullCoordinate(thePoints,1);
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

const int translateCarToTwoOneRotatePi() {
  const Eigen::Vector2d carLocationOnMap(2.,1.);
  const double carRotation=pi();
  cout << "translateCarToTwoOneRotatePi-carLocationOnMap:" << carLocationOnMap.matrix() << std::endl << "transformMapToCar:" << mapToCarTransformation(carLocationOnMap, carRotation).matrix() << std::endl;
  const Eigen::Vector2d mapLocation0(0.,0.);
  const Eigen::Vector2d mapLocationFromCar0=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation0);
  cout << "translateCarToTwoOneRotatePi-mapLocationFromCar0:" << mapLocationFromCar0.matrix() << std::endl;
  assert(areSame(mapLocationFromCar0[0],2.) && areSame(mapLocationFromCar0[1],1.));
  const Eigen::Vector2d mapLocation1(1.,1.);
  const Eigen::Vector2d mapLocationFromCar1=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation1);
  cout << "translateCarToTwoOneRotatePi-mapLocationFromCar1:" << mapLocationFromCar1.matrix() << std::endl;
  assert(areSame(mapLocationFromCar1[0],1.) && areSame(mapLocationFromCar1[1],0.));
  const Eigen::Vector2d mapLocation2(-1.,1.);
  const Eigen::Vector2d mapLocationFromCar2=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation2);
  assert(areSame(mapLocationFromCar2[0],3.) && areSame(mapLocationFromCar2[1],0.));
  const Eigen::Vector2d mapLocation3(-1.,-1.);
  const Eigen::Vector2d mapLocationFromCar3=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation3);
  assert(areSame(mapLocationFromCar3[0],3.) && (mapLocationFromCar3[1],2.));
  const Eigen::Vector2d mapLocation4(1.,-1.);
  const Eigen::Vector2d mapLocationFromCar4=transformMapToCar(mapToCarTransformation(carLocationOnMap, carRotation), mapLocation4);
  assert(areSame(mapLocationFromCar4[0],1.) && areSame(mapLocationFromCar4[1],2.));
  return 0;
}

const int nonZeroCoefficients(const Eigen::VectorXd theCoefficients) {
  int nonZeros=0;
  for (int c=0; c<theCoefficients.size(); c++) {
    if (theCoefficients[c]!=0.) {
      nonZeros++;
    }
  }
  return nonZeros;
}

const int runPolyfit() {
  
  const int polynomialOrder = 3;
  
  const vector<double> globalPtsX_1 = {-93.05002,-107.7717,-123.3917,-134.97,-145.1165,-158.3417};
  const vector<double> globalPtsY_1 = {65.34102,50.57938,33.37102,18.404,4.339378,-17.42898};
  const Eigen::VectorXd globalCoeffs1 = trimmedPolyfit(globalPtsX_1, globalPtsY_1, polynomialOrder);// polynomial fit for midline of road
  cout << "globalCoeffs:" << globalCoeffs1.size() << ":" << std::endl << globalCoeffs1 << std::endl;
  cout << std::endl;
  assert(globalCoeffs1.size()==(polynomialOrder+1));// 1 more than the polynomial order of 3
  assert(nonZeroCoefficients(globalCoeffs1)==polynomialOrder+1);

  const vector<double> globalPtsX_2 = {-107.7717,-123.3917,-134.97,-145.1165,-158.3417,-164.3164};
  const vector<double> globalPtsY_2 = {50.57938,33.37102,18.404,4.339378,-17.42898,-30.18062};
  const Eigen::VectorXd globalCoeffs2 = trimmedPolyfit(globalPtsX_2, globalPtsY_2, polynomialOrder);// polynomial fit for midline of road
  cout << "globalCoeffs:" << globalCoeffs2.size() << ":" << std::endl << globalCoeffs2 << std::endl;
  cout << std::endl;
  assert(globalCoeffs2.size()==(polynomialOrder+1));// 1 more than the polynomial order of 3
  assert(nonZeroCoefficients(globalCoeffs2)==polynomialOrder+1);

  const vector<double> globalPtsX_3 = {-169.3365,-175.4917,-176.9617,-176.8864,-175.0817,-170.3617};
  const vector<double> globalPtsY_3 = {-42.84062,-66.52898,-76.85062,-90.64063,-100.3206,-115.129};
  const Eigen::VectorXd globalCoeffs3 = trimmedPolyfit(globalPtsX_3, globalPtsY_3, polynomialOrder);// polynomial fit for midline of road
  cout << "globalCoeffs:" << globalCoeffs3.size() << ":" << std::endl << globalCoeffs3 << std::endl;
  cout << std::endl;
  assert(globalCoeffs3.size()==(polynomialOrder+1));// 1 more than the polynomial order of 3
  assert(nonZeroCoefficients(globalCoeffs3)==3);// polnomial order is 2 -> 3 coefficients

  const vector<double> globalPtsX_4 = {-176.9617,-176.8864,-175.0817,-170.3617,-164.4217,-158.9417};
  const vector<double> globalPtsY_4 = {-76.85062,-90.64063,-100.3206,-115.129,-124.5206,-131.399};
  const Eigen::VectorXd globalCoeffs4 = trimmedPolyfit(globalPtsX_4, globalPtsY_4, polynomialOrder);// polynomial fit for midline of road
  cout << "globalCoeffs:" << globalCoeffs4.size() << ":" << std::endl << globalCoeffs4 << std::endl;
  cout << std::endl;
  assert(globalCoeffs4.size()==(polynomialOrder+1));// 1 more than the polynomial order of 3
  assert(nonZeroCoefficients(globalCoeffs4)==polynomialOrder+1);
  
  const vector<double> globalPtsX_5 = {129.1083,129.1283,126.8983,122.3383,117.2083,98.34827};
  const vector<double> globalPtsY_5 = {-108.669,-100.349,-89.95898,-79.97897,-69.827,-42.02898};
  const Eigen::VectorXd globalCoeffs5 = trimmedPolyfit(globalPtsX_5, globalPtsY_5, polynomialOrder);// polynomial fit for midline of road
  cout << "globalCoeffs:" << globalCoeffs5.size() << ":" << std::endl << globalCoeffs5 << std::endl;
  cout << std::endl;
  assert(globalCoeffs5.size()==(polynomialOrder+1));// 1 more than the polynomial order of 3
  assert(nonZeroCoefficients(globalCoeffs5)==polynomialOrder+1);
  
  const vector<double> globalPtsX_6 = {129.1083,129.1283,126.8983,122.3383,117.2083,98.34827};
  const vector<double> globalPtsY_6 = {-108.669,-100.349,-89.95898,-79.97897,-69.827,-42.02898};
  const Eigen::VectorXd globalCoeffs6 = trimmedPolyfit(globalPtsX_6, globalPtsY_6, polynomialOrder);// polynomial fit for midline of road
  cout << "globalCoeffs:" << globalCoeffs6.size() << ":" << std::endl << globalCoeffs6 << std::endl;
  cout << std::endl;
  assert(globalCoeffs6.size()==(polynomialOrder+1));// 1 more than the polynomial order of 3
  assert(nonZeroCoefficients(globalCoeffs6)==polynomialOrder+1);

  return 0.;
}

const int testMappingToCar() {
  testCarSameAsMap();// simplest car coordinates same as map
  rotateCarAtOriginByPiOver2();// rotate to the north
  rotateCarAtOriginByPi();// rotate to the west
  translateCarToOneOne();
  translateCarToOneOneRotatePi();
  translateCarToTwoOneRotatePi();
  runPolyfit();
  return 0;
}

const bool RUNTESTS = true;

int main() {
  
  if (RUNTESTS) {
    testMappingToCar();
    return 0;
  }
  std::cout << "mapToCarTransform:" << matrixMapToCarTransformation(-40., 108., 3.733).matrix() << std::endl;
  
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
          vector<double> globalPtsX = j[1]["ptsx"];
          vector<double> globalPtsY = j[1]["ptsy"];
          double globalPx = j[1]["x"];
          double globalPy = j[1]["y"];
          double globalPsi = j[1]["psi"];
          double v = j[1]["speed"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          const Eigen::VectorXd globalCoeffs = trimmedPolyfit(globalPtsX, globalPtsY, 3);// polynomial fit for midline of road
          cout << "globalCoeffs:" << globalCoeffs.size() << ":" << globalCoeffs << std::endl;

          const double globalYDesired=polyeval(globalCoeffs, globalPx);
          const double cte = globalYDesired-globalPy;// difference between current y and the path y at the current x
          cout << "globalYDesired:" << globalYDesired << ", py:" << globalPy << ", cte:" << cte << std::endl;
          // TODO: calculate the orientation error
          //double epsi = ? ;
          // eψ(t+1)  =  ψ(t)−ψdes(t)+((v(t)/Lf)*δ(t)*dt)
          //    ψdes(t) is desired orientation angle psi/ψ,
          //      the desired ψ is the tangent of the path f at t (i.e. f'(t) for path f(t)),
          //      also called: tangential angle, it can be calcualted as: arctan(f'(x(t))).
          //      the path f(t) is calculated by polynomial fitting a series of point x,y to get coefficients: c,
          //      so f'(t) can be calulated from c's, by hand, depending on the degree of the fitted polynomial
          //
          //      there should be only 2 coefficients for a 1 degree polynomial: c0 & c1
          //      therefore f'(t) should be just c1.
          double psiDesiredCalculated = atan(globalCoeffs[1]+2.*globalCoeffs[2]*globalPx+3.*globalCoeffs[3]*pow(globalPx,2));
          const double globalPsiDesired = tangentialAngle(globalCoeffs, globalPx);
          cout << "globalPsiDesired:" << globalPsiDesired << ", psiDesiredCalculated:" << psiDesiredCalculated << std::endl;
          assert(areSame(psiDesiredCalculated, globalPsiDesired));
          const double epsi = globalPsi-globalPsiDesired;
          cout << "epsi:" << epsi << std::endl;
          
          Eigen::VectorXd state(6);
          state << globalPx,globalPy,globalPsi,v,cte,epsi;
          const vector<double> globalSolution = mpc.Solve(state, globalCoeffs);

          const double sX = globalSolution[0];
          const double sY = globalSolution[1];
          const double sPsi = globalSolution[2];
          const double sV = globalSolution[3];
          const double sCte = globalSolution[4];
          const double sEpsi = globalSolution[5];
          double sSteerValue = globalSolution[6];
          double sThrottleValue = globalSolution[7];
          cout << "sSteerValue:" << sSteerValue << ", sThrottleValue:" << sThrottleValue << std::endl;

          const double steer_value=sSteerValue;
          const double throttle_value=0.1;

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          cout << "globalPsi:" << globalPsi << " = " << rad2deg(globalPsi) << " degrees, sine:" << sin(globalPsi) << ", cosine:" << cos(globalPsi) << std::endl;
          const vector<Eigen::Vector2d> nextXYCarValues=transformMapToCar(globalPx, globalPy, -globalPsi, globalPtsX, globalPtsY);

          //Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
          
          mpc_x_vals.push_back(10.);mpc_x_vals.push_back(20.);mpc_x_vals.push_back(30.);mpc_x_vals.push_back(30.);
          mpc_y_vals.push_back(0.1);mpc_y_vals.push_back(0.2);mpc_y_vals.push_back(0.3);mpc_y_vals.push_back(0.4);
          //mpc_x_vals=pullXCoordinate(nextXYCarValues);
          //mpc_y_vals=pullYCoordinate(nextXYCarValues);
          
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          next_x_vals=pullXCoordinate(nextXYCarValues);
          next_y_vals=pullYCoordinate(nextXYCarValues);
          cout << "next_x_vals:" << toString(next_x_vals) << std::endl;
          cout << "next_y_vals:" << toString(next_y_vals) << std::endl;
          //throw 122;

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
          //this_thread::sleep_for(chrono::milliseconds(100));
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


