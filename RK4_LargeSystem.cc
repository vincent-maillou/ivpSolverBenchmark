/**
 * @file RK4_LargeSystem.cc
 * @author Vincent MAILLOU @ AMOLF
 * @brief 
 * @version 0.1
 * @date 2022-09-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <iostream>
#include <stdio.h>

#include <math.h> 
#include <random>
#include <vector>
#include <fstream>
#include <chrono>

#include <Eigen/Dense>

#define PI 3.14159265359
using reel = double;

typedef Eigen::Matrix<reel, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::DiagonalMatrix< reel, Eigen::Dynamic> DiagonalMatrixXd;
typedef Eigen::Matrix<reel, 4, 4> Matrix4d;
typedef Eigen::Matrix<reel, Eigen::Dynamic, 1> VectorXd;



/**
 * @brief M, B and K Matrix defining the following motion equation:
 * M.X_dotdot + B.X_dot + K.X = F(t)
 * 
 */
struct MBKMatrix{
  MatrixXd M;
  MatrixXd B;
  MatrixXd K;

  std::ostream& affiche(std::ostream& out) const;
};

std::ostream& MBKMatrix::affiche(std::ostream& out) const{
  std::cout << "M Matrix:" << std::endl;
  std::cout << M << std::endl << std::endl;

  std::cout << "B Matrix:" << std::endl << std::endl;
  std::cout << B << std::endl;

  std::cout << "K Matrix:" << std::endl << std::endl;
  std::cout << K << std::endl;

  return out;
}

std::ostream& operator<<(std::ostream& sortie, MBKMatrix const& MBK){
    return MBK.affiche(sortie);
}



/**
 * @brief Contain F vector, evolution of the excitation over time and S "Shape" matrix
 * wich define the part of the lattice were the force applie
 * 
 */
struct FShapeMatrix{
  VectorXd F; 
  VectorXd S;

  std::ostream& affiche(std::ostream& out) const;
};

std::ostream& FShapeMatrix::affiche(std::ostream& out) const{
  std::cout << "S vector:" << std::endl;
  std::cout << S.transpose() << std::endl;

  return out;
}

std::ostream& operator<<(std::ostream& sortie, FShapeMatrix const& FS){
    return FS.affiche(sortie);
}



/**
 * @brief Define a square lattice of masses connected with there 
 * closer neigbourgs by a spring. Here we only considere first order/range
 * interaction between the elements of the lattice.
 */
class Lattice{
 private:
  uint excitationFrequency; // Number of initial excitation points, DEFAULT: 72kHz
  uint sampleRate; // [Hz], Number of F(t) points generated per second, DEFAULT: 2*4*72kHz = 576000
  uint simDuration; // TODO: Implement the integer variation of simulation
  uint latticeSize; // Size of the lattice across 1 dimmension
  uint numberOfElements; // Total number of elements in the lattice

  MBKMatrix MBK;
  FShapeMatrix FS;

  void InitialiseFSvector();
  void InitialiseMBKMatrix();

 public:
  Lattice(uint latticeSize_, uint excitationFrequency_ = 72000, uint pointsPerPeriode_ = 4, uint simDuration_ = 1) : 
    latticeSize(latticeSize_),
    excitationFrequency(excitationFrequency_),
    sampleRate(2 * pointsPerPeriode_ * excitationFrequency_ /* [Hz], F(t) upSampled 2 times for RK4 need*/),
    simDuration(simDuration_) 
    {
      numberOfElements = latticeSize*latticeSize;

      InitialiseFSvector();
      InitialiseMBKMatrix();
    }

  std::ostream& affiche(std::ostream& out) const;
  void toFile() const;

  uint getSize() const { return latticeSize; }
  uint getNDOF() const { return numberOfElements; }
  uint getSampleRate() const {
    /* Return half of the SampleRate, we upSample x2 the F() excitation vector since 
    well need to access F(t+0.5*h) during RK4 algorithm. */
    return sampleRate/2;
  }
  uint getSimDuration() const { return simDuration; }
  MBKMatrix getMBK() const { return MBK; }
  FShapeMatrix getFS() const { return FS; }

};

void Lattice::InitialiseFSvector(){
  FS.F.resize((sampleRate*simDuration)+1);

  // Fill the F (excitation) vector
  for(size_t i(0); i < (sampleRate*simDuration)+1; ++i){
    // FS.F[i] = std::cos(2*PI*excitationFrequency*((reel)i/(sampleRate)));
    FS.F[i] = 0.;
  }

  // Initialize the Shape matrix as well, this matrix define on wich elements
  // the excitation forces will be applied.

  FS.S.resize(numberOfElements);

  // Filling the S "Shape" Matrix with a "border" patern
  for(size_t i(0); i < latticeSize; ++i){
    if(i == 0 || i == latticeSize-1){
      // Filling the top and bottom level of the lattice
      for(size_t j(0); j < latticeSize; ++j){
        FS.S[i*latticeSize+j] = 1;
      }
    }
    else{
      // Filling the interleaved layer with border excitation
      FS.S[i*latticeSize] = 2;
      FS.S[i*latticeSize+latticeSize-1] = 2;
    }
  }
}

void Lattice::InitialiseMBKMatrix(){

  // Filling M Matrix
  MBK.M.resize(numberOfElements, numberOfElements);
  MBK.M.setIdentity();

  // Filling B Matrix
  MBK.B.resize(numberOfElements, numberOfElements);
  MBK.B.setIdentity();
  MBK.B = 0.01 * MBK.B;

  // Filling K Matrix
  MBK.K.resize(numberOfElements, numberOfElements);

  
  // TODO: HERE REAL VALUE FROM THE JUPYTER OF MARC
  /* reel minLocalK = 4.07e9;
  reel maxLocalK = 4.4e9;
  reel minCouplingK = 6.9e7;
  reel maxCouplingK = 88e7; */

  // Hope just for debug ?
  reel minLocalK = 8;
  reel maxLocalK = 12;
  reel minCouplingK = 0.08;
  reel maxCouplingK = 0.12;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> disLocal(minLocalK, maxLocalK);

  // TODO: Re-implement random values after debuging?
  // Will be needed when we'll have a look at the entire optimization process

  // 1. Filling the diagonal (local stifness)
  for(size_t i(0); i < numberOfElements; ++i){
    // MBK.K(i, i) = disLocal(gen);
    MBK.K(i, i) = 1;
  }

  // 2. Filling the remaining coef of the D matrix (coupling)
  if(numberOfElements > 2){
    std::uniform_real_distribution<> disCoupling(minCouplingK, maxCouplingK);
    for(size_t i(0); i < numberOfElements-1; ++i){
      // MBK.K(i, i+1) = disCoupling(gen);
      MBK.K(i, i+1) = 0.05;
    }

    for(size_t j(0); j < numberOfElements-1; ++j){
      // MBK.K(j+1, j) = disCoupling(gen);
      MBK.K(j+1, j) = 0.05;
    }

    // 3. Filling the identity by bloc coef (coupling)
    for(size_t i(0); i < numberOfElements-latticeSize; ++i){ //latticeSize
      // MBK.K(i, i+latticeSize) = disCoupling(gen);
      MBK.K(i, i+latticeSize) = 0.05;
    }

    for(size_t j(0); j < numberOfElements-latticeSize; ++j){ //latticeSize
      // MBK.K(j+latticeSize, j) = disCoupling(gen);
      MBK.K(j+latticeSize, j) = 0.05;
    }
  }
}

std::ostream& Lattice::affiche(std::ostream& out) const{
  std::cout << "Lattice of " << numberOfElements << " elements:" << std::endl;
  std::cout << MBK << std::endl;
  std::cout << FS << std::endl;

  return out;
}

std::ostream& operator<<(std::ostream& sortie, Lattice const& lattice){
  return lattice.affiche(sortie);
}

void Lattice::toFile() const{
  std::ofstream ofs("SIM_F.csv");

  ofs << "F(t) = cos(2*PI*f*t)" << std::endl;

  for(size_t i(0); i < (sampleRate*simDuration)+1; ++i){
    ofs << FS.F[i] << std::endl;
  }
}




/**
 * @brief Define the Runge-Kutta 4rst order algorithm.
 * 
 */
class RK4{
 private:
  uint numberOfSteps;
  reel h; // RK4 time-step

  // RK4 Internal parameter
  MatrixXd X_output;

  VectorXd Q1; // Vector of first-part decomposed ODE (for Krk matrix computation)
  VectorXd Q2; // Vector of second-part decomposed ODE (for Mrk matrix computation)

  MatrixXd Mrk; // RK4 M matrix coefficients
  MatrixXd Krk; // RK4 K Matrix coefficients
  Matrix4d ButcherRK4; // RK4 Butcher tab

  // Simulation parameter
  MBKMatrix MBK; // WARNING: At initialization we invert the M Matrix since it's how it appears in the calculation
  FShapeMatrix FS;
  uint DOF;
  
  void RKStep(size_t t);

 public:
  RK4(Lattice lattice_) : 
    MBK(lattice_.getMBK()),
    FS(lattice_.getFS()),
    DOF(lattice_.getNDOF()),
    numberOfSteps(lattice_.getSampleRate() * lattice_.getSimDuration()),
    h(1.*lattice_.getSimDuration()/numberOfSteps)
    {

      ButcherRK4 << 0,    0,    0, 0,
                    1./2, 0,    0, 0,
                    0,    1./2, 0, 0,
                    0,    0,    1, 0;

      //std::cout << ButcherRK4.row(0) << std::endl;

      MBK.M.inverse();

      X_output.resize(DOF, numberOfSteps);

      Q1.resize(DOF); // Initialized with 0
      Q2.resize(DOF);

      Q1(0) = 1;

      Mrk.resize(4, DOF);
      Krk.resize(4, DOF);
    } 

  void Solve();
  void toFile() const;

};

void RK4::Solve(){

  uint completed = 0;
   
  std::cout << "\r" << "RK4 PROGRESSION: " << completed * 10 <<"%   " << std::flush;

  for(size_t i(0); i < numberOfSteps; ++i){
    RKStep(i);
    
    if(i%((numberOfSteps/10)+1) == 0){
      std::cout << "\r" << "RK4 PROGRESSION: " << completed * 10 <<"%" << std::flush;
      ++completed;
    }
    
  }

  std::cout << "\r" << "RK4 PROGRESSION: " << completed * 10 <<"%" << std::endl; 
 }

void RK4::RKStep(size_t t){

  // Compute the 4 coeficient stages of RK4 algorithm
  // Stage 1
  Mrk.row(0) = Q2.transpose() + ButcherRK4.row(0) * Krk ;
  Krk.row(0) = (MBK.M * (
              - MBK.B * (Q2.transpose() + ButcherRK4.row(0) * Mrk ).transpose() 
              - MBK.K * (Q1.transpose() + ButcherRK4.row(0) * Krk ).transpose() 
              + FS.S*FS.F[2*t])).transpose();
  /* Krk.row(0) = (MBK.M * (
              - MBK.B * Q2 
              - MBK.K * Q1 
              + FS.S*FS.F[2*t])).transpose(); */

  // Stage 2
  Mrk.row(1) = Q2.transpose() + ButcherRK4.row(1) * Krk;
  Krk.row(1) = (MBK.M /*Already inverted*/ * (
              - MBK.B * (Q2.transpose() + ButcherRK4.row(1) * Mrk).transpose() 
              - MBK.K * (Q1.transpose() + ButcherRK4.row(1) * Krk).transpose() 
              + FS.S*FS.F[2*t+1])).transpose();
  
  // Stage 3
  Mrk.row(2) = Q2.transpose() + ButcherRK4.row(2) * Krk;
  Krk.row(2) = (MBK.M /*Already inverted*/ * (
              - MBK.B * (Q2.transpose() + ButcherRK4.row(2) * Mrk).transpose() 
              - MBK.K * (Q1.transpose() + ButcherRK4.row(2) * Krk).transpose() 
              + FS.S*FS.F[2*t+1])).transpose();
  
  
  Mrk.row(3) = Q2.transpose() + ButcherRK4.row(3) * Krk;
  Krk.row(3) = (MBK.M /*Already inverted*/ * (
              - MBK.B * (Q2.transpose() + ButcherRK4.row(3) * Mrk).transpose() 
              - MBK.K * (Q1.transpose() + ButcherRK4.row(3) * Krk).transpose() 
              + FS.S*FS.F[2*t+2])).transpose();
   
  // Compute next Q1 and Q2 vectors
  Q1 += h/6 * ( Mrk.row(0) + 2*Mrk.row(1) + 2*Mrk.row(2) + Mrk.row(3) );
  Q2 += h/6 * ( Krk.row(0) + 2*Krk.row(1) + 2*Krk.row(2) + Krk.row(3) );
  
  X_output.col(t) = Q1;
}

void RK4::toFile() const{
  std::ofstream ofs("SIM_X.csv");

  for(size_t i(0); i < DOF; ++i){
    ofs << "X_" << i << ","; 
  }
  ofs << std::endl;
  
  for(size_t i(0); i < numberOfSteps; ++i){
    for(size_t j(0); j < DOF; ++j){
      ofs << X_output(j, i) << ",";
    }
    ofs << std::endl;
  }
}



/****************************
            MAIN
 ****************************/

int main(int argc,char**argv)
{
  bool verbose(false);
  bool writeOutput(false);
  bool benchmarking(false);

  uint latticeSize(1);
  uint fFrequency(72000);
  uint pointsPerPeriode(4);
  uint simDuration(1);

  // Parsing the input parameter
  for(int i=0;i<argc;i++)
  {
    if(argv[i] == std::string("-verbose"))
      verbose = true;

    if(argv[i] == std::string("-out"))
      writeOutput = true;

    if(argv[i] == std::string("-benchmarking"))
      benchmarking = true;

    if(argv[i] == std::string("-latticesize"))
      latticeSize = std::stol(argv[i+1]);

    if(argv[i] == std::string("-ffrequency"))
      fFrequency = std::stol(argv[i+1]);

    if(argv[i] == std::string("-pointsPerPeriode"))
      pointsPerPeriode = std::stol(argv[i+1]);

    if(argv[i] == std::string("-simDuration"))
      simDuration = std::stol(argv[i+1]);
  } 

  // Creating the lattice
  Lattice lattice(latticeSize, fFrequency, pointsPerPeriode, simDuration);
  if(verbose){
    std::cout << lattice << std::endl;
  }

  // Performe the simulation
  RK4 solver(lattice);

  if(benchmarking){
    auto begin = std::chrono::high_resolution_clock::now(); 

    solver.Solve(); 

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
  
    std::cout << "Benchmarked for:" << std::endl;
    std::cout << "  latticesize: " << latticeSize << std::endl;
    std::cout << "  fFrequency: " << fFrequency << std::endl;
    std::cout << "  pointsPerPeriode: " << pointsPerPeriode << std::endl;
    std::cout << "  time: " << elapsed.count()/1000000000.0 << " [s]" << std::endl;
  }
  else{
    solver.Solve();
  }

  if(writeOutput){
    lattice.toFile();
    solver.toFile();
  }
  
  return 0;
}