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
#include <string>

#include <Eigen/Dense>

#define PI 3.14159265359
using reel = double;

typedef Eigen::Matrix<reel, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::DiagonalMatrix< reel, Eigen::Dynamic> DiagonalMatrixXd;
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
 public:
  Lattice(uint latticeSize_, reel sampleRate_, uint numberOfSteps_, reel frequencyOfExcitationSource_): 
    latticeSize(latticeSize_),
    sampleRate(sampleRate_),
    numberOfSteps(numberOfSteps_),
    frequencyOfExcitationSource(frequencyOfExcitationSource_) 
    {
      numberOfElements = latticeSize*latticeSize;

      InitialiseFSvector();
      InitialiseMBKMatrix();
    }

  std::ostream& affiche(std::ostream& out) const;
  void toFile() const;

  uint getSize() const { return latticeSize; }
  uint getNDOF() const { return numberOfElements; }
  reel getTimeStep() const { return 1./sampleRate; }
  uint getNumberOfRKStepsToPerfome() const { return numberOfSteps/2; }
  MBKMatrix getMBK() const { return MBK; }
  FShapeMatrix getFS() const { return FS; }

 private:
  reel sampleRate; // [Hz]
  uint numberOfSteps; // [N]
  reel frequencyOfExcitationSource; // [Hz]

  uint latticeSize; // Size of the lattice across 1 dimmension
  uint numberOfElements; // Total number of elements in the lattice

  MBKMatrix MBK;
  FShapeMatrix FS;

  void InitialiseFSvector();
  void InitialiseMBKMatrix();
};

void Lattice::InitialiseFSvector(){
  FS.F.resize(numberOfSteps+2);

  reel numberOfPeriodeToSimulate( (numberOfSteps+2)*frequencyOfExcitationSource/sampleRate );

  // Fill the F (excitation) vector
  for(size_t i(0); i < numberOfSteps+2; ++i){
    FS.F[i] = std::cos(2.*PI*frequencyOfExcitationSource*( (reel)i*numberOfPeriodeToSimulate/((numberOfSteps+2)*frequencyOfExcitationSource) ));
    // FS.F[i] = 0.;
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
  MBK.B = 0.1 * MBK.B;

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

  for(size_t i(0); i < numberOfSteps; ++i){
    ofs << FS.F[i] << std::endl;
  }
}




/**
 * @brief Define the Runge-Kutta 4rst order algorithm.
 * 
 */
class RK4{
 public:
  RK4(Lattice lattice_) : 
    MBK(lattice_.getMBK()),
    FS(lattice_.getFS()),
    DOF(lattice_.getNDOF()),
    numberOfSteps(lattice_.getNumberOfRKStepsToPerfome()),
    h(lattice_.getTimeStep()),
    h6(h/6.),
    h2(h/2.)
    {
      std::cout << h << std::endl;
      // Pre-compute M^-1 -> B, K F matrix
      MBK.B = MBK.M.inverse()*MBK.B;
      MBK.K = MBK.M.inverse()*MBK.K;
      FS.S = MBK.M.inverse()*FS.S;

      X_output.resize(DOF, numberOfSteps);

      Q1.resize(DOF); // Initialized with 0
      Q2.resize(DOF);

      // Q1(0) = 1;

      Mrk_1.resize(DOF);
      Mrk_2.resize(DOF);
      Mrk_3.resize(DOF);
      Mrk_4.resize(DOF);
      
      Krk_1.resize(DOF);
      Krk_2.resize(DOF);
      Krk_3.resize(DOF);
      Krk_4.resize(DOF);

      Mrk_i.resize(DOF);
      Krk_i.resize(DOF);
    } 

  void Solve();
  void toFile() const;

 private:
  uint numberOfSteps;
  reel h; // RK4 time-step
  reel h6; // = h/6
  reel h2; // = h/2

  // RK4 Internal parameter
  MatrixXd X_output;

  VectorXd Q1; // Vector of first-part decomposed ODE (for Krk matrix computation)
  VectorXd Q2; // Vector of second-part decomposed ODE (for Mrk matrix computation)

  // RK4 M matrix coefficients
  VectorXd Mrk_1;
  VectorXd Mrk_2;
  VectorXd Mrk_3;
  VectorXd Mrk_4;
  
  // RK4 K Matrix coefficients
  VectorXd Krk_1;
  VectorXd Krk_2;
  VectorXd Krk_3;
  VectorXd Krk_4;

  VectorXd Mrk_i;
  VectorXd Krk_i;

  // Simulation parameter
  MBKMatrix MBK; // WARNING: At initialization we invert the M Matrix since it's how it appears in the calculation
  FShapeMatrix FS;
  uint DOF;
  
  void RKStep(size_t t);
  void derivatives(VectorXd& targetK, VectorXd& targetM, size_t t);
};

void RK4::Solve(){

  uint completed = 0;
   
  std::cout << "\r" << "RK4 PROGRESSION: " << completed * 10 <<"%   " << std::flush;

  // std::cout << "GEFE: " << numberOfSteps << std::endl;

  for(size_t i(0); i < numberOfSteps; ++i){
    RKStep(i);
    
    if(i%((numberOfSteps/10)+1) == 0){
      std::cout << "\r" << "RK4 PROGRESSION: " << completed * 10 <<"%" << std::flush;
      ++completed;
    }
    
  }

  std::cout << "\r" << "RK4 PROGRESSION: " << completed * 10 <<"%" << std::endl; 
 }

inline void RK4::derivatives(VectorXd & targetM, VectorXd & targetK, size_t t) {
  targetM = Krk_i;
  targetK = -MBK.B * Krk_i -MBK.K * Mrk_i + FS.S*FS.F[t];
}

inline void RK4::RKStep(size_t t){

  // Compute the 4 coeficient stages of RK4 algorithm
  derivatives(Mrk_1, Krk_1, 2*t);

  Mrk_i = Q1 + h2*h*Mrk_1;
  Krk_i = Q2 + h2*h*Krk_1;

  derivatives(Mrk_2, Krk_2, 2*t+1);

  Mrk_i = Q1 + h2*h*Mrk_2;
  Krk_i = Q2 + h2*h*Krk_2;

  derivatives(Mrk_3, Krk_3, 2*t+1);

  Mrk_i = Q1 + h2*Mrk_2;
  Krk_i = Q2 + h2*Krk_2;

  derivatives(Mrk_4, Krk_4, 2*t+2);

  // Compute next Q1 and Q2 vectors
  Q1 += h6 * ( Mrk_1 + 2*Mrk_2 + 2*Mrk_3 + Mrk_4 );
  Q2 += h6 * ( Krk_1 + 2*Krk_2 + 2*Krk_3 + Krk_4 );
  
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
  reel sampleRate(100);
  reel frequencyOfExcitationSource(1);
  uint numberOfPeriodeToSimulate(1);
  uint numberOfSteps(sampleRate*(numberOfPeriodeToSimulate/frequencyOfExcitationSource));

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

    if(argv[i] == std::string("-sampleRate"))
      sampleRate = std::stod(argv[i+1]);

    if(argv[i] == std::string("-numberOfSteps"))
      numberOfSteps = std::stoul(argv[i+1]);

    if(argv[i] == std::string("-frequencyOfExcitationSource"))
      frequencyOfExcitationSource = std::stod(argv[i+1]);
  } 
  //std::cout << "ee " << frequencyOfExcitationSource << std::endl;
  // Creating the lattice
  Lattice lattice(latticeSize, sampleRate, numberOfSteps, frequencyOfExcitationSource);
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
    std::cout << "  sampleRate: " << sampleRate << std::endl;
    std::cout << "  numberOfSteps: " << numberOfSteps << std::endl;
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