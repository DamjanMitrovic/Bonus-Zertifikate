#include <boost/lambda/lambda.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include <stdio.h> // For printf
#include <stdbool.h> // For boolian
#include <cstdlib>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

double Evaluate_BonusDerivat_loop(const double& S_value, const int& Timesteps_left, bool TouchBarrier);

double Mean_and_discount(const double& Price_Sum, const int& Timesteps_left); 

double Evaluate_Certificate_Price(const double& S_value, const int& Timesteps_left, bool TochBarrier); 

double Propagate_S(const double& S_value, bool& TouchBarrier); 

double Get_StdNormDistr_RandomNumber();

template<std::size_t size>
double Std_Dev(const std::array<double, size>& Arr, const double& Mean); 

template<std::size_t size>
double Ar_Mean(const std::array<double, size>& Arr);

std::random_device rd{};     // Random seed from the computer; Syntax is equal to: std::random_device rd;
//std::mt19937 gen{rd()};    // Method with which the random numbers are generated, Therefore, one needs to provide the seed at which to start (here rd or arbitrary number), Syntax equal to std::mt19937 gen = 2;
std::mt19937 gen{2};         // Method with which the random numbers are generated, Therefore, one needs to provide the seed at which to start (here 2), Syntax equal to std::mt19937 gen = 2;
boost::math::normal dist(0.0,1.0); 

//global constant variables
const int S_0_Min_ = 50;
const int S_0_Max_ = 140;
const int S_0_Stepsize_ = 1;
const int N_Paths_S_      = 100;              // Number of Paths for S
const int N_Paths_ 	     = 1000;      	   		// Number of paths for MC evaluating the Certificate Price
const int N_Evaluations_ = 250;            // Number of Value evaluations (rate of return timesteps)
const int TS_per_eval_	= 1;              // Timesteps per evaluation
const int N_Timesteps_    = N_Evaluations_ * TS_per_eval_; // Number of time steps that are calculated
const double Maturity_       = 1.0;    		// Maturity of the derivative in [a]
const double Timestep_size_ = Maturity_ / N_Timesteps_; // time intervall for each evaluation
//	const double S_0_;        // = 72.93;   // Initial value of the underlying // https://www.onvista.de/aktien/BASF-Aktie-DE000BASF111 - not needed as const
//	const double Z_0_;	//= 93.93; //given from calculation in group task - not needed as const
const double Barrier_     = 50;      // Barrier of the Bonus certificate // https://www.onvista.de/derivate/bonus-zertifikate/BONUS-ZERTIFIKAT-AUF-BASF-DE000DC3YKS7
const double Bonus_level_ = 100;     // Bonuslevel of the Bonus certificate // https://www.onvista.de/derivate/bonus-zertifikate/BONUS-ZERTIFIKAT-AUF-BASF-DE000DC3YKS7
const double mu_ = 0.0855;
const double r_f_         = -0.0055; // Estimated risk free interest rate
														    // Einjährige deutsche staatsanleihe mit Laufzeit 1a und Ausgabedatum 16/04/19 r_f = -0.550% //https://de.investing.com/rates-bonds/germany-1-year-bond-yield
														    // LIBOR r_f = -0,2113 % 	// https://www.finanzen.net/zinsen/historisch/libor/libor-eur-12-monate
														    // Mittelwert aus Staatsanleihe und LIBOR: r_f = -0.38065%
const double Sigma_       = 0.2429;  // Estimated volatility of the underlying BASF11 aktie //https://www.onvista.de/aktien/BASF-Aktie-DE000BASF111
														 		// Standartabweichung geschätzt aus Daten ab 03.02.2010. Aus Daten stetige Renditen ln(S_t1/S_t) berechnet. 
														 		// Varianz über stetige Tagesrenditen berechnet und dann mit 250 tagen multipliziert


int main(int argc, char const *argv[])
{
  double S_0;
  double S_old;
  double S_new;

  double Z_0;
  double Z_old;
  double Z_new;
  double Z_sum_new;

  double Rendite;
  double Mean_Renditen;
	double Volatility;
	double Sharp_Ratio;

  int N_Timesteps_left;

	bool TouchBarrier;

	std::array<double, N_Evaluations_> Renditen;
	
	//only works on linux...creating output directories
	if (system("mkdir Output") == -1) {
		printf("Output Directories not created\n");
		exit(1);
	}

    
	std::ofstream volatility;
	std::ofstream sharp_ratio;
	volatility.open("Output/Volatility.csv");
	sharp_ratio.open("Output/Sharp_Ratio.csv");

	for (int start_price = S_0_Min_; start_price <= S_0_Max_; start_price = start_price+S_0_Stepsize_){
		volatility << start_price << ",";
		sharp_ratio << start_price << ",";
	}
	volatility << std::endl;
	sharp_ratio << std::endl;

	//Different paths for S-Propagation
	for (int path_s = 1; path_s <= N_Paths_S_; path_s++){
  	std::cout << "Path_S: " << path_s << std::endl;
	
	  for (int start_price = S_0_Min_; start_price <= S_0_Max_; start_price = start_price+S_0_Stepsize_){
			std::cout << "S_0: " << start_price << std::endl;
			S_0 = double(start_price);

			//Set TouchBarrier = false upfront (before 1. calculation of the value) 
			//If starting underlying price is already hitting barrier set it to true
	  	TouchBarrier = false;
			if (start_price <= Barrier_){
				TouchBarrier = true;
			}

	    Z_0 = Evaluate_Certificate_Price(S_0, N_Timesteps_, TouchBarrier);

	    S_old = S_0;
	    Z_old = Z_0;


			for (int ts = 1; ts <= N_Evaluations_; ts++){

	    	S_new  = Propagate_S(S_old, TouchBarrier);
				N_Timesteps_left = N_Timesteps_ - ts * TS_per_eval_;
       
	      Z_new = Evaluate_Certificate_Price(S_new, N_Timesteps_left, TouchBarrier);
        
	      //Log returns on daily basis...can be modified by using less timesteps
	      Rendite = log(Z_new / Z_old);
				//Fill array with log returns (steady_rendity)
				Renditen[ts-1] = Rendite;

	      //Preparing for next Timestep
				Z_old = Z_new;
	      S_old = S_new;

			} // ending loop over evaluation timesteps
			//Calculate Mean of Rate of Return
			Mean_Renditen = Ar_Mean(Renditen);
			//Annualization of the volatiliy by multiplying by the square-root of days
			Volatility = sqrt(N_Evaluations_) * Std_Dev(Renditen, Mean_Renditen);
			//Annualisieren of the rate of return by multiplying by days (Evaluation_Points)
			Mean_Renditen = N_Evaluations_ * Mean_Renditen;
			//Sharp_Ratio
			Sharp_Ratio = Mean_Renditen / Volatility;
			//Output to files
			volatility << Volatility << ",";
			sharp_ratio << Sharp_Ratio << ",";
		} // endling loop over outer paths
		volatility << std::endl;
		sharp_ratio << std::endl;
	} // ending loop over starting price
	volatility.close();
	sharp_ratio.close();
	return 0;
}


//Important that the bool is passed by value here
//As it is changed to true if the barrier is hit but only for this path
//For the outer path it should remain as it was
double Evaluate_BonusDerivat_loop(const double& S_value, const int& Timesteps_left, bool TouchBarrier){
	
	double S_tmp = S_value;
	double random_number;

	for(int tm_step = 1; tm_step <= Timesteps_left; tm_step++){
		random_number = Get_StdNormDistr_RandomNumber(); // Get the standard normal distributed random number                                                                             
		S_tmp = S_tmp * exp( (r_f_ - (Sigma_*Sigma_)/2.) * Timestep_size_ + Sigma_ * sqrt(Timestep_size_) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX

		//Check if barrier was touched
		if(S_tmp <= Barrier_) TouchBarrier = true;                                        // If barrier was touched set the boolian to kTRUE 
 	}

  //Determine the correct price of the bonus certificate
	if(TouchBarrier || (S_tmp > Bonus_level_) )  return S_tmp;                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
	else return Bonus_level_;																												  // If the value of the underlying has never touched or undershoot the Barrier and if the value of the underlying is smaller than the Bonuslevel than return the Bonuslevel = Price of Bonus certificate
}


double Mean_and_discount(const double& Price_Sum, const int& Timesteps_left) {
  double discount_time = Timesteps_left * Timestep_size_;
  double Discounted_Price_Zertifikat = Price_Sum * exp(-r_f_ * discount_time) / N_Paths_;

	return Discounted_Price_Zertifikat;
}


double Evaluate_Certificate_Price(const double& S_value, const int& Timesteps_left, bool TouchBarrier) {

				double Z_sum = 0.0; //set them to zero before adding values for monte_carlo
	      for (int path = 0; path < N_Paths_; path++){

	        Z_sum += Evaluate_BonusDerivat_loop(S_value, Timesteps_left, TouchBarrier);

	      } //ending loop over inner paths
       
	      double Z = Mean_and_discount(Z_sum, Timesteps_left);

				return Z;
}


//Important that the bool is not a const reference here...as variable should be changed if underlying hits barrier
//Propagate TS_per_eval timesteps until next evaluation
double Propagate_S(const double& S_value, bool& TouchBarrier) {
  double random_number;
	double S_t;
  double S_t_old = S_value;
  for (int inter_tm_step = 1; inter_tm_step <= TS_per_eval_; inter_tm_step++){
    
		random_number = Get_StdNormDistr_RandomNumber(); // Get the standard normal distributed random number
    //Taking into account that estimated mu_ comes from steady returns so that "mu_" - (Sigma^2)/2 = mu_ in this formula
	  S_t = S_t_old * exp( mu_ * Timestep_size_ + Sigma_ * sqrt(Timestep_size_) * random_number );
    S_t_old = S_t;
	  //Check if barrier was touched
		if(S_t <= Barrier_) TouchBarrier = true;
  }  
  return S_t;
}


double Get_StdNormDistr_RandomNumber(){
  double random_number;
  std::normal_distribution<> nd{0,1}; // Define normal distribution with name nd and mean = 0 and sigma =1 which is equal to standard normal distribution; Return type is double
	random_number = nd(gen);           // Generate random number which is standard normal distributed

	return random_number;
}


template<std::size_t size>
double Ar_Mean(const std::array<double, size>& Arr) {
	double Mean = 0;
	
	for (auto& Val : Arr) {
		Mean += Val;
	}
	Mean /= size;
	return Mean;
}

//template to calculate the error...currently unused
template<std::size_t size>
double Std_Dev(const std::array<double, size>& Arr, const double& Mean) {
	
	double Var = 0;
  double Err;

	for (auto& Val : Arr) {
    Var += pow(Val - Mean, 2);
	}

  Var /= (size - 1);
    
  Err = sqrt(Var);
	return Err;
}


