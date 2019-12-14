#include <boost/lambda/lambda.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include <stdio.h> // For printf
#include <stdbool.h> // For boolian
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>



std::random_device rd{};     // Random seed from the computer; Syntax is equal to: std::random_device rd;
//std::mt19937 gen{rd()};    // Method with which the random numbers are generated, Therefore, one needs to provide the seed at which to start (here rd or arbitrary number), Syntax equal to std::mt19937 gen = 2;
std::mt19937 gen{2};         // Method with which the random numbers are generated, Therefore, one needs to provide the seed at which to start (here 2), Syntax equal to std::mt19937 gen = 2;
boost::math::normal dist(0.0,1.0); 




int main(int argc, char const *argv[])
{
	const int N_Paths_S_      = 10;              // Number of Paths for S
	const int N_paths_ 	     = 10;      	   		// Number of paths for MC evaluating the Certificate Price
	const int N_timeSteps_    = 100;        	 		// Number of time steps that are calculated
	const int maturity_       = 1;    		// Maturity of the derivative in [a]
	const double TimeStepSize_ = (double)maturity/N_timeSteps; // time intervall for each evaluation
	const double S_0_         = 72.93;   // Initial value of the underlying // https://www.onvista.de/aktien/BASF-Aktie-DE000BASF111
	const double Barrier_     = 50;      // Barrier of the Bonus certificate // https://www.onvista.de/derivate/bonus-zertifikate/BONUS-ZERTIFIKAT-AUF-BASF-DE000DC3YKS7
	const double Bonus_Level_ = 100;     // Bonuslevel of the Bonus certificate // https://www.onvista.de/derivate/bonus-zertifikate/BONUS-ZERTIFIKAT-AUF-BASF-DE000DC3YKS7
	const double mu_ = 0.0855;
	const double epsilon_ = 0.0001;
	const double r_f_         = -0.0055; // Estimated risk free interest rate
														    // Einjährige deutsche staatsanleihe mit Laufzeit 1a und Ausgabedatum 16/04/19 r_f = -0.550% //https://de.investing.com/rates-bonds/germany-1-year-bond-yield
														    // LIBOR r_f = -0,2113 % 	// https://www.finanzen.net/zinsen/historisch/libor/libor-eur-12-monate
														    // Mittelwert aus Staatsanleihe und LIBOR: r_f = -0.38065%
	const double Sigma_       = 0.2429;  // Estimated volatility of the underlying BASF11 aktie //https://www.onvista.de/aktien/BASF-Aktie-DE000BASF111
														 		// Standartabweichung geschätzt aus Daten ab 03.02.2010. Aus Daten stetige Renditen ln(S_t1/S_t) berechnet. 
														 		// Varianz über stetige Tagesrenditen berechnet und dann mit 250 tagen multipliziert
  double S;
  double S_old;
  double S_new;
  double S_new_epsilon;

  double Z;
  double Z_eps;
  double Z_sum;
  double Z_eps_sum;

  int N_timeSteps_left;
	bool TouchBarrier;

	for (int path_s = 0; path_s < N_Paths_S_; path_s++){

		TouchBarrier_ = false;
    S = S_0_;

		for (int ts = 1; ts <= N_timeSteps_; ts++){

			  S_old = S;
        S_new  = Propagate_S(S, Sigma_, mu_, Barrier_, 
        	               Bonus_Level_, TimeStepSize_, TouchBarrier);
        S_new_epsilon = S_new_ + epsilon_;
				N_timeSteps_left = N_timeSteps_ - ts;

				Z_sum = 0; //set them to zero before adding values for monte_carlo
        Z_eps_sum = 0;
      for (int path = 0; path < N_paths_; path++){

      	Z_sum += Evaluate_BonusDerivat_loop(S_new_, Sigma_, r_f_, Barrier_, Bonus_Level_, 
      		                         TimeStepSize_, N_timeSteps_left, TouchBarrier);
      	Z_eps_sum += Evaluate_BonusDerivat_loop(S_new_epsilon_, Sigma_, r_f_, Barrier_, Bonus_Level_, 
      		                         TimeStepSize_, N_timeSteps_left, TouchBarrier);
        //Calculate Mean
      	Z = N_paths_;
      	Z_eps/= N_paths_;
        
        //Discount
      }
        
      //Delta berechnen
      //Renditen berechnen
      //Output
		}
		
	}


	return 0;
}

double Mean_and_discount(const double& Price_Sum, const double& r_f const int& N_timeSteps, 
												 const int& N_timeSteps_left, const) {
	double dicount_time = N_timeSteps_left * 
  double Discounted_Price_Zertifikat = Price_Sum * exp(-r_f * do) / N_paths_;

	return Discounted_Price_Zertifikat;
}

//Important that the bool is passed by value here
double Evaluate_BonusDerivat_loop(const double& S_value, const double& Sigma, 
																	const double& r_f, const double& Barrier, 
																	const double& Bonus_Level, const double& delta_t,
																	const int& TimeSteps_left, bool TouchBarrier){
	
	double S_tmp = 0;
	double S_tmp_initial = S_value;
	double random_number;

	for(int tm_step = 1; tm_step <= TimeSteps_left; tm_step++){
		random_number = Get_StdNormDistr_RandomNumber(); // Get the standard normal distributed random number                                                                             
		S_tmp = S_tmp_initial * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
		S_tmp_initial = S_tmp;

		//Check if barrier was touched
		if(S_tmp <= Barrier) TouchBarrier = true;                                        // If barrier was touched set the boolian to kTRUE 
 	}

  //Determine the correct price of the bonus certificate
	if(TouchBarrier || (S_tmp > Bonus_Level) )  return S_tmp;                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
	else return Bonus_Level;																												  // If the value of the underlying has never touched or undershoot the Barrier and if the value of the underlying is smaller than the Bonuslevel than return the Bonuslevel = Price of Bonus certificate
}

//Important that the bool is not a const reference here
double Propagate_S(const double& S_value, const double& Sigma, 
									 const double& mu, const double& Barrier, 
									 const double& Bonus_Level, const double& delta_t,
									 bool& TouchBarrier) {

	random_number = Get_StdNormDistr_RandomNumber(); // Get the standard normal distributed random number
  double S_t = S_value * exp( (mu - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );
  //Check if barrier was touched
	if(S_t <= Barrier) TouchBarrier = true;

  return S_t;

}

double Get_StdNormDistr_RandomNumber(){

  std::normal_distribution<> nd{0,1}; // Define normal distribution with name nd and mean = 0 and sigma =1 which is equal to standard normal distribution; Return type is double
	double random_number = nd(gen);           // Generate random number which is standard normal distributed

	return random_number;
}