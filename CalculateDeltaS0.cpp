
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


struct PriceBonusDerivat_pls_Delta2{
	double priceDerivat_S0;
	double priceDerivat_S0_pls_epsilon;
};


double Get_StdNormDistr_RandomNumber();

void Calculate_Delta(const double S_0, const double Sigma, const double Barrier, const double Bonus_Level, const double mu, const double r_f, const double epsilon,
										 const double maturity, const int N_paths, const int N_timeSteps);
PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta21(const double S_0, const double S_0_pls_epsilon, const double Sigma, 
																																		const double r_f, const double Barrier, const double Bonus_Level, const double maturity, const int N_timeSteps,
																																		const int N_inner_paths, int current_time, bool Barrier_already_touched_S0,
																																		bool Barrier_already_touched_S0_pls_epsilon,
																																		bool DividendsPaid, const double Dividend );
// Global variables

int tmp_step_noX_R = 0;	     // Current time step that is calculated 
std::random_device rd{};     // Random seed from the computer; Syntax is equal to: std::random_device rd;
//std::mt19937 gen{rd()};    // Method with which the random numbers are generated, Therefore, one needs to provide the seed at which to start (here rd or arbitrary number), Syntax equal to std::mt19937 gen = 2;
std::mt19937 gen{2};         // Method with which the random numbers are generated, Therefore, one needs to provide the seed at which to start (here 2), Syntax equal to std::mt19937 gen = 2;
boost::math::normal dist(0.0,1.0); 

int main(){

	//printf("Hello World\n");

	int N_paths 	 ;      	   		// Number of paths that are calculated
	int N_timeSteps;        	 		// Number of time steps that are calculated
	double maturity       = 2.;    		// Maturity of the derivative in [a]
	double S_0         = 72.93;   // Initial value of the underlying // https://www.onvista.de/aktien/BASF-Aktie-DE000BASF111

	double Barrier     = 50;      // Barrier of the Bonus certificate // https://www.onvista.de/derivate/bonus-zertifikate/BONUS-ZERTIFIKAT-AUF-BASF-DE000DC3YKS7
	double Bonus_Level = 100;     // Bonuslevel of the Bonus certificate // https://www.onvista.de/derivate/bonus-zertifikate/BONUS-ZERTIFIKAT-AUF-BASF-DE000DC3YKS7
	double K 					 = 0.0;     // Strike price for the Call which is used to rebild the Bonus certificate
	double CoC_q 			 = 0.0;     // Discountfactor for dividends "Cost of carry"
	double r_f         = -0.0055;   // Estimated risk free interest rate
														    // Einjährige deutsche staatsanleihe mit Laufzeit 1a und Ausgabedatum 16/04/19 r_f = -0.550% //https://de.investing.com/rates-bonds/germany-1-year-bond-yield
														    // LIBOR r_f = -0,2113 % 	// https://www.finanzen.net/zinsen/historisch/libor/libor-eur-12-monate
														    // Mittelwert aus Staatsanleihe und LIBOR: r_f = -0.38065%
	double Sigma       = 0.2429;  // Estimated volatility of the underlying BASF11 aktie //https://www.onvista.de/aktien/BASF-Aktie-DE000BASF111
														 		// Standartabweichung geschätzt aus Daten ab 03.02.2010. Aus Daten stetige Renditen ln(S_t1/S_t) berechnet. 
														 		// Varianz über stetige Tagesrenditen berechnet und dann mit 250 tagen multipliziert

	double mu = 0.0855;
	double epsilon = 0.05;

  Calculate_Delta( S_0,  Sigma,  Barrier,  Bonus_Level,  mu,  r_f,  epsilon,  maturity,  1000, 500 );// M =1 --> Nt = 250, M=2 --> Nt = 500, M=0.5 --> Nt = 125



}//end main function



//_____________________________________________________________________________________________________________________________________________________________________________________________________
//

double Get_StdNormDistr_RandomNumber(){

  std::normal_distribution<> nd{0,1}; // Define normal distribution with name nd and mean = 0 and sigma =1 which is equal to standard normal distribution; Return type is double
	double random_number = nd(gen);           // Generate random number which is standard normal distributed

	return random_number;
}

void Calculate_Delta(const double S_0, const double Sigma, const double Barrier, const double Bonus_Level, const double mu, const double r_f, const double epsilon, const double maturity, const int N_paths, const int N_timeSteps){

	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation

	//only works on linux...creating output directories
	if (system("mkdir Delta_S0_Mat2") == -1) {
		printf("Output Directories not created\n");
		exit(1);
	}

	std::ofstream file_delta;
	file_delta.open("Delta_S0_Mat2/DeltasMat2.csv");

	//Write first row
 	std::vector<double> S0_prices = {40., 45., 46., 47., 48., 48.5, 49., 49.2, 49.4, 49.6, 49.70, 49.8, 49.9, 50.000, 50.10, 50.20, 50.30, 50.40, 50.50, 50.60, 50.70, 50.80, 50.90,
																		51.00, 51.5, 52.00, 52.50, 53.00, 53.50, 54.00, 54.50, 55.00, 56.00, 57.00, 58.00, 59.00, 60.00, 65.0, 70.0, 75., 80., 85., 90., 95., 100.0, 105.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0};

	for(int start_price = 0; start_price <= S0_prices.size() -1; start_price++){

  	file_delta << S0_prices[start_price] << ",";
	}
	file_delta << std::endl;

	for(int path = 1; path <= N_paths; path++){

	double S_tmp_mu = 0.0;
	double S_tmp_mu_pls_epsilon = 0.0;

	double D_tmp = 0.0;
	double D_tmp_pls_epsilon = 0.0;
	double delta = 0.0;

	struct PriceBonusDerivat_pls_Delta2 Struct_2;
	Struct_2.priceDerivat_S0 = 0.0;
	Struct_2.priceDerivat_S0_pls_epsilon = 0.0;

	//bool DividendsPaid = false; // use when dividends are supposed to pay
	bool DividendsPaid = true; // use if dividends should not be used
	double Dividend = 3.30;

	double random_number;
	double time = 0.0;

	printf("___Path: %d   _____________________\n", path);


		for(int start_price = 0; start_price <= S0_prices.size() -1; start_price++){
		  
			bool TouchBarrier_D1             = false;// Boolian for checking whether barrier was touched or undershoot by the value of the certificate
			bool TouchBarrier_D1_pls_epsilon = false;// Boolian for checking whether barrier was touched or undershoot by the value of the certificate plus epsilon

			S_tmp_mu = S0_prices[start_price];
			S_tmp_mu_pls_epsilon = S_tmp_mu + epsilon;


		  //Check if barrier was touched
			if(S_tmp_mu <= Barrier){             
				TouchBarrier_D1 = true;                                        // If barrier was touched set the boolian to kTRUE 
				//printf("Äussere Barrier S_tmp_mu touched ! Path = %d, S_tmp_mu = %f, Bool = %i\n", path, S_tmp_mu, TouchBarrier_D1);
				}
			if(S_tmp_mu_pls_epsilon <= Barrier){
				TouchBarrier_D1_pls_epsilon = true;  
			  //printf("____Äussere Barrier S_tmp_mu_pls_epsilon touched ! Path = %d, S_tmp_mu_pls_epsilon = %f, Bool = %i\n", path, S_tmp_mu_pls_epsilon, TouchBarrier_D1_pls_epsilon);
			}


			Struct_2 = Evaluate_BonusDerivat_loop_for_delta21(S_tmp_mu, S_tmp_mu_pls_epsilon, Sigma, r_f, Barrier, Bonus_Level, maturity, N_timeSteps, 5000, 1, TouchBarrier_D1, TouchBarrier_D1_pls_epsilon, DividendsPaid, Dividend);
			
			D_tmp 						         = Struct_2.priceDerivat_S0;
			D_tmp_pls_epsilon          = Struct_2.priceDerivat_S0_pls_epsilon;
			
			delta = (D_tmp_pls_epsilon - D_tmp )/epsilon;
			file_delta << delta << ",";

		}// end loop through S0s
		file_delta  << "\n";

	}//end loop through paths
	file_delta.close();
}

//_____________________________________________________________________________________________________________________________________________________________________________________________________
//

PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta21(const double S_0, const double S_0_pls_epsilon, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const double maturity, const int N_timeSteps, const int N_inner_paths, int current_time, bool Barrier_already_touched_S0, bool Barrier_already_touched_S0_pls_epsilon, bool DividendsPaid, const double Dividend ){
	
	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation
	double price_S0 = 0.0;
	double price_S0_pls_epsilon = 0.0;
	double tmp_sum_price_S0 = 0.0;
	double tmp_sum_price_S0_pls_epsilon = 0.0;
	// Hier wird der bool value zur Berücksichtigung der Historie auf den Wert bei t_1 gesetzt
	bool TouchBarrier_S0 = Barrier_already_touched_S0;// Boolian for checking whether barrier was touched or undershoot
	bool TouchBarrier_S0_pls_epsilon = Barrier_already_touched_S0_pls_epsilon;// Boolian for checking whether barrier was touched or undershoot
	struct PriceBonusDerivat_pls_Delta2 tmp_struct;

	//printf("1.  S0_init    = %f, current_time = %i, boolian_S0 = %i\n", S_tmp_initial_S0, current_time,TouchBarrier_S0);
	//printf("1.  S0_init_ep = %f, current_time = %i, boolian_ep = %i\n", S_tmp_initial_S0_pls_epsilon, current_time,TouchBarrier_S0_pls_epsilon);

	if(current_time > N_timeSteps){
		printf("Error in function Evaluate_BonusDerivat_loop_for_delta21, current_time > N_timeSteps");
	}
	
	for(int path = 1; path <= N_inner_paths; path++){

		//printf("_______________________Inner Path:  %d_____________\n",path);

		double S_tmp_S0 = S_0;
		double S_tmp_S0_pls_epsilon = S_0_pls_epsilon;
		double S_tmp_initial_S0 = S_0;
		double S_tmp_initial_S0_pls_epsilon = S_0_pls_epsilon;
		double random_number;
		// Hier wird der bool value zur Berücksichtigung der Historie auf den Wert bei t_1 gesetzt
		// Diese Initialisierung wird für jeden pfad durchlauf vorgenommen

		TouchBarrier_S0 = Barrier_already_touched_S0;// Boolian for checking whether barrier was touched or undershoot
		TouchBarrier_S0_pls_epsilon = Barrier_already_touched_S0_pls_epsilon;// Boolian for checking whether barrier was touched or undershoot


		double tmp_price_S0 = 0.0;
		double tmp_price_S0_pls_epsilon = 0.0;
//printf("Innerer Pfad = %d\n", path);
		for(int tmp_step_noX = current_time; tmp_step_noX < N_timeSteps; tmp_step_noX++){

			if(tmp_step_noX > N_timeSteps/2. && !DividendsPaid){
				DividendsPaid = true;
				S_tmp_initial_S0 -= Dividend;
				S_tmp_initial_S0_pls_epsilon -= Dividend;
			}

			random_number = Get_StdNormDistr_RandomNumber();    // Get the standard normal distributed random number                                                                    
			S_tmp_S0             = S_tmp_initial_S0 * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
			S_tmp_S0_pls_epsilon = S_tmp_initial_S0_pls_epsilon * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
			//printf("Verhältnis = %f\n", S_tmp_S0/S_tmp_S0_pls_epsilon);
			

			S_tmp_initial_S0 = S_tmp_S0;
			S_tmp_initial_S0_pls_epsilon = S_tmp_S0_pls_epsilon;


			//Check if barrier was touched
			if(!TouchBarrier_S0 && S_tmp_S0 <= Barrier){
			  TouchBarrier_S0 = true;                                        // If barrier was touched set the boolian to kTRUE 
			  //printf("INNERE Barrier S_tmp_S0 touched ! Path = %d, Time step = %d, S_tmp_S0 = %f, Bool = %i\n", path,tmp_step_noX, S_tmp_S0, TouchBarrier_S0);
			}

			if(!TouchBarrier_S0_pls_epsilon && S_tmp_S0_pls_epsilon <= Barrier) {
				TouchBarrier_S0_pls_epsilon = true;                                        // If barrier was touched set the boolian to kTRUE 
			  //printf("INNERE Barrier S_tmp_S0_pls_epsilon touched ! Path = %d, Time step = %d, S_tmp_S0_pls_epsilon = %f, Bool = %i\n", path,tmp_step_noX, S_tmp_S0_pls_epsilon, TouchBarrier_S0_pls_epsilon);
			  //printf("                                      Innere tmpStep = %d  S_tmp_S0 = %3.5f, %3.5f = S_pls_eps, bool_S0 = %i, %i = bool_eps\n",tmp_step_noX,S_tmp_S0, S_tmp_S0_pls_epsilon, TouchBarrier_S0, TouchBarrier_S0_pls_epsilon);
			}

		}

		//printf("2.  S_tmp_S0 = %f, Boolian_S0 = %i\n", S_tmp_S0, TouchBarrier_S0);
		//printf("2.  S_tmp_ep = %f, Boolian_ep = %i\n", S_tmp_S0_pls_epsilon, TouchBarrier_S0_pls_epsilon);
		//printf("Diskontierungszeitraum = %f\n", (1.-(double)current_time/N_timeSteps)*maturity);


	  //Determine the correct price of the bonus certificate
		if(TouchBarrier_S0 || (S_tmp_S0 > Bonus_Level) )  tmp_price_S0 = S_tmp_S0*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
		else tmp_price_S0 = Bonus_Level*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);																												// If the value of the underlying has never touched or undershoot the Barrier and if the value of the underlying is smaller than the Bonuslevel than return the Bonuslevel = Price of Bonus certificate

		if(TouchBarrier_S0_pls_epsilon || (S_tmp_S0_pls_epsilon > Bonus_Level) )  tmp_price_S0_pls_epsilon = S_tmp_S0_pls_epsilon*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
		else tmp_price_S0_pls_epsilon = Bonus_Level*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);		

		tmp_sum_price_S0 += tmp_price_S0;
		tmp_sum_price_S0_pls_epsilon += tmp_price_S0_pls_epsilon;

	}


	price_S0 = tmp_sum_price_S0/N_inner_paths;
	price_S0_pls_epsilon = tmp_sum_price_S0_pls_epsilon/N_inner_paths;

	tmp_struct.priceDerivat_S0 = price_S0;
	tmp_struct.priceDerivat_S0_pls_epsilon = price_S0_pls_epsilon;


	return tmp_struct;

}






