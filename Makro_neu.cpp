
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

struct weights{
	double d1;
	double d2;
};

struct PriceBonusDerivat_pls_Delta{
	double priceDerivat;
	bool BarrierTouched;
};

struct PriceBonusDerivat_pls_Delta2{
	double priceDerivat_S0;
	double priceDerivat_S0_pls_epsilon;
	bool BarrierTouched_S0;
	bool BarrierTouched_S0_pls_epsilon;
};


double Get_StdNormDistr_RandomNumber();
double Evaluate_BonusDerivat_loop(double S_0, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps );
double Evaluate_BonusDerivat_Rekursiv(double S_0, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps );

weights Get_weight_d1_and_d2(const double S_0, const double K, const double maturity, const double Sigma, const double r_f, const double CoC_q);
double Get_Put(const double S_0, const double K, const double maturity, const double Sigma, const double r_f, const double CoC_q);
double Get_Call(const double S_0, const double K, const double maturity, const double Sigma, const double r_f, const double CoC_q);
double Get_Down_and_In_Put(const double S_0, const double K, const double Barrier, const double maturity, const double Sigma, const double r_f, const double CoC_q);
double Get_Down_and_Out_Put(const double S_0, const double K, const double Barrier, const double maturity, const double Sigma, const double r_f, const double CoC_q);
double Get_Theory_Value(double Put_down_out_value, double Call_value);
double Cumul_StNormFunc(const double x);
void Calculate_Delta(const double S_0, const double Sigma, const double Barrier, const double Bonus_Level, const double mu, const double r_f, const double epsilon, const double maturity, const int N_paths, const int N_timeSteps);
PriceBonusDerivat_pls_Delta Evaluate_BonusDerivat_loop_for_delta(const double S_0, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps, int current_time, bool Barrier_already_touched );
PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta2(const double S_0, const double S_0_pls_epsilon,
																																	const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity,
	 																																const int N_timeSteps, int current_time, bool Barrier_already_touched_S0, 
	 																																bool Barrier_already_touched_S0_pls_epsilon );
PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta21(const double S_0, const double S_0_pls_epsilon,
																																	 const double Sigma, const double r_f, const double Barrier, const double Bonus_Level,
																																	 const int maturity, const int N_timeSteps, const int N_inner_paths, int current_time,
																																	 bool Barrier_already_touched_S0, bool Barrier_already_touched_S0_pls_epsilon );


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
	int maturity       = 1;    		// Maturity of the derivative in [a]
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


	double time 				 ; // Point of time at which 
	double tmp_price     ; // Temporary price loop method 
	double tmp_sum       ; // Sum of the temporary prices recursive methodof all paths calculated
	double mean_price_nd ; // Final price loop method is eqaul to mean value of N_paths results. This one gets not discounted for calculating sigma value.
	double price       	 ; // Final price loop method is eqaul to mean value of N_paths results. This one gets discounted
	double tmp_sigma   	 ; // Sum used for calculating sigma: Sum((S_tmp-S_0)²)
	double sigma       	 ; // Sigma of the final result calculated as Sqrt(Sum((S_tmp-S_0)²)/N_paths)
 
	std::vector<int> Diff_Number_Of_Paths     = {10, 100, 1000, 10000, 100000};
	std::vector<int> Diff_Number_Of_timeSteps = {12, 52, 250, 3000, 180000};

  std::ofstream file_output;
  file_output.open("Prices_and_Sigmas.csv");

  for(int pathStyle =0; pathStyle <= Diff_Number_Of_Paths.size()-1; pathStyle++){
  	for(int timeStep_Style =0; timeStep_Style <= Diff_Number_Of_timeSteps.size()-1; timeStep_Style++){

     	N_paths     = Diff_Number_Of_Paths[pathStyle];
			N_timeSteps = Diff_Number_Of_timeSteps[timeStep_Style];

  		std::ofstream conf_out;
  		std::string out_file = "Prices_Zertifikate" + std::to_string(N_paths) + "_" + std::to_string(N_timeSteps) + ".csv";
  		conf_out.open(out_file);
         	
     	for (int i = 1; i <= 100; i++){

				tmp_price     			= 0.0;
				tmp_sum       			= 0.0;
				mean_price_nd		    = 0.0;
				price       				= 0.0; 
				tmp_sigma   				= 0.0; 
				sigma       				= 0.0; 

				for(int path = 1; path <= N_paths; path++){

					time = (double) path/N_paths;

					//____ IMPORTANT______//
					// You can not evaluate the bonus certificate
					// with the loop method and the recursive method
					// at the same time. Results are then inonsistent
					// because random numbers are different. 
					// Always comment one function call !!

					//Use for the Loop method
					tmp_price = Evaluate_BonusDerivat_loop(S_0, Sigma, r_f, Barrier, Bonus_Level, maturity, N_timeSteps );
					tmp_sum += tmp_price;

					//Use for the recursive method
					//tmp_price = Evaluate_BonusDerivat_Rekursiv(S_0, Sigma, r_f, Barrier, Bonus_Level, maturity, N_timeSteps );
					//tmp_sum += tmp_price;
					//tmp_step_noX_R =0;

				}//end loop over number of paths

				mean_price_nd = tmp_sum/N_paths;
				price 	      = tmp_sum/(N_paths )* exp(-r_f * maturity);

//				printf("#####____________________________________________________##########\n");
//				printf("#####_____Price bonus certificate: Loop = %3.10f\n", price);
//				printf("#####____________________________________________________##########\n");

				conf_out << price << "\n";

		  } //end loop to generate statistics for each configuration
		  conf_out.close();

      // get out mean values and sigma

			// read file line by line
			std::ifstream file(out_file);
			std::string line;
			double mean_value_price = 0.0;
			double sigma = 0.0;
			int number_of_prices = 0;

			while(getline(file, line)){
				number_of_prices++;
				mean_value_price += std::stof(line); 
			}
			std::cout << number_of_prices << std::endl;
			mean_value_price /= number_of_prices;
      file.clear();                 // clear fail and eof bits
      file.seekg(0, std::ios::beg); // back to the start!
      //calculate sigma
			while(getline(file, line)){
				sigma += pow((mean_value_price - std::stof(line)),2);
			}      
      sigma = sqrt(sigma/(number_of_prices*(number_of_prices-1)));  

      file_output << N_paths << "," << N_timeSteps << "," << mean_value_price << "," << sigma <<"\n";

			file.close();

		} // end loop through vector of different number of time steps -> e.g. this loop is just for taking various sizes of time intervalls
	} // end loop through vector of different number of paths -> e.g. this loop is just for taking different numbers of paths

	file_output.close();

/*
					// read file line by line
					std::ifstream file("tmp_Prices_Zertifikate" + std::to_string(N_paths) + "_" + std::to_string(N_timeSteps) + ".csv");
					std::string line;
					double tmp_value_file = 0.0;

					std::size_t a, b;
			
					while(getline(file, line)){
						a = line.find(',', 0);
						b = line.find('\"', a + 1);
						if (b != std::string::npos){
						std::string tmp_string = line.substr(a, 10);
						tmp_value_file = std::stof(tmp_string);
						}
						//a = line.find('', b + 1);
						//tmp_value_file = std::stof(line);
						tmp_sigma += pow((tmp_value_file - mean_price_nd),2);
					}
					file.close();

					sigma = sqrt(tmp_sigma/(N_paths*(N_paths-1)));
					//printf("Sigma deviation of price estimated via MC = %3.12f\n", sigma);

					file_output << N_paths << "," << N_timeSteps << "," << price << "," << sigma <<"\n";
*/


  	double Theory_call_value 						  = Get_Call(S_0, K, maturity, Sigma, r_f, CoC_q);
    double Theory_put_down_and_out_value  = Get_Down_and_Out_Put(S_0, Bonus_Level, Barrier, maturity, Sigma, r_f, CoC_q);
    double Theory_bonus_certificate_value = Get_Theory_Value(Theory_put_down_and_out_value, Theory_call_value);

    printf("Theory value for Bonus certificate = %3.10f\n", Theory_bonus_certificate_value);

//double mu = 0.1;
//double epsilon = 0.01;

  //Calculate_Delta( S_0,  Sigma,  Barrier,  Bonus_Level,  mu,  r_f,  epsilon,  maturity,  1,  1);



}//end main function


//_____________________________________________________________________________________________________________________________________________________________________________________________________
//

double Evaluate_BonusDerivat_loop(const double S_0, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps ){
	
	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation
	double S_tmp = 0;
	double S_tmp_initial = S_0;
	double random_number;
	bool TouchBarrier = false;// Boolian for checking whether barrier was touched or undershoot

	for(int tmp_step_noX = 1; tmp_step_noX <= N_timeSteps; tmp_step_noX++){
		random_number = Get_StdNormDistr_RandomNumber(); // Get the standard normal distributed random number
		if(tmp_step_noX < 5)		//printf("Random number Evaluate BonusDerivatLoop = %f\n", random_number);                                                                              
		S_tmp = S_tmp_initial * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
		S_tmp_initial = S_tmp;

		//Check if barrier was touched
		if(S_tmp <= Barrier) TouchBarrier = true;                                        // If barrier was touched set the boolian to kTRUE 
    //printf("Wert St_Loop = %3.5f and number tmpStep = %d, timeStep*delta_t = %3.3f\n",S_tmp,  tmp_step_noX, tmp_step_noX* delta_t);
	}

  //Determine the correct price of the bonus certificate
	if(TouchBarrier || (S_tmp > Bonus_Level) )  return S_tmp;                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
	else return Bonus_Level;																												// If the value of the underlying has never touched or undershoot the Barrier and if the value of the underlying is smaller than the Bonuslevel than return the Bonuslevel = Price of Bonus certificate
}

//_____________________________________________________________________________________________________________________________________________________________________________________________________
//


//Beware of the stack overflow for too deep recursions!!!
double Evaluate_BonusDerivat_Rekursiv(double S_0, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps ){
	
	bool TouchBarrier_R = false;// Boolian for checking whether barrier was touched or undershoot
	double delta_t = ((double) maturity)/N_timeSteps; // time intervall for each evaluation
  tmp_step_noX_R += 1;

	double random_number = Get_StdNormDistr_RandomNumber();                                                                               // Get the standard normal distributed random number                                                                               // Get the standard normal distributed random number
	double S_tmp = S_0 * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX

	//Check if barrier was touched
	if(S_tmp <= Barrier) TouchBarrier_R = true;                                        // If barrier was touched set the boolian to kTRUE 

	//Check if the end of maturity is reached and determine the correct price of the bonus certificate
	if(tmp_step_noX_R < N_timeSteps){																						// Check if the end of maturity is reached; If not than repeat the procedure; if yes than return the appropriate value
	  //printf("Wert St_Rekursiv = and number tmpStep = %d, timeStep*delta_t = %3.3f\n",  tmp_step_noX_R, tmp_step_noX_R* delta_t);
		return Evaluate_BonusDerivat_Rekursiv(S_tmp, Sigma, r_f, Barrier, Bonus_Level, maturity, N_timeSteps);
	
	}
	else{
		if(TouchBarrier_R || (S_tmp > Bonus_Level) )  return S_tmp;                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
		else return Bonus_Level;																												// If the value o the underlying has never touched or undershoot the Barrier and if the value of the underlying is smaller than the Bonuslevel than return the Bonuslevel = Price of Bonus certificate
	}

}

//_____________________________________________________________________________________________________________________________________________________________________________________________________
//

double Get_StdNormDistr_RandomNumber(){

  std::normal_distribution<> nd{0,1}; // Define normal distribution with name nd and mean = 0 and sigma =1 which is equal to standard normal distribution; Return type is double
	double random_number = nd(gen);           // Generate random number which is standard normal distributed

	return random_number;
}

//_____________________________________________________________________________________________________________________________________________________________________________________________________

double Get_Theory_Value(double Put_down_out_value, double Call_value){

	double theory_value = Put_down_out_value + Call_value;
	return theory_value;
}

double Get_Down_and_Out_Put(const double S_0, const double K, const double Barrier, const double maturity, const double Sigma, const double r_f, const double CoC_q){

	double P_DI = Get_Down_and_In_Put(S_0, K, Barrier, maturity, Sigma, r_f, CoC_q);
	double put  = Get_Put(S_0, K, maturity, Sigma, r_f, CoC_q);
	double P_DO = put - P_DI;
	return P_DO;

}

//_____________________________________________________________________________________________________________________________________________________________________________________________________
double Get_Down_and_In_Put(const double S_0, const double K, const double Barrier, const double maturity, const double Sigma, const double r_f, const double CoC_q){

	double lambda = (r_f - CoC_q + (Sigma*Sigma)/2.)/(Sigma*Sigma);
	double x_1    = log(S_0/Barrier) / (Sigma * sqrt(maturity)) + lambda * Sigma * sqrt(maturity);
	double y      = log(Barrier*Barrier / (S_0 * K)) / (Sigma*sqrt(maturity)) + lambda * Sigma * sqrt(maturity);
	double y_1    = log(Barrier/S_0) / (Sigma*sqrt(maturity)) + lambda * Sigma * sqrt(maturity);

	struct weights tmp_weight = Get_weight_d1_and_d2(S_0, K, maturity, Sigma, r_f, CoC_q);
	double d1 = tmp_weight.d1;
	double d2 = tmp_weight.d2;

	double value_down_and_in_put = -S_0 * Cumul_StNormFunc(-x_1) * exp(-CoC_q * maturity) 
																 + K * exp(-r_f * maturity) * Cumul_StNormFunc(-x_1 + Sigma * sqrt(maturity))
																 + S_0 * exp(-CoC_q * maturity) * pow(Barrier/S_0, 2*lambda) * (Cumul_StNormFunc(y) - Cumul_StNormFunc(y_1))
																 - K * exp(-r_f * maturity) * pow(Barrier/S_0, 2*lambda -2) * (Cumul_StNormFunc(y - Sigma*sqrt(maturity)) - Cumul_StNormFunc(y_1 - Sigma*sqrt(maturity)));
	return value_down_and_in_put;
}
//_____________________________________________________________________________________________________________________________________________________________________________________________________

double Get_Call(const double S_0, const double K, const double maturity, const double Sigma, const double r_f, const double CoC_q){

	struct weights tmp_weight = Get_weight_d1_and_d2(S_0, K, maturity, Sigma, r_f, CoC_q);
	double d1 = tmp_weight.d1;
	double d2 = tmp_weight.d2;
	double call_value = -100.0;

	if(K == 0.){
		call_value = S_0 * exp(-CoC_q * maturity);
	}
	else{
		call_value = S_0 * exp(-CoC_q * maturity) * Cumul_StNormFunc(d1) - K * exp(-r_f * maturity) * Cumul_StNormFunc(d2); 
  }
	return call_value;
}
//_____________________________________________________________________________________________________________________________________________________________________________________________________

double Get_Put(const double S_0, const double K, const double maturity, const double Sigma, const double r_f, const double CoC_q){

	struct weights tmp_weight = Get_weight_d1_and_d2(S_0, K, maturity, Sigma, r_f, CoC_q);
	double d1 = tmp_weight.d1;
	double d2 = tmp_weight.d2;
	double put_value = -1000.0;

	if(K == 0.){
		put_value = 0.0;
	}
	else{
		put_value = K * exp(-r_f * maturity) * Cumul_StNormFunc(-d2) - S_0 * exp(-CoC_q * maturity) * Cumul_StNormFunc(-d1);
	}

	return put_value;
}
//_____________________________________________________________________________________________________________________________________________________________________________________________________

weights Get_weight_d1_and_d2(const double S_0, const double K, const double maturity, const double Sigma, const double r_f, const double CoC_q){
	// K = strike price

	struct weights tmp_struct;
	tmp_struct.d1 = (log(S_0/K) + (r_f - CoC_q + Sigma*Sigma/2.)/maturity)/(Sigma/sqrt(maturity));
	tmp_struct.d2 = tmp_struct.d1 - Sigma * sqrt(maturity);

	return tmp_struct;
}
//_____________________________________________________________________________________________________________________________________________________________________________________________________

double Cumul_StNormFunc(const double d_x){        // Phi(-∞, x) aka N(x)
    return std::erfc(-d_x/std::sqrt(2))/2;				// https://en.cppreference.com/w/cpp/numeric/math/erfc --> https://en.wikipedia.org/wiki/Error_function#Complementary_error_function
}

void Calculate_Delta(const double S_0, const double Sigma, const double Barrier, const double Bonus_Level, const double mu, const double r_f, const double epsilon, const double maturity, const int N_paths, const int N_timeSteps){

	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation

	std::ofstream file_output_StockPrice;
	std::ofstream file_mean_outputStockPrice;
	// std::ofstream file_output_CertificatePrice;
	// std::ofstream file_output_Delta;
	// std::ofstream file_output_RiskFreeAsset;

	file_output_StockPrice.open("StockPrice.csv");


	for(int path = 1; path <= N_paths; path++){

	double S_tmp_rf = 0.0;
	double S_tmp_initial_rf = S_0;
	double S_tmp_rf_pls_epsilon = 0.0;

	double S_tmp_mu = 0.0;
	double S_tmp_mu_pls_epsilon = 0.0;
	double S_tmp_initial_mu = S_0;

	double D_tmp = 0.0;
	double D_tmp_pls_epsilon = 0.0;

	double delta = 0.0;
	double RiskfreeAsset = 0.0;

	double replication_portfolio_vorher = 0.0;
	double replication_portfolio_nachher = 0.0;
	double earnigs_and_losses = 0.0;

	struct PriceBonusDerivat_pls_Delta Struct_price_D_tmp;
	Struct_price_D_tmp.priceDerivat = 0.0;
	Struct_price_D_tmp.BarrierTouched = false; // Initialising for every path with "false". 
	//That means that for tmp_step_noX = 1 the boolian is false. Then after tmp_step_noX = 1 the boolian gets updated for each time step iteration.
	//Therefore, it is most likely that for "early" time steps it is changed to true in the function "Evaluate_BonusDerivat_loop_for_delta" since in the fuction
	//there are N_timeStep - tmp_step_noX iterations in which it can be changed to one. 

	struct PriceBonusDerivat_pls_Delta Struct_price_D_tmp_pls_epsilon;
	Struct_price_D_tmp_pls_epsilon.priceDerivat = 0.0;
	Struct_price_D_tmp_pls_epsilon.BarrierTouched = false;

	struct PriceBonusDerivat_pls_Delta2 Struct_2;
	Struct_2.priceDerivat_S0 = 0.0;
	Struct_2.BarrierTouched_S0 = false;
	Struct_2.priceDerivat_S0_pls_epsilon = 0.0;
	Struct_2.BarrierTouched_S0_pls_epsilon = false;


	bool TouchBarrier_D             = false;// Boolian for checking whether barrier was touched or undershoot by the value of the certificate
	bool TouchBarrier_D_pls_epsilon = false;// Boolian for checking whether barrier was touched or undershoot by the value of the certificate plus epsilon

double kosten_aktien_Kauf_Verkauf = 0.0;
double Sum_kosten_aktien_Kauf_Verkauf = 0.0;
double zinskostenAktien = 0.0;
double Sum_zinskostenAktien = 0.0;
double Kosten_Hedge = 0.0;

double delta_tmp = 0.0;

	double random_number;

	printf("___Path: %d   _____________________\n", path);

		for(int tmp_step_noX = 1; tmp_step_noX <= N_timeSteps; tmp_step_noX++){
			printf("        ÄussererTimeStep:  %d   _____________________\n", tmp_step_noX);

		  random_number = Get_StdNormDistr_RandomNumber();                                                            // Get the standard normal distributed random number

		  // Calculate the Stock price using mu
			S_tmp_mu = S_tmp_initial_mu * exp( (mu - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the stock at the n-st time step: delta_t* tmp_step_noX using mu
		
			if(tmp_step_noX < 2 ) S_tmp_mu = S_tmp_initial_mu;
			else S_tmp_initial_mu = S_tmp_mu;

			S_tmp_mu_pls_epsilon = S_tmp_mu + epsilon;
			//Struct_price_D_tmp_pls_epsilon = Evaluate_BonusDerivat_loop_for_delta(S_tmp_mu_pls_epsilon, Sigma, r_f, Barrier, Bonus_Level, maturity, N_timeSteps,  tmp_step_noX , TouchBarrier_D_pls_epsilon);
			//Struct_price_D_tmp 						 = Evaluate_BonusDerivat_loop_for_delta(S_tmp_mu, Sigma, r_f, Barrier, Bonus_Level, maturity, N_timeSteps,  tmp_step_noX , TouchBarrier_D);

			Struct_2 = Evaluate_BonusDerivat_loop_for_delta21(S_tmp_mu, S_tmp_mu_pls_epsilon, Sigma, r_f, Barrier, Bonus_Level, maturity, 1, 1, tmp_step_noX ,TouchBarrier_D, TouchBarrier_D_pls_epsilon);


			D_tmp 						         = Struct_2.priceDerivat_S0;
			D_tmp_pls_epsilon          = Struct_2.priceDerivat_S0_pls_epsilon;
			TouchBarrier_D 						 = Struct_2.BarrierTouched_S0;
			TouchBarrier_D_pls_epsilon = Struct_2.BarrierTouched_S0_pls_epsilon;

			/*
			D_tmp 						= Struct_price_D_tmp.priceDerivat;
			D_tmp_pls_epsilon = Struct_price_D_tmp_pls_epsilon.priceDerivat;
			TouchBarrier_D 						 = Struct_price_D_tmp.BarrierTouched;
			TouchBarrier_D_pls_epsilon = Struct_price_D_tmp_pls_epsilon.BarrierTouched;
			*/

			delta = (D_tmp_pls_epsilon - D_tmp )/epsilon;
			replication_portfolio_vorher = (delta_tmp * S_tmp_mu)-D_tmp;

			kosten_aktien_Kauf_Verkauf = (delta_tmp -delta)*S_tmp_mu;// in t_0 delta_tmp = 0
			Sum_kosten_aktien_Kauf_Verkauf += kosten_aktien_Kauf_Verkauf;
			//zinskostenAktien = - delta *S_tmp_mu* exp(r_f*delta_t);// Bullshit mit stetiger verzinsung ?

			Sum_zinskostenAktien += zinskostenAktien * (exp(-r_f*delta_t)-1.); // Zinseszinseffekt auf Zinskonto ??
			//Jetzt erst neue Zinsen aufs konto dazuaddieren. In einem jahr wird konto wieder verzinst...
			zinskostenAktien = - delta *S_tmp_mu*(exp(-r_f*delta_t)-1.);// Bullshit mit stetiger verzinsung ?
			Sum_zinskostenAktien += zinskostenAktien;

			delta_tmp = delta;
			replication_portfolio_nachher = (delta * S_tmp_mu)-D_tmp;

			
/*
			if(tmp_step_noX < 2 ){
				delta = (D_tmp_pls_epsilon - D_tmp )/epsilon;
				// delta = (S_tmp_rf_pls_epsilon - S_tmp_rf )/epsilon;
				RiskfreeAsset = (D_tmp - delta*S_tmp_mu);
				replication_portfolio = delta * S_tmp_mu + RiskfreeAsset;
				printf("Path:%d,TimeStep:%d:PrintfNo:__7__delta = %f, RiskfreeAsset = %f, replication_portfolio = %f\n",path,tmp_step_noX, delta, RiskfreeAsset, replication_portfolio );
			}
			else{
				//printf("Path:%d,TimeStep:%d:PrintfNo:__8__delta = %f, RiskfreeAsset = %f, replication_portfolio = %f\n",path,tmp_step_noX, delta, RiskfreeAsset, replication_portfolio );
				//Delta and RiskfreeAsset were not initialized, yet !!
				replication_portfolio = delta * S_tmp_mu + RiskfreeAsset * exp(r_f*delta_t);
	      earnigs_and_losses =+ (D_tmp - replication_portfolio);

	      //No initializing delta and the risk free asset
	      delta = (D_tmp_pls_epsilon - D_tmp)/epsilon;
	      // delta = (S_tmp_rf_pls_epsilon - S_tmp_rf )/epsilon;
				//if(delta > 1 || delta < -1) printf("##_____Error in Function Calculate_Delta !!! Delta larger than 1 or smaller than -1 ! Delta = %3.3f____##\n", delta);
				RiskfreeAsset = (D_tmp - delta*S_tmp_mu)*exp(-r_f*delta_t);

				printf("Path:%d,TimeStep:%d: delta = %f, RiskfreeAsset = %f, repl_portf = %f\n",path,tmp_step_noX, delta, RiskfreeAsset, replication_portfolio );
				//printf("######____delta = %f  ____________######\n",delta);
		  }
*/

				printf("Path:%d,TimeStep:%d: Aktkurs = %f, Derivat = %f, delta = %f, repli_portf_vor = %f, repli_portf_nach = %f, kosten_akti = %f, zinskostAkt = %f \n",path,tmp_step_noX, S_tmp_mu, D_tmp, delta, replication_portfolio_vorher, replication_portfolio_nachher, kosten_aktien_Kauf_Verkauf, zinskostenAktien );
				printf("Sum_kosten_aktien_Kauf_Verkauf = %f, Sum_zinskostenAktien =%f\n", Sum_kosten_aktien_Kauf_Verkauf, Sum_zinskostenAktien);

		  file_output_StockPrice << delta << ",";
		}// end loop through time steps
		file_output_StockPrice << "\n";

		Kosten_Hedge = (Sum_zinskostenAktien + Sum_kosten_aktien_Kauf_Verkauf + delta* S_tmp_mu)*(exp(-r_f*1));// discount for one year since maturity is one year

		printf("Kosten des Hedge = %f, Delta = %f, S_tmp_mu = %f, Sum_zinskostenAktien = %f, Sum_kosten_aktien_Kauf_Verkauf = %f\n", Kosten_Hedge, delta, S_tmp_mu, Sum_zinskostenAktien, Sum_kosten_aktien_Kauf_Verkauf);
	}//end loop through paths
	file_output_StockPrice.close();

			 // read file line by line
			std::ifstream file("StockPrice.csv");

			std::string line;
			double tmp_value_file = 0.0;
			double sum = 0.0;
			double mean = 0.0;
			file_mean_outputStockPrice.open("file_mean_outputStockPrice.csv");

			while(file){
				std::string s;
				getline(file, s );
		    std::istringstream ss( s );

		    while (ss){
		      std::string s;
		      if (!getline( ss, s, ',' )) break;
		      tmp_value_file = std::stof(s);
		      sum += tmp_value_file;
		    }
				mean = sum/N_paths;
				file_mean_outputStockPrice << mean <<"\n";
			}
			file_mean_outputStockPrice.close();

}

//_____________________________________________________________________________________________________________________________________________________________________________________________________
//

PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta21(const double S_0, const double S_0_pls_epsilon, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps, const int N_inner_paths, int current_time, bool Barrier_already_touched_S0, bool Barrier_already_touched_S0_pls_epsilon ){
	
	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation
	double price_S0 = 0.0;
	double price_S0_pls_epsilon = 0.0;
	double tmp_sum_price_S0 = 0.0;
	double tmp_sum_price_S0_pls_epsilon = 0.0;
	bool TouchBarrier_S0 = Barrier_already_touched_S0;// Boolian for checking whether barrier was touched or undershoot
	bool TouchBarrier_S0_pls_epsilon = Barrier_already_touched_S0_pls_epsilon;// Boolian for checking whether barrier was touched or undershoot
	struct PriceBonusDerivat_pls_Delta2 tmp_struct;

	//printf("1.  S0_init    = %f, current_time = %i, boolian_S0 = %i\n", S_tmp_initial_S0, current_time,TouchBarrier_S0);
	//printf("1.  S0_init_ep = %f, current_time = %i, boolian_ep = %i\n", S_tmp_initial_S0_pls_epsilon, current_time,TouchBarrier_S0_pls_epsilon);
	
	for(int path = 1; path <= N_inner_paths; path++){

		//printf("_______________________Inner Path:  %d_____________\n",path);

		double S_tmp_S0 = 0.;
		double S_tmp_S0_pls_epsilon = 0.;
		double S_tmp_initial_S0 = S_0;
		double S_tmp_initial_S0_pls_epsilon = S_0_pls_epsilon;
		double random_number;
		//TouchBarrier_S0 = Barrier_already_touched_S0;// Boolian for checking whether barrier was touched or undershoot
		//TouchBarrier_S0_pls_epsilon = Barrier_already_touched_S0_pls_epsilon;// Boolian for checking whether barrier was touched or undershoot

		//bool TouchBarrier_S0 = false;// Boolian for checking whether barrier was touched or undershoot
		//bool TouchBarrier_S0_pls_epsilon = false;// Boolian for checking whether barrier was touched or undershoot

		double tmp_price_S0 = 0.0;
		double tmp_price_S0_pls_epsilon = 0.0;

		for(int tmp_step_noX = current_time; tmp_step_noX <= N_timeSteps; tmp_step_noX++){
			random_number = Get_StdNormDistr_RandomNumber();    
			//if(tmp_step_noX < 5)		printf("Random number Evaluate DELTA = %f\n", random_number);                                                                              // Get the standard normal distributed random number
			S_tmp_S0             = S_tmp_initial_S0 * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
			S_tmp_S0_pls_epsilon = S_tmp_initial_S0_pls_epsilon * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
			S_tmp_initial_S0 = S_tmp_S0;
			S_tmp_initial_S0_pls_epsilon = S_tmp_S0_pls_epsilon;

			//Check if barrier was touched
			if(S_tmp_S0 <= Barrier)             TouchBarrier_S0 = true;                                        // If barrier was touched set the boolian to kTRUE 
			if(S_tmp_S0_pls_epsilon <= Barrier) TouchBarrier_S0_pls_epsilon = true;                                        // If barrier was touched set the boolian to kTRUE 
//printf("                                      Innere tmpStep = %d  S_tmp_S0 = %3.5f, %3.5f = S_pls_eps, bool_S0 = %i, %i = bool_eps\n",tmp_step_noX,S_tmp_S0, S_tmp_S0_pls_epsilon, TouchBarrier_S0, TouchBarrier_S0_pls_epsilon);
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
	tmp_struct.BarrierTouched_S0 = TouchBarrier_S0;

	tmp_struct.priceDerivat_S0_pls_epsilon = price_S0_pls_epsilon;
	tmp_struct.BarrierTouched_S0_pls_epsilon = TouchBarrier_S0_pls_epsilon;

	//printf("3.  Return D_S0 = %f, Boolian_S0 = %i\n", tmp_struct.priceDerivat_S0, tmp_struct.BarrierTouched_S0);
	//printf("3.  Return D_ep = %f, Boolian_ep = %i\n", tmp_struct.priceDerivat_S0_pls_epsilon, tmp_struct.BarrierTouched_S0_pls_epsilon );

	return tmp_struct;

}



// In this version the boolians get updated and the random number is the same one for the calculation of the derivates using S0 or S0 plus epsilon

/*
PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta2(const double S_0, const double S_0_pls_epsilon, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps, int current_time, bool Barrier_already_touched_S0, bool Barrier_already_touched_S0_pls_epsilon ){
	
	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation
	double S_tmp_S0 = 0.;
	double S_tmp_S0_pls_epsilon = 0.;
	double S_tmp_initial_S0 = S_0;
	double S_tmp_initial_S0_pls_epsilon = S_0_pls_epsilon;
	double random_number;
	bool TouchBarrier_S0 = Barrier_already_touched_S0;// Boolian for checking whether barrier was touched or undershoot
	bool TouchBarrier_S0_pls_epsilon = Barrier_already_touched_S0_pls_epsilon;// Boolian for checking whether barrier was touched or undershoot

	//bool TouchBarrier_S0 = false;// Boolian for checking whether barrier was touched or undershoot
	//bool TouchBarrier_S0_pls_epsilon = false;// Boolian for checking whether barrier was touched or undershoot


	double price_S0 = 0.0;
	double price_S0_pls_epsilon = 0.0;
	struct PriceBonusDerivat_pls_Delta2 tmp_struct;

	//printf("1.  S0_init    = %f, current_time = %i, boolian_S0 = %i\n", S_tmp_initial_S0, current_time,TouchBarrier_S0);
	//printf("1.  S0_init_ep = %f, current_time = %i, boolian_ep = %i\n", S_tmp_initial_S0_pls_epsilon, current_time,TouchBarrier_S0_pls_epsilon);

	for(int tmp_step_noX = current_time; tmp_step_noX <= N_timeSteps; tmp_step_noX++){
		random_number = Get_StdNormDistr_RandomNumber();                                                                               // Get the standard normal distributed random number
		S_tmp_S0             = S_tmp_initial_S0 * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
		S_tmp_S0_pls_epsilon = S_tmp_initial_S0_pls_epsilon * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
		S_tmp_initial_S0 = S_tmp_S0;
		S_tmp_initial_S0_pls_epsilon = S_tmp_S0_pls_epsilon;

		//Check if barrier was touched
		if(S_tmp_S0 <= Barrier)             TouchBarrier_S0 = true;                                        // If barrier was touched set the boolian to kTRUE 
		if(S_tmp_S0_pls_epsilon <= Barrier) TouchBarrier_S0_pls_epsilon = true;                                        // If barrier was touched set the boolian to kTRUE 
    printf("                Innere tmpStep = %d  S_tmp_S0 = %3.5f, %3.5f = S_pls_eps, bool_S0 = %i, %i = bool_eps\n",tmp_step_noX,S_tmp_S0, S_tmp_S0_pls_epsilon, TouchBarrier_S0, TouchBarrier_S0_pls_epsilon);
	}

	//printf("2.  S_tmp_S0 = %f, Boolian_S0 = %i\n", S_tmp_S0, TouchBarrier_S0);
	//printf("2.  S_tmp_ep = %f, Boolian_ep = %i\n", S_tmp_S0_pls_epsilon, TouchBarrier_S0_pls_epsilon);
	//printf("Diskontierungszeitraum = %f\n", (1.-(double)current_time/N_timeSteps)*maturity);


  //Determine the correct price of the bonus certificate
	if(TouchBarrier_S0 || (S_tmp_S0 > Bonus_Level) )  price_S0 = S_tmp_S0*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
	else price_S0 = Bonus_Level*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);																												// If the value of the underlying has never touched or undershoot the Barrier and if the value of the underlying is smaller than the Bonuslevel than return the Bonuslevel = Price of Bonus certificate

	if(TouchBarrier_S0_pls_epsilon || (S_tmp_S0_pls_epsilon > Bonus_Level) )  price_S0_pls_epsilon = S_tmp_S0_pls_epsilon*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
	else price_S0_pls_epsilon = Bonus_Level*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);							


	tmp_struct.priceDerivat_S0 = price_S0;
	tmp_struct.BarrierTouched_S0 = TouchBarrier_S0;

	tmp_struct.priceDerivat_S0_pls_epsilon = price_S0_pls_epsilon;
	tmp_struct.BarrierTouched_S0_pls_epsilon = TouchBarrier_S0_pls_epsilon;

	//printf("3.  Return D_S0 = %f, Boolian_S0 = %i\n", tmp_struct.priceDerivat_S0, tmp_struct.BarrierTouched_S0);
	//printf("3.  Return D_ep = %f, Boolian_ep = %i\n", tmp_struct.priceDerivat_S0_pls_epsilon, tmp_struct.BarrierTouched_S0_pls_epsilon );

	return tmp_struct;

}
*/



// In this function the boolians are not getting updated

/*
PriceBonusDerivat_pls_Delta Evaluate_BonusDerivat_loop_for_delta(const double S_0, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps, int current_time, bool Barrier_already_touched ){
	
	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation
	double S_tmp = 0;
	double S_tmp_initial = S_0;
	double random_number;
	bool TouchBarrier = Barrier_already_touched;// Boolian for checking whether barrier was touched or undershoot
	double price = 0.0;
	struct PriceBonusDerivat_pls_Delta tmp_struct;

	printf("Evaluate_Bon: S_0 = %f, current_time = %i, boolian = %i\n", S_tmp_initial, current_time,TouchBarrier);

	for(int tmp_step_noX = current_time; tmp_step_noX <= N_timeSteps; tmp_step_noX++){
		random_number = Get_StdNormDistr_RandomNumber();                                                                               // Get the standard normal distributed random number
		S_tmp = S_tmp_initial * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
		S_tmp_initial = S_tmp;

		//Check if barrier was touched
		if(S_tmp <= Barrier) TouchBarrier = true;                                        // If barrier was touched set the boolian to kTRUE 
    //printf("______Wert St_Loop = %3.5f and number tmpStep = %d and bool = %i\n",S_tmp,  tmp_step_noX, TouchBarrier);
	}

	printf("2. Evaluate_BonusDerivat_loop_for_delta: S_tmp = %f, Boolian = %i\n", S_tmp, TouchBarrier);


  //Determine the correct price of the bonus certificate
	if(TouchBarrier || (S_tmp > Bonus_Level) )  price = S_tmp*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
	else price = Bonus_Level*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);																												// If the value of the underlying has never touched or undershoot the Barrier and if the value of the underlying is smaller than the Bonuslevel than return the Bonuslevel = Price of Bonus certificate

	tmp_struct.priceDerivat = price;
	tmp_struct.BarrierTouched = TouchBarrier;

	return tmp_struct;

}
*/