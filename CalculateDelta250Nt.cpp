
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
	double priceDerivat_S0_min_epsilon;
};


double Get_StdNormDistr_RandomNumber();

void Calculate_Delta(const double S_0, const double Sigma, const double Barrier, const double Bonus_Level, const double mu, const double r_f, const double epsilon, const double maturity, const int N_paths, const int N_timeSteps);
PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta21(const double S_0, const double S_0_pls_epsilon, const double S_0_min_epsilon, const double Sigma, 
																																		const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps,
																																		const int N_inner_paths, int current_time, bool Barrier_already_touched_S0,
																																		bool Barrier_already_touched_S0_pls_epsilon, bool Barrier_already_touched_S0_min_epsilon,
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

	double mu = 0.0855;
	double epsilon = 0.05;

  Calculate_Delta( S_0,  Sigma,  Barrier,  Bonus_Level,  mu,  r_f,  epsilon,  maturity,  100, 250 );



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

	std::ofstream file_Kosten_Delta_Hedging;
	file_Kosten_Delta_Hedging.open("Kosten_Delta_Hedging.csv");
	file_Kosten_Delta_Hedging << "Outer_PathNo" << "," << "FV_Zinskosten" << "," << "FV_AktienKaufVerkauf" << "," << "FV_delta_S_tmp_muT" << "," << "PV_KostenHedge" << "\n";

	for(int path = 1; path <= N_paths; path++){

	std::ofstream file_output_Delta;
	file_output_Delta.open("Delta_" + std::to_string(path) + ".csv");
	file_output_Delta << "time" << "," << "S_tmp_mu" << "," << "D_tmp" << "," << "D_tmp_pls_epsilon" << "," << "D_tmp_min_epsilon" << "," << "delta" << "," << "gamma" << 
	"," << "gamma1" << "," << "kosten_aktien_Kauf_Verkauf" << "," << "zinskostenAktien" << "," << "RiskFreeAsset_vorUmschichtung" << "," 
	<< "RiskFreeAsset_NachUmschichtung" << "," << "Gap_value_vorUmsch" << "," << "Gap_valueNachUmsch" << "," << "Sum_kosten_aktien_Kauf_Verkauf" << "," << "Sum_zinskostenAktien" << "\n";

	double S_tmp_mu = 0.0;
	double S_tmp_mu_alt = 0.0;
	double S_tmp_mu_pls_epsilon = 0.0;
	double S_tmp_mu_min_epsilon = 0.0;
	double S_tmp_initial_mu = S_0;

	double D_tmp = 0.0;
	double D_tmp_alt = 0.0;
	double D_tmp_pls_epsilon = 0.0;
	double D_tmp_min_epsilon = 0.0;

	double delta = 0.0;
	double gamma = 0.0;
	double gamma1 = 0.0;
	double RiskfreeAsset = 0.0;

	double RiskFreeAsset_beforeRebalancing = 0.0;
	double RiskFreeAsset_afterRebalancing = 0.0;
	double earnigs_and_losses = 0.0;
// Initialising for every path with "false". 
	//That means that for tmp_step_noX = 1 the boolian is false. Then after tmp_step_noX = 1 the boolian gets updated for each time step iteration.
	//Therefore, it is most likely that for "early" time steps it is changed to true in the function "Evaluate_BonusDerivat_loop_for_delta" since in the fuction
	//there are N_timeStep - tmp_step_noX iterations in which it can be changed to one. 


	struct PriceBonusDerivat_pls_Delta2 Struct_2;
	Struct_2.priceDerivat_S0 = 0.0;
	Struct_2.priceDerivat_S0_pls_epsilon = 0.0;
	Struct_2.priceDerivat_S0_min_epsilon = 0.0;


	bool TouchBarrier_D1             = false;// Boolian for checking whether barrier was touched or undershoot by the value of the certificate
	bool TouchBarrier_D1_pls_epsilon = false;// Boolian for checking whether barrier was touched or undershoot by the value of the certificate plus epsilon
	bool TouchBarrier_D1_min_epsilon = false;// Boolian for checking whether barrier was touched or undershoot by the value of the certificate min epsilon
	//bool DividendsPaid = false; // use when dividends are supposed to pay
	bool DividendsPaid = true; // use if dividends should not be used
	double Dividend = 3.30;

	double kosten_aktien_Kauf_Verkauf = 0.0;
	double Sum_kosten_aktien_Kauf_Verkauf = 0.0;
	double zinskostenAktien = 0.0;
	double Sum_zinskostenAktien = 0.0;
	double Kosten_Hedge = 0.0;

	double delta_old = 0.0;
	double gap_Value = 0.0; // Das ist differenz zwischen Hedging portfolio vor der Umschichtung und dem aktuellen Derivate Wert 
													// = (delta_old * S_aktuell + riskfreeAsset_vorher) -D_aktuell
	double check_gap_value_after_rebalancing = 0.0;

	double random_number;
	double time = 0.0;

	printf("___Path: %d   _____________________\n", path);

		for(int tmp_step_noX = 1; tmp_step_noX <= N_timeSteps; tmp_step_noX++){
			time = (double)tmp_step_noX/N_timeSteps;
			//printf("Path:%d,TimeStep:%d: S_tmp_mu = %f %f = S_tmp_mu_pls_epsilon \n",path,tmp_step_noX, S_tmp_mu, S_tmp_mu_pls_epsilon);
			//printf("        ÄussererTimeStep:  %d   _____BoolS0 = %i, BoolS0Ep = %i________________\n", tmp_step_noX,TouchBarrier_D1, TouchBarrier_D1_pls_epsilon);
      //printf("Timestep\n");
		  random_number = Get_StdNormDistr_RandomNumber();                                                            // Get the standard normal distributed random number
		  
		if(tmp_step_noX > N_timeSteps/2. && !DividendsPaid){
				DividendsPaid = true;
				S_tmp_initial_mu -= Dividend;
			}

		  // Calculate the Stock price using mu
			//S_tmp_mu = S_tmp_initial_mu * exp( (mu - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the stock at the n-st time step: delta_t* tmp_step_noX using mu
			S_tmp_mu = S_tmp_initial_mu * exp(mu * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the stock at the n-st time step: delta_t* tmp_step_noX using mu
		
			if(tmp_step_noX < 2 ) S_tmp_mu = S_tmp_initial_mu;
			else S_tmp_initial_mu = S_tmp_mu;

			S_tmp_mu_pls_epsilon = S_tmp_mu + epsilon;
			S_tmp_mu_min_epsilon = S_tmp_mu - epsilon;

	

		  //Check if barrier was touched
			if(S_tmp_mu <= Barrier){             
				TouchBarrier_D1 = true;                                        // If barrier was touched set the boolian to kTRUE 
				//printf("Äussere Barrier S_tmp_mu touched ! Path = %d, Time step = %d, S_tmp_mu = %f, Bool = %i\n", path,tmp_step_noX, S_tmp_mu, TouchBarrier_D1);
				}
			if(S_tmp_mu_pls_epsilon <= Barrier){
				TouchBarrier_D1_pls_epsilon = true;  
			  //printf("____Äussere Barrier S_tmp_mu_pls_epsilon touched ! Path = %d, Time step = %d, S_tmp_mu_pls_epsilon = %f, Bool = %i\n", path,tmp_step_noX, S_tmp_mu_pls_epsilon, TouchBarrier_D1_pls_epsilon);
			}
			if(S_tmp_mu_min_epsilon <= Barrier){
				TouchBarrier_D1_min_epsilon = true;  
			  //printf("____Äussere Barrier S_tmp_mu_pls_epsilon touched ! Path = %d, Time step = %d, S_tmp_mu_pls_epsilon = %f, Bool = %i\n", path,tmp_step_noX, S_tmp_mu_pls_epsilon, TouchBarrier_D1_pls_epsilon);
			}


			Struct_2 = Evaluate_BonusDerivat_loop_for_delta21(S_tmp_mu, S_tmp_mu_pls_epsilon,S_tmp_mu_min_epsilon, Sigma, r_f, Barrier, Bonus_Level, maturity, N_timeSteps, 1000000, tmp_step_noX,TouchBarrier_D1, TouchBarrier_D1_pls_epsilon, TouchBarrier_D1_min_epsilon, DividendsPaid, Dividend);
			
			D_tmp 						         = Struct_2.priceDerivat_S0;
			D_tmp_pls_epsilon          = Struct_2.priceDerivat_S0_pls_epsilon;
			D_tmp_min_epsilon          = Struct_2.priceDerivat_S0_min_epsilon;
			
			//printf("Path:%d,TimeStep:%d: D_tmp = %f %f = D_tmp_pls_epsilon \n",path,tmp_step_noX, D_tmp, D_tmp_pls_epsilon);



			delta = (D_tmp_pls_epsilon - D_tmp )/epsilon;

			gamma = (D_tmp_pls_epsilon - 2.*D_tmp + D_tmp_min_epsilon)/(epsilon*epsilon);
			//gamma1 = (delta-delta_old)/(S_tmp_mu - S_tmp_mu_alt)
			gamma1 = (delta-delta_old)/(S_tmp_mu - S_tmp_mu_alt);

			//RiskFreeAsset_beforeRebalancing = (delta_old * S_tmp_mu)-D_tmp;
			RiskFreeAsset_beforeRebalancing = -(D_tmp_alt -(delta_old * S_tmp_mu_alt))*exp(r_f*delta_t);// Hinzufügen von minus, weil rsik free asset = kredit

			gap_Value = (delta_old * S_tmp_mu) - RiskFreeAsset_beforeRebalancing -D_tmp; // Minus für risk free asset, weil kredit als negativer wert gezählt wurde, Negativer Wert = Höhe der nicht abgesicherten position

			kosten_aktien_Kauf_Verkauf = (delta_old -delta)*S_tmp_mu;// in t_0 delta_old = 0 //  Kosten sind negativ, Erlöse positiv
			Sum_kosten_aktien_Kauf_Verkauf += kosten_aktien_Kauf_Verkauf;


			Sum_zinskostenAktien += zinskostenAktien * (exp(r_f*delta_t)-1.); // Zinseszinseffekt auf Zinskonto ??
			//Jetzt erst neue Zinsen aufs konto dazuaddieren. In einem jahr wird konto wieder verzinst...
			zinskostenAktien = - delta *S_tmp_mu*(exp(r_f*delta_t)-1.); // Kosten sind negativ. Das minus wurde hinzugefuegt, damit Kosten als negative Werte auftreten, wegen des negativen Zinses
			Sum_zinskostenAktien += zinskostenAktien;

			delta_old = delta;
			S_tmp_mu_alt = S_tmp_mu;
			D_tmp_alt = D_tmp;
			RiskFreeAsset_afterRebalancing = -(D_tmp - (delta * S_tmp_mu));// nagative Werte, wenn 

			check_gap_value_after_rebalancing = (delta*S_tmp_mu)-RiskFreeAsset_afterRebalancing - D_tmp;

			//printf("Path:%d,TimeStep:%d: Aktkurs = %f, Derivat = %f, delta = %f, repli_portf_vor = %f, repli_portf_nach = %f, kosten_akti = %f, zinskostAkt = %f \n",path,tmp_step_noX, S_tmp_mu, D_tmp, delta, RiskFreeAsset_beforeRebalancing, RiskFreeAsset_afterRebalancing, kosten_aktien_Kauf_Verkauf, zinskostenAktien );
			//printf("Sum_kosten_aktien_Kauf_Verkauf = %f, Sum_zinskostenAktien =%f\n", Sum_kosten_aktien_Kauf_Verkauf, Sum_zinskostenAktien);

 		 file_output_Delta << time << "," << S_tmp_mu << "," << D_tmp << "," << D_tmp_pls_epsilon << "," << D_tmp_min_epsilon << "," << delta << "," << gamma << "," 
 		 << gamma1 << "," << kosten_aktien_Kauf_Verkauf << "," << zinskostenAktien << "," << RiskFreeAsset_beforeRebalancing << "," 
 		 << RiskFreeAsset_afterRebalancing << "," << gap_Value << "," << check_gap_value_after_rebalancing << "," << Sum_kosten_aktien_Kauf_Verkauf << "," << Sum_zinskostenAktien << "\n";
		}// end loop through time steps
		//file_output_Delta << "\n";

		Kosten_Hedge = (Sum_zinskostenAktien + Sum_kosten_aktien_Kauf_Verkauf - delta* S_tmp_mu)*(exp(-r_f*maturity));// discount for one year since maturity is one year
		file_output_Delta << "FV_Zinskosten" << "," << "FV_AktienKaufVerkauf" << "," << "FV_delta_S_tmp_muT" << "," << "PV_KostenHedge" << "\n";
		file_output_Delta << Sum_zinskostenAktien << "," << Sum_kosten_aktien_Kauf_Verkauf << "," << delta * S_tmp_mu << "," << Kosten_Hedge << "\n";
		file_output_Delta.close();

		file_Kosten_Delta_Hedging << path << "," << Sum_zinskostenAktien << "," << Sum_kosten_aktien_Kauf_Verkauf << "," << delta * S_tmp_mu << "," << Kosten_Hedge << "\n";

		//printf("Kosten des Hedge = %f, Delta = %f, S_tmp_mu = %f, Sum_zinskostenAktien = %f, Sum_kosten_aktien_Kauf_Verkauf = %f\n", Kosten_Hedge, delta, S_tmp_mu, Sum_zinskostenAktien, Sum_kosten_aktien_Kauf_Verkauf);
	}//end loop through paths

	file_Kosten_Delta_Hedging.close();

	//std::ifstream file1("DeltaX1.csv");
	//std::string line1;
/*
	std::size_t a, b, c;
	int counter = 0;
	std::ifstream file2("DeltaX2.csv");
	std::string line2;

	for(int tmp_step_noX1 = 1; tmp_step_noX1 <= N_timeSteps; tmp_step_noX1++){
		printf("Counter1 = %d\n",tmp_step_noX1);
		double tmp_value_file = 0.0;
		double sum_column = 0.0;
		double mean_column = 0.0;
		std::size_t a, b, c;
		for(int path1 = 1; path1 <= N_paths; path1++){
			printf("Hier gewesen\n");
			a = line2.find(',', tmp_step_noX1);
			b = line2.find(',', a + 1);
			c = line2.find('\"', b + 1);
			
			std::string tmp_string;
			if(b != std::string::npos){
			 tmp_string = line2.substr(a+1, b-1-a);
			 printf("%s\n", tmp_string.c_str());
			}
			tmp_value_file = std::stof(tmp_string);
			sum_column += tmp_value_file;
		}
			mean_column = sum_column/N_paths;
		printf("Mean Column = %f\n", mean_column);
	}
		file2.close();
*/

/*

	while(getline(file1, line1)){
 
		printf("Counter = %d\n",counter);
		double tmp_value_file = 0.0;
		double sum_column = 0.0;
		double mean_column = 0.0;
	std::ifstream file2("DeltaX1.csv");
	std::string line2;
		while(getline(file2, line2)){
			printf("Hier gewesen\n");
			a = line2.find(',', counter);
			b = line2.find(',', a + 1);
			c = line2.find('\"', b + 1);
			
			std::string tmp_string;
			if(b != std::string::npos){
			 tmp_string = line2.substr(a+1, b-1-a);
			 printf("%s\n", tmp_string.c_str());
			}
			//else if(b != std::string::npos){
			//  tmp_string = line2.substr(b, 1);
			//  printf("%s\n", tmp_string.c_str());
			//}
			
			tmp_value_file = std::stof(tmp_string);
			
			sum_column += tmp_value_file;
		}
		file2.close();
		mean_column = sum_column/N_paths;
		printf("Mean Column = %f\n", mean_column);
		file2.clear();                 // clear fail and eof bits
    file2.seekg(0, std::ios::beg); // back to the start!
   	counter++;
	}
		file1.close();
	*/

}

//_____________________________________________________________________________________________________________________________________________________________________________________________________
//

PriceBonusDerivat_pls_Delta2 Evaluate_BonusDerivat_loop_for_delta21(const double S_0, const double S_0_pls_epsilon, const double S_0_min_epsilon, const double Sigma, const double r_f, const double Barrier, const double Bonus_Level, const int maturity, const int N_timeSteps, const int N_inner_paths, int current_time, bool Barrier_already_touched_S0, bool Barrier_already_touched_S0_pls_epsilon, bool Barrier_already_touched_S0_min_epsilon, bool DividendsPaid, const double Dividend ){
	
	double delta_t = (double)maturity/N_timeSteps; // time intervall for each evaluation
	double price_S0 = 0.0;
	double price_S0_pls_epsilon = 0.0;
	double price_S0_min_epsilon = 0.0;
	double tmp_sum_price_S0 = 0.0;
	double tmp_sum_price_S0_pls_epsilon = 0.0;
	double tmp_sum_price_S0_min_epsilon = 0.0;
	// Hier wird der bool value zur Berücksichtigung der Historie auf den Wert bei t_1 gesetzt
	bool TouchBarrier_S0 = Barrier_already_touched_S0;// Boolian for checking whether barrier was touched or undershoot
	bool TouchBarrier_S0_pls_epsilon = Barrier_already_touched_S0_pls_epsilon;// Boolian for checking whether barrier was touched or undershoot
	bool TouchBarrier_S0_min_epsilon = Barrier_already_touched_S0_min_epsilon;// Boolian for checking whether barrier was touched or undershoot
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
		double S_tmp_S0_min_epsilon = S_0_min_epsilon;
		double S_tmp_initial_S0 = S_0;
		double S_tmp_initial_S0_pls_epsilon = S_0_pls_epsilon;
		double S_tmp_initial_S0_min_epsilon = S_0_min_epsilon;
		double random_number;
		// Hier wird der bool value zur Berücksichtigung der Historie auf den Wert bei t_1 gesetzt
		// Diese Initialisierung wird für jeden pfad durchlauf vorgenommen

		TouchBarrier_S0 = Barrier_already_touched_S0;// Boolian for checking whether barrier was touched or undershoot
		TouchBarrier_S0_pls_epsilon = Barrier_already_touched_S0_pls_epsilon;// Boolian for checking whether barrier was touched or undershoot
		TouchBarrier_S0_min_epsilon = Barrier_already_touched_S0_min_epsilon;// Boolian for checking whether barrier was touched or undershoot


		double tmp_price_S0 = 0.0;
		double tmp_price_S0_pls_epsilon = 0.0;
		double tmp_price_S0_min_epsilon = 0.0;
//printf("Innerer Pfad = %d\n", path);
		for(int tmp_step_noX = current_time; tmp_step_noX < N_timeSteps; tmp_step_noX++){

			if(tmp_step_noX > N_timeSteps/2. && !DividendsPaid){
				DividendsPaid = true;
				S_tmp_initial_S0 -= Dividend;
				S_tmp_initial_S0_pls_epsilon -= Dividend;
				S_tmp_initial_S0_min_epsilon -= Dividend;
			}

			random_number = Get_StdNormDistr_RandomNumber();    // Get the standard normal distributed random number                                                                    
			S_tmp_S0             = S_tmp_initial_S0 * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
			S_tmp_S0_pls_epsilon = S_tmp_initial_S0_pls_epsilon * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
			S_tmp_S0_min_epsilon = S_tmp_initial_S0_min_epsilon * exp( (r_f - (Sigma*Sigma)/2.) * delta_t + Sigma * sqrt(delta_t) * random_number );   // Calculate the price of the underlying at the n-st time step: delta_t* tmp_step_noX
			//printf("Verhältnis = %f\n", S_tmp_S0/S_tmp_S0_pls_epsilon);
			

			S_tmp_initial_S0 = S_tmp_S0;
			S_tmp_initial_S0_pls_epsilon = S_tmp_S0_pls_epsilon;
			S_tmp_initial_S0_min_epsilon = S_tmp_S0_min_epsilon;


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

			if(!TouchBarrier_S0_min_epsilon && S_tmp_S0_min_epsilon <= Barrier) {
				TouchBarrier_S0_min_epsilon = true;                                        // If barrier was touched set the boolian to kTRUE 
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

		if(TouchBarrier_S0_min_epsilon || (S_tmp_S0_min_epsilon > Bonus_Level) )  tmp_price_S0_min_epsilon = S_tmp_S0_min_epsilon*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);                       // If barrier was (touched or undershoot) or if (the value of the underlying S_t is greater than the Bonuslevel Bonus_Level) than return value of the underlying = Price of Bonus certificate
		else tmp_price_S0_min_epsilon = Bonus_Level*exp(-r_f*(1.-(double)current_time/N_timeSteps)*maturity);		

		tmp_sum_price_S0 += tmp_price_S0;
		tmp_sum_price_S0_pls_epsilon += tmp_price_S0_pls_epsilon;
		tmp_sum_price_S0_min_epsilon += tmp_price_S0_min_epsilon;

	}


	price_S0 = tmp_sum_price_S0/N_inner_paths;
	price_S0_pls_epsilon = tmp_sum_price_S0_pls_epsilon/N_inner_paths;
	price_S0_min_epsilon = tmp_sum_price_S0_min_epsilon/N_inner_paths;

	tmp_struct.priceDerivat_S0 = price_S0;
	tmp_struct.priceDerivat_S0_pls_epsilon = price_S0_pls_epsilon;
	tmp_struct.priceDerivat_S0_min_epsilon = price_S0_min_epsilon;


	//printf("3.  Return D_S0 = %f, Boolian_S0 = %i\n", tmp_struct.priceDerivat_S0, tmp_struct.BarrierTouched_S0);
	//printf("3.  Return D_ep = %f, Boolian_ep = %i\n", tmp_struct.priceDerivat_S0_pls_epsilon, tmp_struct.BarrierTouched_S0_pls_epsilon );

/*
Frage: Welchen Bool wert muss man für den nächsten aüßeren Time step übergeben. Man sollte ja nicht in dieser Function den bool wert eines pfades an den nächsten übergeben
da diese ja unabhängig voneinander sind, sondern man sollte ihn für jeden neuen pfad wieder mit dem übergebenen wert aus t_-1 d.h. aus dem vorherigen Funktionsausruf übergeben
(Um aber die Historie des derivats zu berücksichtigen muss man einen wert als return value haben um ihn für den nächsten aufruf der Funktion
an einem zeitpunkt delta_t später zu haben.)
Vorschlag: Den letzen bool wert des letzten inneren pfades übergeben, Da diese Funktion innerhalb der äußeren Funktion Calculate_Delta() ja N_path mal dieselben Rechnungen durchführt 
kann hat man so unter umständen eine streueng drinne.
Frage was ist das maximale delta das möglich ist
*/

	return tmp_struct;

}






