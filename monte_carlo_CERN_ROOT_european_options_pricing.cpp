#include <iostream>
#include <random>
#include <iomanip>
#include "TFile.h" // .root specific libraries
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"

using std::cout;
using std::endl;


// structure to calculate price of a european call or put option via monte-carlo
struct monte_carlo_options {


        // constructor
        monte_carlo_options(int number_of_events_in, double spot_in, double strike_in, double risk_free_rate_in,
                            double volatility_in, double time_to_maturity_in){

                number_of_events = number_of_events_in;
                So = spot_in;
                St = strike_in;
                r = risk_free_rate_in;
                v = volatility_in;
                t = time_to_maturity_in;
        }


        // function to estimate option(s) value
        void estimate_option(){

                // create root file
                create_root_file();

                // generate random numbers and create a normal distribution (uses Mersenne Twister PRNG), creating (pseudo) brownian motion, Wt
                std::random_device local_rand;
                std::mt19937 local_gen(local_rand());
                std::normal_distribution<double> my_normal_distro(0,1);

                // create .root histograms
                create_root_distributions();

                // loop for all MC events
                for(int i = 0; i < number_of_events; i++) {

                        // update brownian motion
                        Wt = my_normal_distro(local_gen);

                        // calculate drift, diffusion, and then the spot price at maturity
                        drift_component = (r - (0.5 * v * v)) * t;
                        diffusion_component = v * std::sqrt(t) * Wt;
                        estimated_price_at_maturity = So * std::exp(drift_component + diffusion_component);

                        // determine payoff
                        single_call_payoff = get_payoff(estimated_price_at_maturity, St, true);
                        single_put_payoff = get_payoff(estimated_price_at_maturity, St, false);

                        // fill .root histograms
                        fill_root_distributions();

                }

                // verify MC model
                verify_MC();

                // print MC parameters, mean payoff and present values
                print_answer();

                // write to .root file
                fRootHis->Write();

                // manage memory, note: it is better to use smart pointers, however, I have used raw pointers so this can run on older versions of C++
                manage_memory();

        } // estimate_option


        // destructor
        ~monte_carlo_options() = default;


// protected
protected:

        // protected parameters
        char fileRoot[100];
        int number_of_events{};
        double So{}, St{}, r{}, v{}, t{}, drift_component{}, diffusion_component{}, estimated_price_at_maturity{}, single_call_payoff{}, single_put_payoff{}, Wt{};
        double MC_mean{}, MC_mean_error{}, MC_sigma{}, MC_sigma_error{}, MC_chiSquare{}, mean_call_payoff{}, mean_put_payoff{}, present_call_value{}, present_put_value{};
        TH1D * MC_values; TH1D * estimated_price_at_maturity_distro; TH1D * payoff_call_distro; TH1D * payoff_put_distro; // .root
        TF1 * f1; // .root
        TFile * fRootHis; // .root


        // function to create .root file
        void create_root_file(){
                sprintf(fileRoot,"european_options_plots.root");
                fRootHis = new TFile(fileRoot,"RECREATE");
        }


        // function to create distributions in CERN ROOT
        void create_root_distributions(){

                // distributions for monte-carlo values from N(0,1), calculated price at maturity, and payoff
                MC_values = new TH1D("MC_values","MC_values",200,-10,10);
                estimated_price_at_maturity_distro = new TH1D("estimated_price_at_maturity", "estimated_price_at_maturity", 1000000, 0, 100000);
                payoff_call_distro = new TH1D("payoff_call", "payoff_call", 100000, -10, 100000);
                payoff_put_distro = new TH1D("payoff_put", "payoff_put", 100000, -10, 100000);
        }


        // function to fill distributions in CERN ROOT
        void fill_root_distributions(){

                MC_values->Fill(Wt,1);
                estimated_price_at_maturity_distro->Fill(estimated_price_at_maturity,1);
                payoff_call_distro->Fill(single_call_payoff,1);
                payoff_put_distro->Fill(single_put_payoff,1);
        }


        // function to calculate payoff from estimated price
        double get_payoff(double estimated_price_at_maturity, double St, double call){

                if(call) { return std::max(estimated_price_at_maturity - St, 0.0); }
                return std::max(St - estimated_price_at_maturity, 0.0);
        }


        // function to verify MC brownian motion model N(0,1) parameters using CERN ROOT
        void verify_MC(){

                f1 = new TF1("f1","gaus",-10,10);
                f1->SetParameter(1,0);
                MC_values->Fit("f1","QR");

                MC_mean = f1->GetParameter(1);
                MC_mean_error = f1->GetParError(1);
                MC_sigma = f1->GetParameter(2);
                MC_sigma_error = f1->GetParError(2);
                MC_chiSquare = (f1->GetChisquare())/(f1->GetNDF());
        }


        // function to print MC parameters, mean payoff and present values
        void print_answer(){

                mean_call_payoff = payoff_call_distro->GetMean();
                present_call_value = mean_call_payoff * std::exp(-r * t);

                mean_put_payoff = payoff_put_distro->GetMean();
                present_put_value = mean_put_payoff * std::exp(-r * t);

                cout << "\n" << "MC parameters: " << "\n"
                     << " reduced chi-square of fit = " << MC_chiSquare << "\n"
                     << " mean = " << MC_mean << " +/- " << MC_mean_error << "\n"
                     << " sigma = " << MC_sigma << " +/- " << MC_sigma_error << "\n"
                     << "\n" << "Call option details: " << "\n"
                     << " estimated payoff: " << mean_call_payoff << "\n"
                     << " present value: " << present_call_value << "\n"
                     << "\n" << "Put option details: " << "\n"
                     << " estimated payoff: " << present_put_value << "\n"
                     << " present value: " << present_put_value << "\n"
                     << "\n" << "Time to maturity: " << (t*365) << " days" << "\n" << endl;

        }


        // function for memory management
        void manage_memory(){
                delete MC_values; delete estimated_price_at_maturity_distro; delete payoff_call_distro; delete payoff_put_distro;
        }


}; // struct


// main
int main(){

        // ----- ----- -----
        //
        // READ ME
        //
        // ----- ----- -----
        // notes on method:
        //
        // This code allows the price of a european call or put option to be calculated,
        // the model determins the future price of an option using a solution to the Black-Scholes-Merton SDE
        // given by: St = S0 * exp((r - 0.5*v^2)T + v*sqrt(T)*Wt), where,
        // St = underlying stock price at maturity, S0 = current underlying stock price, r = risk free rate,
        // v = volatility, T = time to maturity, Wt = brownian motion,
        // the pay-offs are then then calculated using monte-carlo, averaged, and the present/discount value (option price) is determined.
        // This model therefore assumes a constant risk free rate and a constant volatility.
        //
        // ----- ----- -----
        // notes on code:
        //
        // this code was created using the C++ 11 STL, CERN ROOT 6.14/02 and was tested using LLVM 9.0.0,
        // to run this code, specify the parameters, as below, compile and run.
        // the relevent value will be printed to the terminal and plots will be produced in a .root file for analysis
        //
        // ----- ----- -----
        // parameters:
        //
        double spot = 100; // spot price (current price) of the underlying asset on the market, e.g., 1 share of NVDA
        double strike = 110; // strike price of the option, i.e., the agreed upon price (of the underlying asset) at which the option can be exercised (traded)
        double risk_free_rate = 0.10; // risk free interest rate, expressed as a percentage, i.e., the rate of return on a zero-risk asset, e.g., a goverment bond
        double volatility = 0.30; // volatility of underlying asset, it's essentially the standard deviation of the daily price changes of underlying asset
        double time_to_maturity = 365.0 / 365.0; // time to maturity, i.e., time untill option is exercised/expires, expressed in years, days/365
        int number_of_events = 1000000; // number of monte-carlo events to estimate option price (single thread)
        //
        // ----- ----- -----
        // run code:
        //
        // create structure...
        monte_carlo_options estimate_001 = monte_carlo_options(number_of_events, spot, strike, risk_free_rate, volatility, time_to_maturity);
        //
        // call member function to calculate values
        estimate_001.estimate_option();
        //
        //
        cout << "done." << endl;
        return 0;
}
