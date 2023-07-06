#include <iostream>
#include <random>
#include <iomanip>
#include <thread>

using std::cout;
using std::endl;



// structure to calculate price of a european call or put option via monte-carlo
struct monte_carlo_options {

        // constructor
        monte_carlo_options(int number_of_events_in, double spot_in, double strike_in, double risk_free_rate_in,
                            double volatility_in, double time_to_maturity_in, int thread_in, double * call_payoff_ptr_in,
                            double * call_value_ptr_in, double * put_payoff_ptr_in, double * put_value_ptr_in){

                number_of_events = number_of_events_in;
                So = spot_in;
                St = strike_in;
                r = risk_free_rate_in;
                v = volatility_in;
                t = time_to_maturity_in;
                thread = thread_in;
                call_payoff_ptr_inside = call_payoff_ptr_in;
                call_value_ptr_inside = call_value_ptr_in;
                put_payoff_ptr_inside = put_payoff_ptr_in;
                put_value_ptr_inside = put_value_ptr_in;
        }

        // function call operator - special function for multi threading, this needs to be defined within the strtucture/class,
        // and its scope defines what is executed when the structure is called on a thread
        void operator()(){

                estimate_option();
        }

        // monte-carlo function to estimate option(s) value
        void estimate_option(){

                // // generate random numbers and create a normal distribution (uses Mersenne Twister PRNG), creating (pseudo) brownian motion, Wt
                std::random_device local_rand;
                std::mt19937 local_gen(local_rand());
                std::normal_distribution<double> my_normal_distro(0,1);

                // loop for all MC events
                for(int i = 0; i < number_of_events; i++) {

                        // update brownian motion
                        Wt = my_normal_distro(local_gen);

                        // calculate drift, diffusion, and then the spot price at maturity
                        drift_component = (r - (0.5 * v * v)) * t;
                        diffusion_component = v * std::sqrt(t) * Wt;
                        estimated_price_at_maturity = So * std::exp(drift_component + diffusion_component);

                        // sum payoffs
                        cumulative_call_payoff += get_payoff(estimated_price_at_maturity, St, true);
                        cumulative_put_payoff += get_payoff(estimated_price_at_maturity, St, false);

                }

                // grab data from thread
                grab_data();

        } // estimate_option


        // destructor
        ~monte_carlo_options() = default;


// protected
protected:

        // protected parameters
        int number_of_events{}, thread{};
        double So{}, St{}, r{}, v{}, t{}, drift_component{}, diffusion_component{}, estimated_price_at_maturity{}, Wt{};
        double mean_call_payoff{}, mean_put_payoff{}, present_call_value{}, present_put_value{};
        double cumulative_call_payoff{}, cumulative_put_payoff{};
        double * call_payoff_ptr_inside; double * call_value_ptr_inside; double * put_payoff_ptr_inside; double * put_value_ptr_inside;


        // function to calculate payoff from estimated price
        double get_payoff(double estimated_price_at_maturity, double St, double call){

                if(call) { return std::max(estimated_price_at_maturity - St, 0.0); }
                return std::max(St - estimated_price_at_maturity, 0.0);
        }


        // function to grab data and send to main stack
        void grab_data(){

                // calculate payoff values
                mean_call_payoff = cumulative_call_payoff / number_of_events;
                mean_put_payoff = cumulative_put_payoff / number_of_events;

                // calculate present values
                present_call_value = mean_call_payoff * std::exp(-r * t);
                present_put_value = mean_put_payoff * std::exp(-r * t);

                // assign values to main stack (outside of thread)
                * call_payoff_ptr_inside = mean_call_payoff;
                * call_value_ptr_inside = present_call_value;
                * put_payoff_ptr_inside = mean_put_payoff;
                * put_value_ptr_inside = present_put_value;

        }

}; // struct



// function to activate multi-thread, start monte-carlo process and print results 
void run_multi_thread(double spot_in, double strike_in, double risk_free_rate_in, double volatility_in,
                      double time_to_maturity_in, int number_of_events_per_thread_in, int number_of_threads_in){


        // initialise values which will be updated from threads 
        double call_payoff[number_of_threads_in];
        double call_value[number_of_threads_in];
        double put_payoff[number_of_threads_in];
        double put_value[number_of_threads_in];

        // initialise pointers to access above variables
        double * call_payoff_ptr[number_of_threads_in];
        double * call_value_ptr[number_of_threads_in];
        double * put_payoff_ptr[number_of_threads_in];
        double * put_value_ptr[number_of_threads_in];

        // define variables to sum values grabbed from threads 
        double cumulative_call_payoff{}, cumulative_call_present{}, cumulative_put_payoff{}, cumulative_put_present{};


        // create an array for the monte-carlo objects (1 object per thread) and an array for the threads
        monte_carlo_options * array_of_objects[number_of_threads_in];
        std::thread * array_of_threads[number_of_threads_in];

        // assign pointers to initialised values which will be updated from threads 
        for(int i = 0; i < number_of_threads_in; i++) {

                call_payoff_ptr[i] = &call_payoff[i];
                call_value_ptr[i] = &call_value[i];
                put_payoff_ptr[i] = &put_payoff[i];
                put_value_ptr[i] = &put_value[i];
        }

        // create the monte-carlo object and activate thread for specified number of threads
        for(int i = 0; i < number_of_threads_in; i++) {

                // create monte-carlo object 
                array_of_objects[i] = new monte_carlo_options(number_of_events_per_thread_in, spot_in, strike_in, risk_free_rate_in, volatility_in, time_to_maturity_in,
                                                              i, call_payoff_ptr[i], call_value_ptr[i], put_payoff_ptr[i], put_value_ptr[i]);
                
                // send monte-carlo object to thread and activate
                array_of_threads[i] = new std::thread(*array_of_objects[i]);
        }

        // wait untill all threads have finished processing...
        for(int i = 0; i < number_of_threads_in; i++) { array_of_threads[i]->join(); }

        // grab data from all threads
        for(int i = 0; i < number_of_threads_in; i++) {

                // call payoff
                cumulative_call_payoff += call_payoff[i];
                // call present value
                cumulative_call_present += call_value[i];
                // put payoff
                cumulative_put_payoff += put_payoff[i];
                // put present value
                cumulative_put_present += put_value[i];

        }


        // print results to terminal
        cout << "Call option: " << "\n"
                " payoff " << cumulative_call_payoff/number_of_threads_in << "\n" <<
                " present value " << cumulative_call_present/number_of_threads_in << "\n"
                "Put option: " << "\n" <<
                " payoff " << cumulative_put_payoff/number_of_threads_in << "\n" <<
                " present value " << cumulative_put_present/number_of_threads_in << endl;


        // release memory
        for (int i = 0; i < number_of_threads_in; ++i) {
                delete array_of_objects[i]; delete array_of_threads[i];
        }

}


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
        // v = volatility, T = time to maturity, Wt = brownian motion.
        // The pay-offs are then calculated using monte-carlo, averaged, and the present value (option price) is determined.
        // This model therefore assumes a constant risk free rate and a constant volatility.
        //
        // ----- ----- -----
        // notes on code:
        //
        // This code was created using C++ 11 STL and was tested using LLVM 9.0.0,
        // to run this code, specify the parameters, as below, compile and run.
        // The code uses a multi-threading process in which you can specify how many threads you wish to use via 'number_of_threads',
        // each thread will run a monte-carlo process, with the number of events as specifed in 'number_of_events_per_thread',
        // once the monte-carlo process' finish the results are compiled and the options: call payoff, put payoff, call present value,
        // and put present value are printed to the terminal.
        //
        // ----- ----- -----
        // parameters:
        //
        double spot = 100; // spot price (current price) of the underlying asset on the market, e.g., 1 share of NVDA
        double strike = 110; // strike price of the option, i.e., the agreed upon price (of the underlying asset) at which the option can be exercised (traded)
        double risk_free_rate = 0.10; // risk free interest rate, expressed as a percentage, i.e., the rate of return on a zero-risk asset, e.g., a goverment bond
        double volatility = 0.30; // volatility of underlying asset, it's essentially the standard deviation of the daily price changes of underlying asset
        double time_to_maturity = 365.0 / 365.0; // time to maturity, i.e., time untill option is exercised/expires, expressed in years, days/365
        int number_of_events_per_thread = 1000000; // number of monte-carlo events (per thread) to estimate option price
        int number_of_threads = 4; // number of threads to be used, note: each thread will start a monte-carlo process, hence you will have a process running on each thread
        //
        // ----- ----- ----- 
        // Run code:
        //
        run_multi_thread(spot, strike, risk_free_rate, volatility, time_to_maturity, number_of_events_per_thread, number_of_threads);
        //
        //
        //  
        cout << "done." << endl;
        return 0;
}
