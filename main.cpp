#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <chrono>
#include <cmath>


class Event{
protected:
    std::default_random_engine m_generator;
public:
    Event(std::default_random_engine& generator):m_generator(generator){}
    virtual double rate(const double time, const std::vector<int>& state) = 0;
    virtual void addEventVector(std::vector<int>& state, const double scale) = 0;
};

class ReproduceConnected:public Event{
    double m_b_RC;
    double m_m;
    double m_alpha;
public:
    ReproduceConnected(std::default_random_engine& generator, double b_RC, double m, double alpha):Event(generator),m_b_RC(b_RC),m_m(m), m_alpha(alpha){}
    virtual double rate(const double time, const std::vector<int>& state){
        double sum = 0;
        for(auto it = ++state.begin();it<state.end();it++){sum+=*it;}
        return (1-m_m)*m_b_RC*std::pow(sum,m_alpha);
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        std::discrete_distribution<int> eventSampler(++state.begin(),state.end());
        int event = eventSampler(m_generator);
        state[event]+=scale;
    }
};

class ReproduceFree:public Event{
    double m_b_F;
    double m_m;
    double m_alpha;
public:
    ReproduceFree(std::default_random_engine& generator, double b_F, double m,double alpha):Event(generator),m_b_F(b_F),m_m(m),m_alpha(alpha){}
    virtual double rate(const double time, const std::vector<int>& state){
        return (1-m_m)*m_b_F*std::pow(state[0],m_alpha);
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        state[0]+=scale;
    }
};

class Switch_CF:public Event{
    double m_s_CF;
    double m_alpha;
public:
    Switch_CF(std::default_random_engine& generator,double s_CF,double alpha):Event(generator), m_s_CF(s_CF), m_alpha(alpha){}
    virtual double rate(const double time, const std::vector<int>& state){
        return m_s_CF*std::pow(state[0],m_alpha);
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        state[1]-=scale;
        state[0]+=scale;
    }
};
class Switch_FC:public Event{
    double m_s_FC;
    double m_alpha;
    double m_alpha_C;
public:
    Switch_FC(std::default_random_engine& generator,double s_FC,double alpha, double alpha_C):Event(generator), m_s_FC(s_FC), m_alpha(alpha), m_alpha_C(alpha_C){}
    virtual double rate(const double time, const std::vector<int>& state){
        return m_s_FC*std::pow(state[0],m_alpha_C)*std::pow(state[1],m_alpha);
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        state[0]-=scale;
        state[1]+=scale;

    }
};
class DieConnected:public Event{
    double m_d_C;
    double m_c_CC;
    double m_c_CF;
    double m_alpha_C;
    double m_alpha_F;
public:
    DieConnected(std::default_random_engine& generator,double d_C, double c_CC, double c_CF,double alpha_C,double alpha_F):Event(generator),m_d_C(d_C),m_c_CC(c_CC),m_c_CF(c_CF),m_alpha_C(alpha_C),m_alpha_F(alpha_F){}
    virtual double rate(const double time, const std::vector<int>& state){
        double sum = 0;
        for(auto it = ++state.begin();it<state.end();it++){sum+=*it;}
        return (m_d_C+m_c_CC*std::pow(sum,m_alpha_C)+m_c_CF*std::pow(state[0],m_alpha_F))*sum;
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        std::discrete_distribution<int> eventSampler(++state.begin(),state.end());
        int event = eventSampler(m_generator);
        state[event+1]-=scale;
    }
};
class DieFree:public Event{
    double m_d_F;
    double m_c_FF;
    double m_c_FC;
    double m_alpha_F;
    double m_alpha_C;
public:
    DieFree(std::default_random_engine& generator,double d_F, double c_FF, double c_FC, double alpha_F, double alpha_C):Event(generator),m_d_F(d_F),m_c_FF(c_FF),m_c_FC(c_FC), m_alpha_F(alpha_F), m_alpha_C(alpha_C){}
    virtual double rate(const double time, const std::vector<int>& state){
        double sum = 0;
        for(auto it = ++state.begin();it<state.end();it++){sum+=*it;}
        return (m_d_F+m_c_FF*std::pow(state[0],m_alpha_F)+m_c_FC*std::pow(sum,m_alpha_C))*state[0];
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        state[0]-=scale;
    }
};

class AddConnection:public Event{
    double m_p;
    double m_alpha;
    int m_position;
public:
    AddConnection(std::default_random_engine& generator,double p, double alpha, int position):Event(generator),m_p(p), m_alpha(alpha), m_position(position){}
    virtual double rate(const double time, const std::vector<int>& state){
        return m_p*std::pow(state[m_position],m_alpha);
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        state[m_position]-=scale;
        state[m_position+1]+=scale;
    }
};

class RemoveConnection:public Event{
    double m_p;
    double m_alpha;
    int m_position;
public:
    RemoveConnection(std::default_random_engine& generator,double p, double alpha, int position):Event(generator),m_p(p), m_alpha(alpha), m_position(position){}
    virtual double rate(const double time, const std::vector<int>& state){
        return m_p*std::pow(state[m_position],m_alpha);
    }
    virtual void addEventVector(std::vector<int>& state, const double scale){
        state[m_position]-=scale;
        state[m_position-1]+=scale;
    }
};
std::vector<std::vector<int>> run_simulation(std::default_random_engine& generator, int k_max, double t_max, std::vector<int> n, double step_size, std::vector<std::unique_ptr<Event>>& events, double K){
    int k = 0;
    double t_k =0;
    bool cellAvailable = false;
    for(const auto& el: n){
        cellAvailable = el>0;
        if(cellAvailable) break;
    }
    std::vector<std::vector<int>> nl{n};
    int i=0;
    while(k<k_max and t_k<t_max and cellAvailable){
        double t = 0;
        std::vector<int> n_tmp = nl[i];
        while(t<step_size && k<k_max){
            double R_tot = 0;
            std::vector<double> rates;
            for(const auto& event: events){
                double rate = event->rate(t_k+t,n_tmp);
                rates.push_back(rate);
                R_tot += rate;
            }
            std::exponential_distribution<double> stepSampler(K*R_tot);
            double step = stepSampler(generator);
            if (t+step < step_size){
                std::discrete_distribution<int> eventSampler(rates.begin(),rates.end());
                int event = eventSampler(generator);
                events[event]->addEventVector(n_tmp,K);
            }
            t = t+step;
            k++;
        }
        nl.push_back(n_tmp);
        cellAvailable = false;
        for(const auto& el: n){
            cellAvailable = el>0;
            if(cellAvailable) break;
        }
        t_k=t_k+step_size;
        i++;
    }
    return nl;
}


int main() {
    long long seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    std::ofstream seedFile;
    seedFile.open("./seeds.txt" ,std::ios_base::app);
    auto timepoint = std::chrono::system_clock::now();
    std::time_t timepoint_t = std::chrono::system_clock::to_time_t(timepoint);

    char timeStr[100];
    if (std::strftime(timeStr, sizeof(timeStr), "%d/%m/%Y %T: ", std::localtime(&timepoint_t))) {
        seedFile << timeStr;
    }
    seedFile << seed << "\n";
    seedFile.close();

    std::ifstream i("modelParameters.json");
    nlohmann::json j;
    i >> j;
        int k_max = j["k_max"];
        double t_max = j["t_max"];
        std::vector<int> n;
        for(int i = 0;i<int(j["maxConnections"])+1;i++){
            n.push_back(j["start"][i]);
        }
        double step_size = j["step_size"];
        double K = j["K"];
        std::vector<std::unique_ptr<Event>> events;

        events.push_back(std::move(std::make_unique<ReproduceConnected>(generator, j["ReproduceConnected"]["b"],
                                                                        j["ReproduceConnected"]["m"],
                                                                        j["ReproduceConnected"]["alpha"])));
        events.push_back(std::move(
                std::make_unique<ReproduceFree>(generator, j["ReproduceFree"]["b"], j["ReproduceFree"]["m"],
                                                j["ReproduceFree"]["alpha"])));

    if(j["model"] == "Kraut") {
        events.push_back(
                std::move(std::make_unique<Switch_CF>(generator, j["Switch_CF"]["s_CF"], j["Switch_CF"]["alpha"])));
        events.push_back(std::move(
                std::make_unique<Switch_FC>(generator, j["Switch_FC"]["s_FC"], j["Switch_FC"]["alpha"],
                                            j["Switch_FC"]["alpha_C"])));
    }
        events.push_back(std::move(
                std::make_unique<DieConnected>(generator, j["DieConnected"]["d_C"], j["DieConnected"]["c_CC"],
                                               j["DieConnected"]["c_CF"], j["DieConnected"]["alpha_C"],
                                               j["DieConnected"]["alpha_F"])));
        events.push_back(std::move(
                std::make_unique<DieFree>(generator, j["DieFree"]["d_F"], j["DieFree"]["c_FF"], j["DieFree"]["c_FC"],
                                          j["DieFree"]["alpha_F"], j["DieFree"]["alpha_C"])));
    if(j["model"] == "connections")
    {
        int maxConnections = j["maxConnections"];
        events.push_back(std::move(
                std::make_unique<AddConnection>(generator, j["addConnections"][0]["p"],j["addConnections"][0]["alpha"],0)));

        for(int i = 1;i<maxConnections;i++){
            events.push_back(std::move(
                    std::make_unique<AddConnection>(generator, j["addConnections"][i]["p"],j["addConnections"][i]["alpha"],i)));
            events.push_back(std::move(
                    std::make_unique<RemoveConnection>(generator, j["removeConnections"][i-1]["p"],j["removeConnections"][i-1]["alpha"],i)));
        }
        events.push_back(std::move(
                std::make_unique<RemoveConnection>(generator, j["removeConnections"][maxConnections-1]["p"],j["removeConnections"][maxConnections-1]["alpha"],maxConnections)));
    }
    auto res = run_simulation(generator, k_max, t_max, n, step_size, events, K);
    std::string directory = "./";
    std::string resFileName = directory + "cpp_results.json";
    nlohmann::json ja(res);
    std::ofstream myfile;
    myfile.open(resFileName);
    myfile << ja.dump();
    myfile.close();
    std::cout<<"finished"<<std::endl;
    return 0;
}
