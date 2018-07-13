#include <iostream>
#include "CPGNetworkSimulator.hpp"

int main(){
    std::vector<std::string> musclenames = {"IP", "GM", "VL", "TA", "SO", "BF", "GA", "HAM", "RF"};
    std::vector<std::vector<std::string>> mnnames;
    std::vector<std::string> left = {"MnIP_L_hind", "MnGM_L_hind", "MnVL_L_hind", "MnTA_L_hind", "MnSO_L_hind", "MnBF_L_hind","MnGA_L_hind"};
    std::vector<std::string> right = {"MnIP_R_hind", "MnGM_R_hind", "MnVL_R_hind", "MnTA_R_hind", "MnSO_R_hind", "MnBF_R_hind","MnGA_R_hind"};
    mnnames.push_back(left);
    mnnames.push_back(right);
    
    CPGNetworkSimulator simulator("./models/4CPG9MN_20180614_FB_LPN_FB.txt",musclenames,mnnames);
    simulator.setAlpha(0.3);
    double dt = 0.001;
    double dur  = 10.0;

    for (double t=0.0;t<dur;t+=dt){
        simulator.step(dt);
        std::vector<std::vector<double>> act = simulator.getAct();
        std::cout << act[0][0] << " " <<  act[1][0] << std::endl;
    }
    return 0;
}