//
//  Network.cpp
//  iCPG
//
//  Created by Simon Danner on 22/05/15.
//  Copyright (c) 2015 Simon Danner. All rights reserved.
//


#include <boost/lexical_cast.hpp>
#include "Network.hpp"

std::istream& safeGetline(std::istream& is, std::string& t)
{
    //copied from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
    t.clear();
    
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();
    
    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

template<typename T> bool isValid(const std::string& num) {
   bool flag = true; 
   try { 
      T tmp = boost::lexical_cast<T>(num); 
   } 
   catch (boost::bad_lexical_cast &e) { 
      flag = false; 
   } 
   return flag; 
} 



bool allEqual(myvec v){
    int diff=0;
    for (auto it = v.begin();it!=v.end();++it){
        if(v[0]!=(*it)){
            diff++;
        }
    }
    if (diff==0){
        return true;
    }else{
        return false;
    }
}

Network::Network(){
    names = std::map<int,std::string>();
    std::cout << "called Network()" << std::endl;
}

Network::Network(int N_NP, int N_normal):Network(){
    //Network();
    initialize(N_NP,N_normal);
    
    
};

void Network::initialize(int N_NP, int N_normal){
    N_NaP=N_NP;
    N_norm=N_normal;
    
    N_Neurons =N_NaP+N_norm;
    
    ELeak = create_scalarV(N_NaP+N_norm,c_ELeak);
    ENa = create_scalarV(N_NaP+N_norm,c_ENa);
    ESynE = create_scalarV(N_NaP+N_norm,c_ESynE);
    ESynI = create_scalarV(N_NaP+N_norm,c_ESynI);
    
    
    Cmem = create_scalarV(N_NaP+N_norm,c_Cmem);
    
    mk = create_scalarV(N_NaP+N_norm,c_mk);
    mV12 = create_scalarV(N_NaP+N_norm,c_mV12);
    
    hk = create_scalarV(N_NaP+N_norm,c_hk);
    hV12 = create_scalarV(N_NaP+N_norm,c_hV12);
    htau = create_scalarV(N_NaP+N_norm,c_htau);
    hTauK = create_scalarV(N_NaP+N_norm,c_hTauK);
    hTauV12 = create_scalarV(N_NaP+N_norm,c_hTauV12);
    hTau0 = create_scalarV(N_NaP+N_norm,c_hTau0);
    
    gBarLeak = create_scalarV(N_NaP+N_norm,c_gBarLeak);
    gBarNaP = create_scalarV(N_NaP+N_norm,c_gBarNaP);
    
    outputFunction = std::vector<int>(N_NaP+N_norm,0);
    
    Vmax = create_scalarV(N_NaP+N_norm,c_Vmax);
    Vmin = create_scalarV(N_NaP+N_norm,c_Vmin);
    
    tauOutput = create_scalarV(N_Neurons,c_tauOutput);
    
    sigmaNoise = create_scalarV(N_Neurons,c_sigmaNoise);
    tauNoise = create_scalarV(N_Neurons,c_tauNoise);
    
    connE = connection_matrix(N_NaP+N_norm);
    connI = connection_matrix(N_NaP+N_norm);
    connC = connection_matrix(N_Neurons);
    
    driveE = std::list<drive>();
    driveI = std::list<drive>();
    
    feedbackIa = std::list<feedback>();
    feedbackIb = std::list<feedback>();
    feedbackII = std::list<feedback>();
    
    
    transV = myvec(N_Neurons);
    //randomGen = OrnsteinUhlenbeck(N_NP+N_normal,0.0,1.0,5.0);
}

myvec Network::create_scalarV(int N, double k){
    myvec ret = myvec(N);
    std::fill(ret.begin(),ret.end(),k);
    return ret;
}

myvec Network::genInitialCond() const{
    myvec ret;
    if(initial.size()==N_norm+N_NaP*2){
        ret= initial;
    }else{
        ret = myvec(N_norm+N_NaP*2);
        for (int i=0;i<N_NaP+N_norm;i++){
            ret[i]=ELeak[i];
        }
        for( int i=0;i<N_NaP;++i){
            ret[N_Neurons+i]=(int((i+1)/2)%2)*0.3+0.3;
        }
    }
    return ret;
}


void Network::setConnE(const int from,const int to, double *w){
    setConnE(from,to,w,new double(0.0));
}
void Network::setConnI(const int from,const int to, double *w){
    setConnI(from,to,w,new double(0.0));
}
void Network::setConnE(const int from,const int to, double *w, double *slope){
    connection con;
    con.from=from;
    con.weight=w;
    con.slope=slope;
    connE(to).push_back(con);
}
void Network::setConnI(const int from,const int to, double *w,  double *slope){
    connection con;
    con.from=from;
    con.weight=w;
    con.slope=slope;
    connI(to).push_back(con);
}
void Network::setConnC(const int from,const int to,double* w){
    connection con;
    con.from=from;
    con.weight=w;
    connC(to).push_back(con);
}

void Network::resetConn(connection_matrix& cm, const int from, const int to,  double *w){
    for(auto it = cm(to).begin();it!=cm(to).end();++it){
        if( it->from == from){
            it->weight=w;
        }
    }
}


double* Network::getConn(connection_matrix& cm, const int from, const int to){
    for(auto it = cm(to).begin();it!=cm(to).end();++it){
        if( it->from== from){
            return it->weight;
        }
    }
    return new double(-1.0);
}

double* Network::getConnI(const int from, const int to){
    return getConn(connI,from,to);
}
double* Network::getConnE(const int from, const int to){
    return getConn(connE,from,to);
}

// set excitatory drive to neuron to with weight weight and offset offset
void Network::setDriveE(const int to,  double *weight,  double *offset){
    driveE.push_back(drive(to,weight,offset));
}

// set inhibitors drive to neuron to with weight weight and offset offset
void Network::setDriveI(const int to,  double *weight,  double *offset){
    driveI.push_back(drive(to,weight,offset));
}

void Network::setFeedbackIa(const int to, const int fromleg,  const int frommg,  double* weight){
    feedbackIa.push_back(feedback(to,fromleg,frommg,weight));
}
void Network::setFeedbackIb(const int to, const int fromleg, const int frommg,  double* weight){
    feedbackIb.push_back(feedback(to,fromleg,frommg, weight));
}
void Network::setFeedbackII(const int to, const int fromleg, const int frommg,  double* weight){
    feedbackII.push_back(feedback(to,fromleg,frommg,weight));
}
void Network::setFeedbackCutaneous(const int to, const int fromleg,  double* weight){
    feedbackCutaneous.push_back(feedback(to,fromleg,-1,weight));
}




std::string Network::getName(const int neuronID) const{
    return (*names.find(neuronID)).second;
}

int Network::getIndex(const std::string name) const{
    int ret = -1;
    for (auto it=names.begin(); it!=names.end(); ++it){
        if(it->second==name){
            ret = it->first;
        }
    }
    return ret;
}

std::ostream& operator<<(std::ostream& stream, const Network& net){
    stream << "N_NaP " << net.N_NaP << std::endl;
    stream << "N_Normal " << net.N_norm << std::endl;
    stream << std::endl;
    
    stream << "simDuration " << net.simDuration << std::endl;
    //stream << "scalingFactor " << net.scalingFactor << std::endl;
    
    stream << std::endl;
    
    for (int i=0;i<net.N_Neurons;++i){
        stream << "neuron " << i << ": "<<  net.getName(i) << std::endl;
    }
    stream << std::endl;
    
    for(auto it = net.variableMap.begin();it!=net.variableMap.end();++it){
        stream << "variable " << it->first << " " << std::to_string(static_cast<long long>(*(it->second))) << std::endl;
    }
    

    //stream << "connections:" << std::endl;
    for (int to = 0;to<net.N_Neurons;++to){
        for (auto it=net.connE(to).begin();it!=net.connE(to).end();++it){
            if(*it->weight!=0.0){
                stream << "connectionE " <<  net.getName(it->from) << " -> " << net.getName(to) << " : " << *(it->weight) << " + " << *(it->slope) << " * t" << std::endl;
            }
        }
        for (auto it=net.connI(to).begin();it!=net.connI(to).end();++it){
            if(*it->weight!=0.0){
                stream << "connectionI " << net.getName(it->from) << " -o " << net.getName(to) << " : -" << *(it->weight) << " + " << *(it->slope) << " * t" << std::endl;
            }
        }
    }
    stream << std::endl ;
    for (auto it= net.driveE.begin();it!=net.driveE.end();++it){
        if (*it->weight!=0||*it->offset!=0.0){
            stream << "driveE " <<  *it->weight << " * t + " << *it->offset << " -> "  << net.getName(it->to) << std::endl;
        }
    }
    
    for (auto it= net.driveI.begin();it!=net.driveI.end();++it){
        if (*it->weight!=0||*it->offset!=0.0){
            stream << "driveI " << *it->weight << " * t + " << *it->offset << " -o -"  << net.getName(it->to) << std::endl;
        }
    }
    
    for (auto it= net.feedbackIa.begin();it!=net.feedbackIa.end();++it){
        if (it->weight!=0){
            std::string name;
            switch (it->fromleg){
                case 0:
                    name="L_hind";
                    break;
                case 1:
                    name="R_hind";
                    break;
                case 2:
                    name="L_front";
                    break;
                case 3:
                    name="R_front";
                    break;
            }
            stream << "feedbackIa " << name << " -> " << net.getName(it->to) << " : " << it->weight << std::endl;
        }
    }
    for (auto it= net.feedbackIb.begin();it!=net.feedbackIb.end();++it){
        if (it->weight!=0){
            std::string name;
            switch (it->fromleg){
                case 0:
                    name="L_hind";
                    break;
                case 1:
                    name="R_hind";
                    break;
                case 2:
                    name="L_front";
                    break;
                case 3:
                    name="R_front";
                    break;
            }
            stream << "feedbackIb " << name << " -> " << net.getName(it->to) << " : " << it->weight << std::endl;
        }
    }
    for (auto it= net.feedbackII.begin();it!=net.feedbackII.end();++it){
        if (it->weight!=0){
            std::string name;
            switch (it->fromleg){
                case 0:
                    name="L_hind";
                    break;
                case 1:
                    name="R_hind";
                    break;
                case 2:
                    name="L_front";
                    break;
                case 3:
                    name="R_front";
                    break;
            }
            stream << "feedbackII " << name << " -> " << net.getName(it->to) << " : " << it->weight << std::endl;
        }
    }
    
    stream << std::endl;
    stream << "gLeak \t" << (allEqual(net.gBarLeak) ? myvec(1,net.gBarLeak[0]): net.gBarLeak ) << std::endl;
    stream << "gBarNaP\t " << (allEqual(net.gBarNaP) ? myvec(1,net.gBarNaP[0]): net.gBarNaP ) << std::endl;
    stream << "Eleak\t " << (allEqual(net.ELeak) ? myvec(1,net.ELeak[0]): net.ELeak ) << std::endl;
    stream << "ENa\t " << (allEqual(net.ENa) ? myvec(1,net.ENa[0]): net.ENa ) << std::endl;
    stream << "ESynE\t " << (allEqual(net.ESynE) ? myvec(1,net.ESynE[0]): net.ESynE ) << std::endl;
    stream << "ESynI\t " << (allEqual(net.ESynI) ? myvec(1,net.ESynI[0]): net.ESynI ) << std::endl;
    stream << "Cmem\t " << (allEqual(net.Cmem) ? myvec(1,net.Cmem[0]): net.Cmem ) << std::endl;
    stream << "mk\t " << (allEqual(net.mk) ? myvec(1,net.mk[0]): net.mk ) << std::endl;
    stream << "mV12\t " << (allEqual(net.mV12) ? myvec(1,net.mV12[0]): net.mV12 ) << std::endl;
    
    stream << "hk\t " << (allEqual(net.hk) ? myvec(1,net.hk[0]): net.hk ) << std::endl;
    stream << "hV12\t " << (allEqual(net.hV12) ? myvec(1,net.hV12[0]): net.hV12 ) << std::endl;
    
    stream << "htau\t " << (allEqual(net.htau) ? myvec(1,net.htau[0]): net.htau ) << std::endl;
    stream << "hTauK\t " << (allEqual(net.hTauK) ? myvec(1,net.hTauK[0]): net.hTauK ) << std::endl;
    stream << "hTauV12\t " << (allEqual(net.hTauV12) ? myvec(1,net.hTauV12[0]): net.hTauV12 ) << std::endl;
    
    stream << "hTau0\t " << (allEqual(net.hTau0) ? myvec(1,net.hTau0[0]): net.hTau0 ) << std::endl;
    
    stream << "Vmax\t " << (allEqual(net.Vmax) ? myvec(1,net.Vmax[0]): net.Vmax ) << std::endl;
    stream << "Vmin\t " << (allEqual(net.Vmin) ? myvec(1,net.Vmin[0]): net.Vmin ) << std::endl;
    stream << "sigmaNoise\t " << (allEqual(net.sigmaNoise) ? myvec(1,net.sigmaNoise[0]): net.sigmaNoise ) << std::endl;
    stream << "tauNoise\t " << (allEqual(net.tauNoise) ? myvec(1,net.tauNoise[0]): net.tauNoise ) << std::endl;
    stream << "outputFunction\t " << net.outputFunction  << std::endl;
    stream << "initialConditions\t" << net.genInitialCond() << std::endl;
    
    return stream;
}

void Network::set_para(myvec &to,double value,int start,int end){
    if((start<=end)&&(end<to.size())){
        for (int i = 0;i<to.size();i++){
            if (i>=start&&i<=end){
                to[i]=value;
            }else{
                to[i]=-1234567890;
            }
        }
    }else{
        std::cout << "error in index for parameter" << start << " " << end <<std::endl;
        return;
    }
}
void Network::assign_para(myvec &to,myvec from){
    if(to.size()==from.size()){
        for(int i=0;i<to.size();i++){
            if(from[i]!=-1234567890){
                to[i]=from[i];
            }
        }
    }
}

Network::Network(std::string filename,std::vector<std::string> musclenames, std::vector<std::vector<std::string>> mnnames){
    std::ifstream myfile (filename);
    std::string line;
    
    if (myfile.is_open())
    {
        int NNaP=-1;
        int NNorm=-1;
        for (int i = 0;i<=1;++i){
            std::string::size_type sz;
            getline(myfile,line);
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of("\t "));
            
            if (strs[0]=="N_NaP"){
                NNaP = std::stoi(strs[1],&sz);
                std::cout << "NNaP = " << NNaP << std::endl;
            }else if(strs[0]=="N_Normal"){
                NNorm = std::stoi(strs[1],&sz);
                std::cout << "NNorm = " << NNorm << std::endl;
            }
        }
        if(NNaP==-1||NNorm==-1){
            std::cout << "N_NaP and/or N_Normal not defined in the first two lines" << std::endl;
            return;
        }
        initialize(NNaP,NNorm);
        lscond=std::vector<LimbSensorCondition>(2);
        for(int i = 0;i<musclenames.size();i++){
            for(int k = 0;k<2;k++){
                lscond[k].Ia.push_back(0.0);
                lscond[k].Ib.push_back(0.0);
                lscond[k].II.push_back(0.0);
            }
        }
        myvec ofun_double=create_scalarV(N_NaP+N_norm,0.0);
        int no_neuron = -1;
        while ( safeGetline (myfile,line) )
        {
            
            std::string::size_type sz;
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of("\t "));
            
            std::vector<std::string>::iterator i = strs.begin();
            while(i != strs.end())
            {
                if((*i)==std::string("")){
                    strs.erase(i);
                }else{
                    ++i;
                }
            }
            if(strs.size()<2){
                continue;
            }
            if (strs[0]=="neuron"){
                //int no = std::stoi(strs[1],&sz);
                if(++no_neuron<N_Neurons){
                    std::string neuronname;
                    if(strs.size()==2){
                        neuronname=strs[1];
                    }else{
                        neuronname=strs[2];
                    }
                    names[no_neuron]=neuronname;
                    std::cout << "adding neuron " << neuronname << " nr "  << no_neuron   << std::endl;
                }else{
                    std::cout << "max nr of neurons exceeded with " << no_neuron << "neurons"   << std::endl;
                }
            }else if (strs[0]=="variable"){
                double* value = new double;
                *value = std::stod(strs[2],&sz);
                variableMap[strs[1]]=value;
            }else if (strs[0]=="connectionE"||strs[0]=="connectionI"||strs[0]=="connectionC"){
                int from = -1;
                int to = -1;
                double* weight;
                double* slope;
                for (auto it=names.begin(); it!=names.end(); ++it){
                    if(it->second==strs[1]){
                        from = it->first;
                    }
                    if(it->second==strs[3]){
                        to = it->first;
                    }
                }
                if(isValid<double>(strs[5])){
                    weight = new double;
                    *weight = std::stod(strs[5],&sz);
                }else if(variableMap.find(strs[5])!=variableMap.end()){
                    auto it = variableMap.find(strs[5]);
                    weight = it->second;
                }
                if(isValid<double>(strs[7])){
                    slope = new double;
                    *slope = std::stod(strs[7],&sz);
                }else if(variableMap.find(strs[7])!=variableMap.end()){
                    auto it = variableMap.find(strs[7]);
                    slope = it->second;
                }
                if(from==-1||to==-1){
                    std::cout << "at line " << line << " neuron names were not recognized" << std::endl;
                }else{
                    if (strs[0]=="connectionE"){
                        setConnE(from,to,weight,slope);
                        std::cout << "adding excitatory connection from " << from << strs[1] << " to " << to << strs[3] << " w " << *weight << " slope " << *slope << std::endl;
                    }
                    if (strs[0]=="connectionI"){
                        setConnI(from,to,weight,slope);
                        std::cout << "adding inhibitory connection from " << from << strs[1] << " to " << to << strs[3] << " w " << *weight <<  " slope " << *slope << std::endl;
                    }
                    if (strs[0]=="connectionC"){
                        setConnC(from,to,weight);
                        std::cout << "adding cholinergi connection from " << from << strs[1] << " to " << to << strs[3] << " w " << *weight << std::endl;
                    }
                }
            }else if(strs[0]=="driveE"||strs[0]=="driveI"){
                double *weight = nullptr;
                double *offset = nullptr;
                int to=-1;
                for (auto it=names.begin(); it!=names.end(); ++it){
                    if(it->second==strs[7]){
                        to = it->first;
                    }
                }
                if(isValid<double>(strs[1])){
                    weight = new double;
                    *weight = std::stod(strs[1],&sz);
                }else if(variableMap.find(strs[1])!=variableMap.end()){
                    auto it = variableMap.find(strs[1]);
                    weight = it->second;
                }
                if(isValid<double>(strs[5])){
                    offset = new double;
                    *offset = std::stod(strs[5],&sz);
                }else if(variableMap.find(strs[5])!=variableMap.end()){
                    auto it = variableMap.find(strs[5]);
                    offset = it->second;
                }
                if(strs[0]=="driveE"){
                    setDriveE(to, weight, offset);
                    std::cout << "adding excitatory drive to " << to << strs[7] << " weight " << *weight << " offset " << *offset << std::endl;
                }else if (strs[0]=="driveI"){
                    setDriveI(to, weight, offset);
                    std::cout << "adding inhibitory drive to " << to << strs[7] << " weight " << *weight << " offset " << *offset << std::endl;
                }
            }else if (strs[0]=="feedbackIa"||strs[0]=="feedbackIb"||strs[0]=="feedbackII"||strs[0]=="feedbackCutaneous"){
                double* weight;
                if(isValid<double>(strs[6])){
                    weight = new double;
                    *weight = std::stod(strs[6],&sz);
                }else if(variableMap.find(strs[6])!=variableMap.end()){
                    auto it = variableMap.find(strs[6]);
                    weight = it->second;
                }
                int to=-1;
                int fromleg=-1;
                int frommg=-1;
                for (auto it=names.begin(); it!=names.end(); ++it){
                    if(it->second==strs[4]){
                        to = it->first;
                    }
                }
                if(strs[1]=="L_hind"){
                    fromleg=0;
                }
                if(strs[1]=="R_hind"){
                    fromleg=1;
                }
                if(strs[1]=="L_front"){
                    fromleg=2;
                }
                if(strs[1]=="R_front"){
                    fromleg=3;
                }
                for(std::size_t m_in=0;m_in<musclenames.size();++m_in){
                    if(strs[2]==musclenames[m_in]){
                        frommg=static_cast<int>(m_in);
                    }
                }
                
                if(strs[0]=="feedbackIa"){
                    setFeedbackIa(to, fromleg, frommg, weight);
                    std::cout << "adding feedbackIa from "<< fromleg << strs[1] << " " << strs[2] << frommg << " to " << to << strs[3] << " weight " << *weight << std::endl;
                }else if(strs[0]=="feedbackIb"){
                    setFeedbackIb(to, fromleg,frommg, weight);
                    std::cout << "adding feedbackIb from "<< fromleg << strs[1] << " " << strs[2] << frommg << " to " << to << strs[3] << " weight " << *weight << std::endl;
                }else if(strs[0]=="feedbackII"){
                    setFeedbackII(to, fromleg,frommg, weight);
                    std::cout << "adding feedbackII from "<< fromleg << strs[1] << " " << strs[2] << frommg << " to " << to << strs[3] << " weight " << *weight << std::endl;
                }else if(strs[0]=="feedbackCutaneous"){
                    setFeedbackCutaneous(to, fromleg, weight);
                    std::cout << "adding feedbackCutaneous from "<< fromleg << strs[1] << " to " << to << strs[3] << " weight " << *weight << std::endl;
                }
            }else if (strs[0]=="feedbackBodyTilt") {
                feedback_body_tilt bt;
                
                if(isValid<double>(strs[6])){
                    bt.weight = new double;
                    *(bt.weight) = std::stod(strs[6],&sz);
                }else if(variableMap.find(strs[6])!=variableMap.end()){
                    auto it = variableMap.find(strs[6]);
                    bt.weight = it->second;
                }
                if(strs.size()>8){
                    if(isValid<double>(strs[8])){
                        bt.cutoff = new double;
                        *(bt.cutoff) = std::stod(strs[8],&sz);
                    }else if(variableMap.find(strs[8])!=variableMap.end()){
                        auto it = variableMap.find(strs[8]);
                        bt.cutoff = it->second;
                    }
                }else{
                    bt.cutoff = new double(0.0);
                }
                for (auto it=names.begin(); it!=names.end(); ++it){
                    if(it->second==strs[4]){
                        bt.to = it->first;
                    }
                }
                if(strs[1]=="ANTERIOR"){
                    bt.direction=0;
                }
                if(strs[1]=="POSTERIOR"){
                    bt.direction=1;
                }
                if(strs[1]=="LEFT"){
                    bt.direction=2;
                }
                if(strs[1]=="RIGHT"){
                    bt.direction=3;
                }
                if(strs[2]=="ANGLE"){
                    bt.type=0;
                }
                if(strs[2]=="VELOCITY"){
                    bt.type=1;
                }
                feedbackBodyTilt.push_back(bt);
                std::cout << "adding feedbackBodyTilt " <<  strs[1] << bt.direction << " " << strs[2] << bt.type << " to " << strs[4] << bt.to << " with weight " << *bt.weight << " cutoff " << *bt.cutoff << std::endl;
            }else if (strs[0]=="simDuration") {
                simDuration=std::stod(strs[1],&sz);
                std::cout << "setting simDuration to " << simDuration << std::endl;
            }else if (strs[0]=="scalingFactor"){
                sf=std::stod(strs[1],&sz);
                std::cout << "setting scalingFactor to " << sf << std::endl;
            }else if (strs[0]=="stepwise"){
                stepwise=true;
                nSteps = std::stoi(strs[1]);
            }else if(strs[0][0]=='/'){
                // do nothing
                
            }else if(strs[0]=="initialConditions"){
                initial =myvec(N_norm+N_NaP*2);
                for (int j = 1;j<=N_norm+N_NaP*2;++j){
                    if (j>=strs.size()){
                        std::cout<< "not all initial conditions specified" << std::endl;
                        break;
                    }
                    initial[j-1]=std::stod(strs[j],&sz);
                }
            }else{
                myvec para;
                if (strs.size()==2){
                    para = create_scalarV(N_NaP+N_norm,std::stod(strs[1],&sz));
                }else if (strs.size()==3){
                    std::vector<std::string> substrs;
                    boost::split(substrs, strs[1], boost::is_any_of("[]:"));
                    auto iter = substrs.begin();
                    while(iter != substrs.end())
                    {
                        if((*iter)==std::string("")){
                            substrs.erase(iter);
                        }else{
                            ++iter;
                        }
                    }
                    if(substrs.size()==2){
                        para = myvec(N_Neurons);
                        set_para(para,std::stod(strs[2],&sz),std::stod(substrs[0],&sz),std::stod(substrs[1],&sz));
                        std::cout<< "update " << strs[0] << " for neurons " << std::stod(substrs[0],&sz) << " : "
                        << std::stod(substrs[1],&sz) << " to " << std::stod(strs[2],&sz) << std::endl;
                    }else{
                        std::cout<< "syntax error in specification for " << strs[0] << std::endl;
                        break;
                    }
                }else if (strs.size()==4){
                    if(strs[1]=="neurons"){
                        para = myvec(N_Neurons,-1234567890);
                        std::string to;
                        for (auto it=names.begin(); it!=names.end(); ++it){
                            if(it->second.find(strs[2])!=std::string::npos){
                                para[it->first]=std::stod(strs[3],&sz);
                                to+=std::to_string(static_cast<long long>(it->first))+", ";
                            }
                        }
                        std::cout << "update " << strs[0] << " for neurons " << to << "(*" << strs[2] << "*) to "
                        << std::stod(strs[3],&sz) << std::endl;
                    }else{
                        std::cout<< "syntax error in specification for " << strs[0] << std::endl;
                        break;
                    }
                }else {
                    para = myvec(N_Neurons);
                    for (int j = 1;j<=N_NaP+N_norm;j++){
                        if (j>strs.size()){
                            std::cout<< "not enough values specified for parameter " << strs[0] << std::endl;
                            break;
                        }else{
                            para(j-1)=std::stod(strs[j],&sz);
                        }
                    }
                }
                    if (strs[0]=="gLeak"){
                        assign_para(gBarLeak,para);
                        std::cout << "set gLeak to " << (allEqual(gBarLeak) ? myvec(1,gBarLeak[0]): gBarLeak) << std::endl;
                    }else if(strs[0] == "gBarNaP"){
                        assign_para(gBarNaP,para);
                        std::cout << "set gBarNaP to " << (allEqual(gBarNaP) ? myvec(1,gBarNaP[0]): gBarNaP) << std::endl;
                    }else if(strs[0] == "Eleak"){
                        assign_para(ELeak, para);
                        std::cout << "set Eleak to " << (allEqual(ELeak) ? myvec(1,ELeak[0]): ELeak) << std::endl;
                    }else if(strs[0] == "ENa"){
                        assign_para(ENa,para);
                        std::cout << "set ENa to " << (allEqual(ENa) ? myvec(1,ENa[0]): ENa) << std::endl;
                    }else if(strs[0] == "ESynE"){
                        assign_para(ESynE,para);
                        std::cout << "set ESynE to " << (allEqual(ESynE) ? myvec(1,ESynE[0]): ESynE) << std::endl;
                    }else if(strs[0] == "ESynI"){
                        assign_para(ESynI,para);
                        std::cout << "set ESynI to " << (allEqual(ESynI) ? myvec(1,ESynI[0]): ESynI) << std::endl;
                    }else if(strs[0] == "Cmem"){
                        assign_para(Cmem,para);
                        std::cout << "set Cmem to " << (allEqual(Cmem) ? myvec(1,Cmem[0]): Cmem) << std::endl;
                    }else if(strs[0] == "mk"){
                        assign_para(mk,para);
                        std::cout << "set mk to " << (allEqual(mk) ? myvec(1,mk[0]): mk) << std::endl;
                    }else if(strs[0] == "mV12"){
                        assign_para(mV12,para);
                        std::cout << "set mV12 to " << (allEqual(mV12) ? myvec(1,mV12[0]): mV12) << std::endl;
                    }else if(strs[0] == "hk"){
                        assign_para(hk,para);
                        std::cout << "set hk to " << (allEqual(hk) ? myvec(1,hk[0]): hk) << std::endl;
                    }else if(strs[0] == "hV12"){
                        assign_para(hV12,para);
                        std::cout << "set hV12 to " << (allEqual(hV12) ? myvec(1,hV12[0]): hV12) << std::endl;
                    }else if(strs[0] == "htau"){
                        assign_para(htau,para);
                        std::cout << "set htau to " << (allEqual(htau) ? myvec(1,htau[0]): htau) << std::endl;
                    }else if(strs[0] == "hTauK"){
                        assign_para(hTauK,para);
                        std::cout << "set hTauK to " << (allEqual(hTauK) ? myvec(1,hTauK[0]): hTauK) << std::endl;
                    }else if(strs[0] == "hTau0"){
                        assign_para(hTau0,para);
                        std::cout << "set hTau0 to " << (allEqual(hTau0) ? myvec(1,hTau0[0]): hTau0) << std::endl;
                    }else if(strs[0] == "Vmax"){
                        assign_para(Vmax,para);
                        std::cout << "set Vmax to " << (allEqual(Vmax) ? myvec(1,Vmax[0]): Vmax) << std::endl;
                    }else if(strs[0] == "Vmin"){
                        assign_para(Vmin,para);
                        std::cout << "set Vmin to " << (allEqual(Vmin) ? myvec(1,Vmin[0]): Vmin) << std::endl;
                    }else if(strs[0] == "hTauV12"){
                        assign_para(hTauV12,para);
                        std::cout << "set hTauV12 to " << (allEqual(hTauV12) ? myvec(1,hTauV12[0]): hTauV12) << std::endl;
                    }else if(strs[0] == "sigmaNoise"){
                        assign_para(sigmaNoise,para);
                        std::cout << "set sigmaNoise to " << (allEqual(sigmaNoise) ? myvec(1,sigmaNoise[0]): sigmaNoise) << std::endl;
                    }else if(strs[0] == "tauNoise"){
                        assign_para(tauNoise,para);
                        std::cout << "set tauNoise to " << (allEqual(tauNoise) ? myvec(1,tauNoise[0]): tauNoise) << std::endl;
                    }else if(strs[0] == "outputFunction"){
                        assign_para(ofun_double,para);
                        std::cout << "set linearOutput to " << (allEqual(ofun_double) ? myvec(1,ofun_double[0]): ofun_double) << std::endl;
                    }
                
            }
        }
        while(no_neuron<N_Neurons-1){
            names[++no_neuron]=std::string("not initialized");
        }
        myfile.close();
        randomGen = OrnsteinUhlenbeck(NNaP+NNorm,0.0,1.0,tauNoise[0]);
        for(int in = 0; in<outputFunction.size();in++){
            outputFunction[in]=std::lround(ofun_double[in]);
        }
    }
    
    
    for (int i = 0;i<mnnames.size();++i){
        std::vector<int> v_;
        for(int j = 0 ; j<mnnames[0].size();++j){
            v_.push_back(getIndex(mnnames[i][j]));
            std::cout << mnnames[i][j] << " is neuron nr " << getIndex(mnnames[i][j]) << std::endl;
        }
        in_Act.push_back(v_);
    }
}

inline double sech(double z){return 2/(exp(z)+exp(-z));};
//inline double pos(double d){return d>0.0?d:0.0;};
inline double pos(double d){return std::signbit(d)?0.0:d;};

//integration step
void Network::step(const myvec &x, myvec &dxdt, double t){
    myvec randV = randomGen.get(t);
    for (int i = 0;i<N_Neurons;++i){
        switch(outputFunction[i]){
            case 0:
                transV(i) = std::min(1.0,pos(x[i]-Vmin[i])/(Vmax[i]-Vmin[i]));
                break;
            case 1:
                transV(i) = std::min(1.0,pos(pow(x[i]-Vmin[i],2)/pow(Vmax[i]-Vmin[i],2)));
                break;
        }
    }
    for (int i = 0;i<N_Neurons;++i){
        // cholinergic connections need to be listed in order and can't contain circles
        for (auto it=connC(i).begin();it!=connC(i).end();++it){
            if((1.0 + (*it->weight)*transV(it->from))>1.0){
                transV[i] *= (1.0 + (*it->weight)*transV(it->from));
            }
        }
    }
    for (int i = 0;i<N_Neurons;++i){
        dxdt[i]=(-gBarLeak[i]*(x[i]-ELeak[i]));
        for (auto it=connE(i).begin();it!=connE(i).end();++it){
            dxdt[i]-=pos(*it->weight)*transV(it->from)*(x[i]-ESynE[i]);
        }
        for (auto it=connI(i).begin();it!=connI(i).end();++it){
            dxdt[i]-=pos(*it->weight)*transV(it->from)*(x[i]-ESynI[i]);
        }
        dxdt[i]+=-(randV[i]*sigmaNoise[i]);
        dxdt[i]/=Cmem[i];
        if (i < N_NaP){
            double mpInf = 1./(1.+exp((x[i]-mV12[i])/mk[i]));
            double hpInf = 1./(1.+exp((x[i]-hV12[i])/hk[i]));
            double tau_inf = hTau0[i]+(htau[i]-hTau0[i])/(cosh((x[i]-hTauV12[i])/hTauK[i]));
            
            dxdt[i]-=gBarNaP[i]*x[N_NaP+N_norm+i]*mpInf*(x[i]-ENa[i])/Cmem[i];
            
            dxdt[N_NaP+N_norm+i]=(hpInf-x[N_NaP+N_norm+i])/tau_inf;
        }
    }
    
    for(auto it=feedbackIa.begin();it!=feedbackIa.end();++it){
        if(*it->weight>=0){
            dxdt[it->to]-=lscond[it->fromleg].Ia[it->frommg]*(*it->weight)*(x[it->to]-ESynE[it->to])/Cmem[it->to];
        }else{
            dxdt[it->to]-=lscond[it->fromleg].Ia[it->frommg]*-1.0*(*it->weight)*(x[it->to]-ESynI[it->to])/Cmem[it->to];
        }
    }
    for(auto it=feedbackIb.begin();it!=feedbackIb.end();++it){
        if(*it->weight>=0){
            dxdt[it->to]-=lscond[it->fromleg].Ib[it->frommg]*(*it->weight)*(x[it->to]-ESynE[it->to])/Cmem[it->to];
        }else{
            dxdt[it->to]-=lscond[it->fromleg].Ib[it->frommg]*-1.0*(*it->weight)*(x[it->to]-ESynI[it->to])/Cmem[it->to];
        }
    }
    for(auto it=feedbackII.begin();it!=feedbackII.end();++it){
        if(*it->weight>=0){
            dxdt[it->to]-=lscond[it->fromleg].II[it->frommg]*(*it->weight)*(x[it->to]-ESynE[it->to])/Cmem[it->to];
        }else{
            dxdt[it->to]-=lscond[it->fromleg].II[it->frommg]*-1.0*(*it->weight)*(x[it->to]-ESynI[it->to])/Cmem[it->to];
        }
    }
    for(auto it=feedbackCutaneous.begin();it!=feedbackCutaneous.end();++it){
        if(*it->weight>=0){
            dxdt[it->to]-=lscond[it->fromleg].cutaneous*(*it->weight)*(x[it->to]-ESynE[it->to])/Cmem[it->to];
        }
    }
    for(auto it=feedbackBodyTilt.begin();it!=feedbackBodyTilt.end();++it){
        double angle;
        double velocity;
        if (it->direction==0||it->direction==1){
            angle = body_tilt.anterior_posterior_angle;
            velocity = body_tilt.anterior_posterior_velocity;
        }else{
            angle = body_tilt.left_right_angle;
            velocity = body_tilt.left_right_velocity;
        }
        double value;
        if(it->type==0){
            value = angle;
        }else{
            value = velocity;
        }
        value-=*it->cutoff;
        if (it->direction==1||it->direction==3){
            value *= -1.0;
        }
        if(value>0.0 && *it->weight){
            //std::cout << "adding tiltfb" <<std::endl;
            dxdt[it->to]-=value*(*it->weight)*(x[it->to]-ESynE[it->to])/Cmem[it->to];
        }
    }
    for (auto it = driveE.begin(); it != driveE.end(); ++it){
        dxdt[it->to] -= ((alpha * *it->weight * 1e5) + *it->offset) * (x[it->to] - ESynE[it->to]) / Cmem[it->to];
    }

    for (auto it = driveI.begin(); it != driveI.end(); ++it){
        dxdt[it->to] -= ((alpha * *it->weight * 1e5) + *it->offset) * (x[it->to] - ESynI[it->to]) / Cmem[it->to];
    }
}
OrnsteinUhlenbeck::OrnsteinUhlenbeck(int n, double m, double s, double t):N(n),mu(m),sigma(s),tau(t){
    typedef std::chrono::high_resolution_clock myclock;
    seed=(unsigned)myclock::now().time_since_epoch().count();
    generator.seed(seed);
    data = boost::circular_buffer<myvec>(Nhist);
    for(int i = 0;i<Nhist;++i)
        data.push_back(myvec(N));
    distribution=std::normal_distribution<double>(0.0,1.0);
};

void OrnsteinUhlenbeck::calculateUntil(double t){
    for(int i = 1; i<=t+1;++i){
        myvec x=data[i-1];
        data[i]=myvec(N);
        for(int j=0;j<N;++j){
            (data[i])[j]=x[j]+(mu-x[j])/tau+sigma*(std::sqrt((2./tau)))*distribution(generator);
        }
    }
}
myvec OrnsteinUhlenbeck::get(double t){
    if(std::ceil(time0-t)+1>Nhist){
        time0 = t-2.0; // reset buffer in case too old time points are requested (usually happens only when simulation is reset)
    }
    while(std::ceil(t)>=time0){
        myvec temp(N);
        for(int j=0;j<N;++j){
            double x = data.back()[j];
            temp[j]=x+(mu-x)/tau+sigma*(std::sqrt((2./tau)))*distribution(generator);
        }
        data.push_back(temp);
        time0+=1.0;
    }
    int nb = std::ceil(time0-t);
    myvec r1 = *(data.rend()+nb);
    myvec r2 = *(data.rend()+(nb-1));
    myvec ret(N);
    double k = t-std::floor(t);
    for(int i = 0; i<N;i++){
        ret[i] = r1[i] + (r2[i] - r1[i])*k;
    }
    return ret;
}

bool Network::updateVariable(const std::string name,const double value){
    for(auto it = variableMap.begin();it!=variableMap.end();++it){
        if(name==it->first){
            (*it->second)=value;
            std::cout << "set variable " << name << " to " << value << std::endl;
            return true;
        }
    }
    return false;
}

double Network::getVariableValue(std::string name){
    return *variableMap[name];
};

double* Network::getVariablePointers(std::string name){
    return variableMap[name];
};

void Network::setPara(std::string parameter, std::string neuron_substr,double value){
    for (auto it=names.begin(); it!=names.end(); ++it){
        if(it->second.find(neuron_substr)!=std::string::npos){
            if(parameter=="ELeak"){
                ELeak[it->first]=value;
            }
            std::cout << "reset " << parameter << " for neuron " << it->second << " to " << value << std::endl;
        }
    }
}
