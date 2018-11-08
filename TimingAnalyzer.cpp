///// c++ Program for a simple Timing Analyzer
/// Author : Ali Abyaneh 810193626

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include<vector>
#include <stdlib.h>
#include <unordered_map>
using namespace std;

double findMax(std::vector<float> &f){
  double max = f[0];
  for (size_t i = 1; i < f.size(); i++) {
    if(f[i] > max)
      max = f[i];
  }
  return max;
}
double findMin(std::vector<float> &f){
  double min = f[0];
  for (size_t i = 1; i < f.size(); i++) {
    if(f[i] < min)
      min = f[i];
  }
  return min;
}
void STOF(std::vector<string> s, std::vector<float> &f){
  for (size_t i = 0; i < s.size(); i++) {
    f.push_back(atof(s[i].c_str()));
    // atof(splitedLine[1].c_str())
  }
}
void split_string(string line, vector <string> &result)
{
    int index = 0;
    result.push_back(string());
    int index_ini = line.find(".") - 1;
    for (int i = index_ini; i < line.size(); i++)
    {
        if(line[i] == ':' || line[i] == ')' || line[i] == ' ' || line[i] == '(')
        {
            index++;
            while(line[i] == ':' || line[i] == ')' || line[i] == ' ' || line[i] == '(')
                i++;
            if(i < line.size() - 1 )
              result.push_back(string());
        }
        if(i < line.size() - 1 )
          result[index].push_back(line[i]);
    }
}

struct DFF{
  double minDelay;
  double maxDelay;
  double Hold;
  double Setup;
};
class DelayData{
public:
  map<std::string, map<std::string, double> > minDelay;
  map<std::string, map<std::string, double> > maxDelay;
  map<std::string, map<std::string, DFF> > DFFdelay;
};
class SDFParser{
  ifstream fp;
  DelayData dData;
public:
  SDFParser(const char* file_name){
    fp.open(file_name);
  }
  void Parse();
  void getDFF(std::string INSTANCE, std::string CELLTYPE);
  bool getCOMB(std::string INSTANCE, std::string CELLTYPE);
  DelayData getData();
};
class Graph{
private:
    unordered_map<string, unordered_map<string, float> > nodes;  // key-vlue indicate Edges
public:
    void add_node(string node_name, const unordered_map<string, float> &edge);
    map<string, float> Diakstra(string source);
    float calcLongestPath(string source, string Dest, float pre_dist, float &max_dist,map<string, bool> &visited);
    float LongestPath(string source, string Dest);
    string calcMinDis(map<string, float> distance, map<string, bool> visited);
};

float max(map<string,float> a){
    float m = 0.0;
    for(auto & i : a){
        if(i.second > m)
            m = i.second;
    }
    return m;
}
float Graph::LongestPath(string source, string Dest){
    float dist = 0;
    map<string, bool> visited;
    for(auto & i : nodes){
        visited[i.first] = false;
    }
    if(source != Dest)
        visited[source] = true;
    bool v = calcLongestPath(source,Dest,0.0,dist,visited);
    return dist;
}
 map<string,float>  AND( map<string,float>a,  map<string,bool>b){
     for(auto & i : a){
         i.second = i.second * (float)b[i.first];
     }
     return a;
}
float Graph::calcLongestPath(string source, string Dest, float pre_dist, float &max_dist, map<string, bool> &visited){
    map<string,float>distance;
    map<string,bool>valid;
    for(auto & i : nodes){
        distance[i.first] = 0.0;
        valid[i.first] = false;
    }
    for(auto & i : nodes){
        if(visited[i.first] == false)
            if(nodes[source].find(i.first) != nodes[source].end())
            {
                visited[i.first] = true;
                distance[i.first] = nodes[source][i.first] + pre_dist;
                valid[i.first] = calcLongestPath(i.first,Dest,distance[i.first], max_dist,visited);

                if(i.first == Dest)
                    valid[i.first] = true;
                if(!valid[i.first])
                    distance[i.first] = pre_dist;
                if(valid[i.first])
                    visited[i.first] = false;

            }
             map<string,float> x = AND(distance,valid);
            if(max(x) > max_dist)
                max_dist = max(distance);
    }
    for(auto & i : nodes){
        if(valid[i.first])
            return 1;
    }
    return 0;

}
string Graph::calcMinDis(map<string, float> distance, map<string, bool> visited){
    float minDistance = 1e10;
    string minNode;
    for(auto & i : visited){
        if(i.second == false && distance[i.first] <= minDistance)
        {
            minNode = i.first;
            minDistance = distance[i.first];
        }
    }
    return minNode;
}
void Graph::add_node(string node_name, const unordered_map<string, float> &edge){
    nodes.insert(unordered_map<string, unordered_map<string, float> >::value_type(node_name, edge));
}
 map<string, float> Graph::Diakstra(string source){
    map<string, float> distance;
    map<string, bool> visited;
    for(auto & i : nodes){
        distance[i.first] = 1e10;
        visited[i.first] = false;
    }
    distance[source] = 0;
     for(auto & i : nodes){
        string s = Graph::calcMinDis(distance, visited);
        visited[s] = true;
        for(auto & j : nodes){
            if(visited[j.first] == false && nodes[s].find(j.first) != nodes[s].end() )
                if(nodes[s].find(j.first)->second + distance[s] < distance[j.first])
                    distance[j.first] =  nodes[s][j.first] + distance[s];
        }
     }
     return distance;
}
DelayData SDFParser::getData(){
    return this->dData;
}
void SDFParser::getDFF(std::string INSTANCE, std::string CELLTYPE){
    std::string line;
    std::getline(fp, line);
    line.clear();
    std::getline(fp, line);
    line.clear();
    std::getline(fp, line);
    std::vector<string> splitedLine;
    std::vector<float> delays;
    split_string(line, splitedLine);
    STOF(splitedLine, delays);
    double maxQ = findMax(delays);
    double minQ = findMin(delays);
    std::getline(fp, line);
    std::vector<string> splitedLine_2;
    std::vector<float> delays_2;
    split_string(line, splitedLine_2);
    STOF(splitedLine_2, delays_2);
    double maxN = findMax(delays_2);
    double minN = findMin(delays_2);
    // cout << splitedLine_2.size() << endl;
    // cout << maxQ << endl;
    std::getline(fp, line);
    std::getline(fp, line);
    std::getline(fp, line);
    vector<double> max(4);
    vector<double> min(4);
    for (size_t i = 0; i < 4; i++) {
      std::getline(fp, line);
      std::vector<string> splitedLine;
      std::vector<float> delays;
      split_string(line, splitedLine);
      STOF(splitedLine, delays);
      max[i] = findMax(delays);
      min[i] = findMin(delays);
    }
    dData.DFFdelay[CELLTYPE][INSTANCE].minDelay = minQ;
    dData.DFFdelay[CELLTYPE][INSTANCE].maxDelay = maxQ;
    dData.DFFdelay[CELLTYPE][INSTANCE].Hold = max[0];
    dData.DFFdelay[CELLTYPE][INSTANCE].Setup = max[2];
//    cout << CELLTYPE << " " << INSTANCE << dData.DFFdelay[CELLTYPE][INSTANCE].minDelay << endl;// = minQ;


}

bool SDFParser::getCOMB(std::string INSTANCE, std::string CELLTYPE){
  std::string line;
  std::getline(fp, line);
  std::getline(fp, line);
  vector<float> max;
  vector<float> min;
  int CELL_Key = 0;
  while(!fp.eof())
  {
    std::getline(fp, line);
    int pos = line.find("CELL");
    if(pos != -1)
      CELL_Key = 1;
    if(CELL_Key == 1 || fp.eof())
      break;
    if(line != ")\r" && line != "  )\r" && line != "    )\r"){
        std::vector<string> splitedLine;
        std::vector<float> delays;
        split_string(line, splitedLine);
        STOF(splitedLine, delays);
        max.push_back(findMax(delays));
        min.push_back(findMin(delays));
    }
  }
  dData.minDelay[CELLTYPE][INSTANCE] = findMin(min);
  dData.maxDelay[CELLTYPE][INSTANCE] = findMax(max);
  return CELL_Key;
}
void SDFParser::Parse()
{
  string s;
  char c;
//  cout << s << endl;
  bool cell_key = 0;
  while(!fp.eof()){
    std::string line;
    if(cell_key == 0)
      std::getline(fp, line);
    else{
      line = "CELL";
      cell_key = 0;
    }
    if(line.find("CELL") != -1)
    {
        line.clear();
        std::getline(fp, line);
        int pos = line.find("CELLTYPE");
        if(pos != -1){
          int ind1 = line.find("\"");
          std::string CELLType = line.substr(ind1 + 1);
          CELLType = CELLType.substr(0, CELLType.size() - 3);
          line.clear();
          std::getline(fp, line);
          ind1 = line.find("INSTANCE");
          std::string INSTANCE = line.substr(ind1 + 8);
          INSTANCE = INSTANCE.substr(0, INSTANCE.size() - 2);
          if (CELLType.find("DFF") != -1) {
            getDFF(INSTANCE, CELLType);
          }
          else{
            cell_key  = getCOMB(INSTANCE, CELLType);
          }
        }
    }
  }

}
int main() {
  SDFParser ps("simplified_sdf.sdf");
  ps.Parse();
  DelayData dData;
  dData = ps.getData();

  string DFF_Model = "DFF_X1";
  string bank_selection = " bank_selection";
  string dec = " dec";
  string ir = " ir";
  string pc = " pc";
  string MUX_Model = "AOI22_X1";
  string  mux2_1 = " mux2_1";
  string  mux2_2 = " mux2_2";
  string  mux2_3 = " mux2_3";
  string  mux2_4 = " mux2_4";
  string  mux2_5 = " mux2_5";
  string  FA_Model = "FA_X1";
  string  adder_4 = " adder_4";
  string HA_Model = "HA_X1";
  string  add_75 = " add_75";
  string NOR_Model = "NOR4_X1";
  string  ins_mem = " ins_mem";
  Graph Max;
  float pc_to_ir = dData.DFFdelay[DFF_Model][pc].maxDelay + dData.maxDelay[NOR_Model][ins_mem] + dData.DFFdelay[DFF_Model][ir].Setup ;
  float pc_to_pc = dData.DFFdelay[DFF_Model][pc].maxDelay + dData.DFFdelay[DFF_Model][pc].Setup + dData.maxDelay[FA_Model][adder_4] + dData.maxDelay[MUX_Model][mux2_1] + dData.maxDelay[MUX_Model][mux2_5] ;

  Max.add_node(pc, {{ir, pc_to_ir}, {pc, pc_to_pc}});

  float ir_to_dec = dData.DFFdelay[DFF_Model][ir].maxDelay + dData.DFFdelay[DFF_Model][dec].Setup + std::max(std::max( dData.maxDelay[MUX_Model][mux2_3],  dData.maxDelay[MUX_Model][mux2_4]), dData.maxDelay[MUX_Model][mux2_2]) ;
  float ir_to_pc =  dData.maxDelay[MUX_Model][mux2_5] + dData.DFFdelay[DFF_Model][ir].maxDelay + dData.DFFdelay[DFF_Model][pc].Setup + std::max( dData.maxDelay[FA_Model][adder_4] + dData.maxDelay[MUX_Model][mux2_1] + dData.maxDelay[MUX_Model][mux2_5], dData.maxDelay[MUX_Model][mux2_1] + dData.maxDelay[MUX_Model][mux2_5]) ;
  Max.add_node(ir, {{dec, ir_to_dec}, {pc, ir_to_pc}});
  Max.add_node(dec, {});
  float b_to_dec = dData.DFFdelay[DFF_Model][bank_selection].maxDelay + dData.DFFdelay[DFF_Model][dec].Setup + std::max(dData.maxDelay[MUX_Model][mux2_3], dData.maxDelay[MUX_Model][mux2_4]) ;
  Max.add_node(bank_selection, {{dec, b_to_dec}});
  cout << "MAXIMUM DELAY ANALYSIS RESULT \n";
  cout << "---------------------------------------------------\n";
  cout << "PC To PC Maximum Delay                                     : " << Max.LongestPath(pc,pc) << endl;
  cout << "PC To IR Maximum Delay                                     : " << Max.LongestPath(pc,ir) << endl;
  cout << "IR To Decode Maximum Delay                                 : " << Max.LongestPath(ir, dec) << endl;
  cout << "IR To PC Maximum Delay                                     : " << Max.LongestPath(ir, pc) << endl;
  cout << "Bank Selection To Decode Maximum Delay                     : " << Max.LongestPath(bank_selection, dec) << endl;


  Graph Min;
  float pc_to_ir_min = dData.DFFdelay[DFF_Model][pc].minDelay + dData.minDelay[NOR_Model][ins_mem] - dData.DFFdelay[DFF_Model][ir].Hold   ;
  float pc_to_pc_min = dData.DFFdelay[DFF_Model][pc].minDelay - dData.DFFdelay[DFF_Model][pc].Hold + dData.minDelay[HA_Model][add_75] + dData.minDelay[MUX_Model][mux2_1]  ;

  Min.add_node(pc, {{ir, pc_to_ir_min}, {pc, pc_to_pc_min}});

  float ir_to_dec_min = dData.DFFdelay[DFF_Model][ir].minDelay - dData.DFFdelay[DFF_Model][dec].Hold + std::min(std::min( dData.minDelay[MUX_Model][mux2_3],  dData.minDelay[MUX_Model][mux2_4]), dData.minDelay[MUX_Model][mux2_2])   ;
  float ir_to_pc_min =  dData.minDelay[MUX_Model][mux2_5] + dData.DFFdelay[DFF_Model][ir].minDelay - dData.DFFdelay[DFF_Model][pc].Hold + std::min( dData.minDelay[FA_Model][adder_4] + dData.minDelay[MUX_Model][mux2_1] + dData.minDelay[MUX_Model][mux2_5], dData.minDelay[MUX_Model][mux2_1] + dData.minDelay[MUX_Model][mux2_5])   ;
  Min.add_node(ir, {{dec, ir_to_dec_min}, {pc, ir_to_pc_min}});
  Min.add_node(dec, {});
  float b_to_dec_min = dData.DFFdelay[DFF_Model][bank_selection].minDelay - dData.DFFdelay[DFF_Model][dec].Hold + std::min(dData.minDelay[MUX_Model][mux2_3], dData.minDelay[MUX_Model][mux2_4])  ;
  Min.add_node(bank_selection, {{dec, b_to_dec_min}});
  map<string, float> PC_ShortPaths = Min.Diakstra(pc);
  map<string, float> IR_ShortPaths = Min.Diakstra(ir);
  map<string, float> BS_ShortPaths = Min.Diakstra(bank_selection);
  // Diakstra can not calculate path for Crown Graph
  if(PC_ShortPaths[pc] == 0)
      PC_ShortPaths[pc] = std::min(pc_to_pc_min, pc_to_ir_min + ir_to_pc_min);
  cout << "\nMINIMUM DELAY ANALYSIS RESULT \n";
  cout << "---------------------------------------------------\n";
  cout << "PC To PC Minimum Delay                  -    hold time     : " << PC_ShortPaths[pc] << endl;
  cout << "PC To IR Minimum Delay                  -    hold time     : " << PC_ShortPaths[ir] << endl;
  cout << "IR To Decode Minimum Delay              -    hold time     : " << IR_ShortPaths[dec] << endl;
  cout << "IR To PC Minimum Delay                  -    hold time     : " << IR_ShortPaths[pc] << endl;
  cout << "Bank Selection To Decode Minimum Delay  -    hold time     : " << BS_ShortPaths[dec] << endl;

  return 0;
}
