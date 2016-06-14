#include "model.cpp"
#include "compare.cpp"
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using std::vector;
using std::string;

int main(int argc, char* argv[])
{
    char eulerdef[3];
    double euler[3];
    double diff;
    para p;//parameters
    string groundtruth_file;//patterns from real experiment
    string protein_file;//read protein file path
    string filepatt;
    string savepath;
    ifstream ipt;
    vector<string> savefile;
    stringstream temp;
    clock_t start,end;

    ipt.open("task");
	ipt >> groundtruth_file;
	ipt >> protein_file;
    ipt >> savepath;
    ipt >> eulerdef[0] >>  eulerdef[1] >>  eulerdef[2]; 

    start = clock();
    int k=0;
    
    printf("initial model...\n");
    model mdl = model(p);
    printf("reading protein...\n");
    mdl.read_protein(protein_file.c_str());
    printf("reading groundtruth...\n");
    int flag = mdl.read_groundtruth(groundtruth_file.c_str());
    while(!ipt.eof()){
        ipt >> euler[0] >>  euler[1] >> euler[2];
        if(ipt.eof())break;
        printf("setting euler angle as:%4.2lf,%4.2lf,%4.2lf;%c,%c,%c respectively\n",\
                euler[0],euler[1],euler[2],eulerdef[0],eulerdef[1],eulerdef[2]);
        mdl.set_euler(euler,eulerdef);
        printf("generating pattern\n");
        mdl.gene_pattern();
        temp.clear();
        temp << savepath <<  k++ << ':' << euler[0] << ',' <<  euler[1] << ',' << euler[2] << ".dat";
        temp >> filepatt ;
        std::cout << filepatt << endl;
        savefile.push_back(filepatt);
        printf("save pattern,filename:%s\n",filepatt.c_str());
        mdl.save_pattern(filepatt.c_str());
    }
    end = clock();
    printf("USE TIME:%6.2fsec\n",(double)(end-start)/CLOCKS_PER_SEC);
}
