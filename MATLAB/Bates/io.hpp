//
//  io.hpp
//  Jul2022
//
//  Created by Antonio Santos on 19/07/2022.
//  Copyright © 2022 Antonio Santos. All rights reserved.
//

#ifndef io_hpp
#define io_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;
//

class IO{
public:
    IO(){}
    
    vector<double> readVector(const char* file,int *n)
    {
        FILE *fp;
        fp = fopen(file, "r");
        char buff[64];
        int nrow = -1;
        while (!feof(fp)) {
            fscanf(fp, "%s",buff);
            nrow++;
        }
        *n = nrow;
        rewind(fp);
        vector<double> result = vector<double>(nrow);
        for(int i=0;i<nrow;i++){
            fscanf(fp, "%s",buff);
            result[i] = atof(buff);
        }
        fclose(fp);
        return result;
    }
    //
    vector<int> readVectorInt(const char* file,int *n)
    {
        FILE *fp;
        fp = fopen(file, "r");
        char buff[64];
        int nrow = -1;
        while (!feof(fp)) {
            fscanf(fp, "%s",buff);
            nrow++;
        }
        *n = nrow;
        rewind(fp);
        vector<int> result = vector<int>(nrow);
        for(int i=0;i<nrow;i++){
            fscanf(fp, "%s",buff);
            result[i] = atoi(buff);
        }
        fclose(fp);
        return result;
    }

};



#endif /* io_hpp */
