#ifndef io_cuh
#define io_cuh
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
//
//
template <typename T>
class IO{
	public:
		IO(){}
		//
		//
		T* readVector(const char *file,int *n)
		{
			FILE *fp;
			fp = fopen(file,"r");
			char buff[64];
			int nrow = -1;
			while(!feof(fp)){
				int tt = fscanf(fp,"%s",buff);
				nrow++;
			}
			*n = nrow;
			rewind(fp);
			T* result = new T[nrow];
			for(int i=0;i<nrow;i++){
				int tt = fscanf(fp,"%s",buff);
				result[i] = (T)atof(buff);
			}
			fclose(fp);
			return result;
		}
		//
		//
		T* readVectorCUDA(const char *file,int *n)
		{
			FILE *fp;
			fp = fopen(file,"r");
			char buff[64];
			int nrow = -1;
			while(!feof(fp)){
				int tt = fscanf(fp,"%s",buff);
				nrow++;
			}
			*n = nrow;
			rewind(fp);
			T* result;
			cudaMallocManaged(&result,nrow*sizeof(T));
			for(int i=0;i<nrow;i++){
				int tt = fscanf(fp,"%s",buff);
				result[i] = (T)atof(buff);
			}
			fclose(fp);
			return result;
		}

		//
		void writeVector(T *x,const char *file,int n,int ndigits=4)
		{
			FILE *fp;
			char str[6];
			sprintf(str,"%s.%df\\n","%",ndigits);
			fp = fopen(file,"wa");
			for(int i=0;i<n;i++) fprintf(fp,str,x[i]);
			fclose(fp);
		}

};

#endif
