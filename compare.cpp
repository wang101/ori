#include<stdio.h>

double compare(const char* groundtruth_file,const char* test_file,int scr_size)
{
    FILE* grd;
    FILE* comp;
    int nrow,ncol;
    double t1,t2;
    double sum=0;
    grd = fopen(groundtruth_file,"r");
    comp = fopen(test_file,"r");
    if(grd==NULL || comp==NULL){
        printf("can not open file to compare!\n");
        return -2;
    }
    nrow = scr_size;
    ncol = scr_size;
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            fscanf(grd,"%lf,",&t1);
            fscanf(comp,"%lf,",&t2);
            sum += (t1-t2)*(t1-t2)/(nrow*ncol);
            if(feof(grd) ^ feof(comp))
            {
                printf("image size doesn't match, please check your ground-truth file\n");
                fclose(grd);
                fclose(comp);
                return -1;
            }
        }
    }
    fclose(grd);
    fclose(comp);
    return sum;
}

