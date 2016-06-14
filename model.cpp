#include<math.h> 
#include<vector>
#include<stdio.h>
#include<complex>
#include<stdlib.h>
#include<string.h>
#include"parameter.cpp"
#include "vec3.cpp"
#define PI 3.1415926535897
#define eI (complex<double>(0,1.0))

using namespace std;

const double Angstrom = 1.0e-10;

class model
{
	public:
		double pix_size;//perimeter of each pixel
		int scr_size;//screen size
        double wave_number;//wavenumber = 1/lambda
        double distance;//distance from crystal to screen
        double epsilon;//parameter reserved for calculation of crystal
        double euler[3];
        char eulerdef[3];
        int natom;
        double RT[3][3];
        model(para p)
		{
            wave_number = 1/p.lambda;
			scr_size = p.scr_size;
            ncol = scr_size;
            nrow = scr_size;
            pix_size = p.pix_size;
            distance = p.distance;
            epsilon = p.epsilon;
            natom = p.natom;
            atom_crd = (double**)calloc(natom,sizeof(double *));
            atom_crd[0] = (double *)calloc(natom*3,sizeof(double));
            for(int i=1;i<natom;i++)
                atom_crd[i] = atom_crd[0] + i*3;
            scr_crd = (vec**)calloc(scr_size,sizeof(vec *));
            scr_crd[0] = (vec*)calloc(scr_size*scr_size,sizeof(vec));
            for(int i=0;i<scr_size;i++)
                scr_crd[i] = scr_crd[0] + scr_size*i;
            trans_crd = (vec*)calloc(natom,sizeof(vec));
            patt = (complex<double>*) calloc(scr_size*scr_size,sizeof(complex<double>));
            groundtruth = (double*) calloc(scr_size*scr_size,sizeof(double));
            printf("checking parameters:\n");
            printf("screensize:%d pixels, %4.2lf meters\n",scr_size,scr_size*pix_size);
            printf("lambda:%4.2e\n distance to screen:%4.2f\n",p.lambda,distance);
		}
        void set_euler(double a[3],char def[3]);
        int read_groundtruth(const char* groundtruth_file);
        void read_protein(const char* protein_file);
        void gene_pattern();
        double comparenow();
        int save_pattern(const char* savefile);
	private:
        complex<double>* patt;
        double* groundtruth;
        void Rot(double t,char a,double Rn[3][3]);
		void rotate();
        int ncol,nrow;
        double** atom_crd;
        vec* trans_crd;
        vec** scr_crd;
};

void model::Rot(double t,char a,double Rn[3][3])
{
    double rx[3][3] = {{1,0,0},{0,cos(t),-sin(t)},{0,sin(t),cos(t)}};
    double ry[3][3] = {{cos(t),0,sin(t)},{0,1,0},{-sin(t),0,cos(t)}};
    double rz[3][3] = {{cos(t),-sin(t),0},{sin(t),cos(t),0},{0,0,1}};
    if(a=='x'){
        memcpy(Rn,rx,sizeof(double)*9);
    }else if(a=='y'){
        memcpy(Rn,ry,sizeof(double)*9);
    }else if(a=='z'){
        memcpy(Rn,rz,sizeof(double)*9);
    }
}

void model::set_euler(double a[3],char def[3])
{
    euler[0] = a[0];
    euler[1] = a[1];
    euler[2] = a[2];
    eulerdef[0] = def[0];
    eulerdef[1] = def[1];
    eulerdef[2] = def[2];
}



void model::rotate()
{
	double R1[3][3];
    double R2[3][3];
    double R3[3][3];
    double Rall[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            RT[i][j]=0;

    for(int i=0;i<3;i++)
    {
        euler[i] = euler[i]*PI/180.0;
    }
    Rot(euler[0],eulerdef[0],R1);
    Rot(euler[1],eulerdef[1],R2);
    Rot(euler[2],eulerdef[2],R3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                Rall[i][j]+=R1[i][k]*R2[k][j];
            }
        }
    }
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                RT[i][j]+=Rall[i][k]*R3[k][j];
            }
        }
    }
    for(int k=0;k<natom;k++)
    {
        double rr1,rr2,rr3;
        rr1 = atom_crd[k][0]*Angstrom;
        rr2 = atom_crd[k][1]*Angstrom;
        rr3 = atom_crd[k][2]*Angstrom;
        double nr1 = (RT[0][0]*rr1 + RT[0][1]*rr2 + RT[0][2]*rr3);
        double nr2 = (RT[1][0]*rr1 + RT[1][1]*rr2 + RT[1][2]*rr3);
        double nr3 = (RT[2][0]*rr1 + RT[2][1]*rr2 + RT[2][2]*rr3);
        trans_crd[k] = vec(nr1,nr2,nr3);
    }
}
 
void model::read_protein(const char* proteinfile)
{
    FILE* file = fopen(proteinfile,"r");
    fscanf(file,"%d",&natom);
    for(int i=0;i<natom;i++)
        fscanf(file,"%lf,%lf,%lf",&atom_crd[i][0],&atom_crd[i][1],&atom_crd[i][2]);
    fclose(file);
    //init scr
    double ssc = scr_size/2.0-0.5;
    for(int i=0;i<scr_size;i++){
        for(int j=0;j<scr_size;j++){
            double x = (i-ssc)*pix_size;
            double y = (j-ssc)*pix_size;
            double z = distance;
            double ex,ey,ez;
            double r = sqrt(x*x+y*y+z*z);
            ex = wave_number*x/r;
            ey = wave_number*y/r;
            ez = wave_number*z/r;
            scr_crd[i][j] = vec(ex,ey,ez-wave_number);
        }
    }
}

int model::read_groundtruth(const char* groundtruth_file)
{
    FILE* file = fopen(groundtruth_file,"r");
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            if(feof(file)){
                fclose(file);
                printf("groundtruth file invalid,compare aborted\n");
                return 0;
            }
            fscanf(file,"%lf,",(groundtruth+i*ncol+j));
        }
    }
    fclose(file);
    return 1;
}


void model::gene_pattern()
{
    rotate();

    for(int j=0;j<scr_size;j++){
        for(int i=0;i<scr_size;i++){
            *(patt+j*scr_size+i) = 0;
            for(int k=0;k<natom;k++)
            {
                complex<double> sglAtm = exp(2.0*PI*scr_crd[j][i].dot(trans_crd[k])*eI);
                *(patt+j*scr_size+i) += sglAtm ;
            }
        }
    }
}


int model::save_pattern(const char* savefile)
{
    FILE* prpat=NULL;
    prpat = fopen(savefile,"w");
	for(int j=0;j<scr_size;j++){
		for(int k=0;k<scr_size;k++){
		    fprintf(prpat,"%f,", abs(*(patt+scr_size*j+k)));
		}
		fprintf(prpat,"\n");
    }
    fclose(prpat);
    return 0;
}

double model::comparenow()
{
    double sum=0;
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            double t1 = *(groundtruth+i*ncol+j);
            double t2 = abs(*(patt+i*ncol+j));
            sum += (t1-t2)*(t1-t2)/(ncol*nrow);
        }
    }
    return sum;
}
/*    
double compare(const char* gt,const char* cmp,int scr_size)
{
    FILE* grd;
    FILE* comp;
    int nrow,ncol;
    double t1,t2;
    double sum=0;
    grd = fopen(gt,"r");
    comp = fopen(cmp,"r");
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
*/
