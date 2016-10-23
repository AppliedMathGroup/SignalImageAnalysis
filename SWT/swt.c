//
//  swt.c
//  
//
//  Created by Ricardo Monge on 10/23/16.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <clapack.h>
#endif
/* global variables */
int f; /* frequency */
double Dt;/* Delta Time */
int n; /* size of input data */
double* V; /* holds input data */
double* V1; /* holds copy of input data */
double* C; /* holds coefficient matrix */
double* SF; /* holds frequencies for each data point */
double* CV; /* holds approximation for verification*/
double Dm;

/* function prototypes */
void show_help();
void load_file(char*);
void create_coefficient_matrix();
void solve_matrix();
void compute_frequencies();
void compute_approximation();
void compute_dm();
void output_data();

/* main program code */
int main(int argc, char** argv) {
    /* Verify input arguments */
    
    if(argc!=4)   show_help();
    
    load_file(argv[3]);
    f=atoi(argv[1]);
    Dt=atof(argv[2]);
    printf("Frequency: %d\n",f);
    printf("Delta Time: %f\n",Dt);
    printf("Number of samples: %d\n",n);
    /* create matrices */
    C = (double*)malloc(sizeof(double)*n*n);
    SF = (double*)malloc(sizeof(double)*n);
    CV = (double*)malloc(sizeof(double)*n);
    
    
    create_coefficient_matrix();
    solve_matrix();
    compute_frequencies();
    compute_approximation();
    compute_dm();
    output_data();
}

/* implementation of SWT specific commands */
void create_coefficient_matrix(){
    int i,j,lj, pos;
    int q,r,qmod2;
    for(i=1;i<=n;i++)
        for(j=1;j<=n;j++) {
            lj=n-j+1;
            pos=(i)+(j-1)*n;
            q=i/lj;
            r=i%lj;
            qmod2=q%2;
            if(qmod2==0 && r==0) C[pos-1]=-1;
            if(qmod2==0 && r!=0) C[pos-1]=1;
            if(qmod2!=0 && r==0) C[pos-1]=1;
            if(qmod2!=0 && r!=0) C[pos-1]=-1;
        }
    
}
void solve_matrix(){
    //double A[9] = {76, 27, 18, 25, 89, 60, 11, 51, 32};
    //double b[3] = {10, 7, 43};
    
    int N = n;
    int nrhs = 1;
    int lda = n;
    int ipiv[n];
    int ldb = n;
    int info;
    
    dgesv_(&N, &nrhs, C, &lda, ipiv, V, &ldb, &info);
    if(info!=0) {
        fprintf(stderr, "dgesv_ fails %d\n", info);
        exit(5);
    }
    /* solution is in V */
}
void compute_frequencies(){
    int i;
    double t=1/(2*Dt);
    for(i=1;i<=n;i++) {
        SF[i-1]=t*(n/(n-i+1.0));
    }
}
void compute_approximation(){
    double SQW[n];
    int i,j;
    double ti;
    for(i=1;i<=n;i++) SQW[i-1]=2*Dt*(float)i/n;
    double de=Dt/n;
    for(i=1;i<=n;i++) {
        ti=i*de+de/2;
        CV[i-1]=0;
        for(j=1;j<n;j++)
            CV[i-1]+=V[i-1]*pow(-1,floor(2*ti/SQW[i-1]));
    }
}
void compute_dm(){
    int i,j;
    Dm=0;
    double k;
    for(i=1;i<n;i++) {
        k=fabs(V1[i-1]-CV[i-1]);
        if(k>Dm) Dm=k;
    }
}

void output_data(){
    int i;
    for(i=1;i<n;i++) {
        printf("%d: (%f,%f) \n",i,SF[i-1],V[i-1]);
    }
    //printf("\nDm=%f\n\n",Dm);
}
/* implementation of helper functions */
void show_help(){
    printf("SWT frequency time datafile\n\n");
    exit(1);
}

void load_file(char* fname) {
    int lines_allocated = 128;
    int max_line_len = 100;
    /* read data into buffer */
    char **words = (char **)malloc(sizeof(char*)*lines_allocated);
    if (words==NULL)
    {
        fprintf(stderr,"Out of memory (1).\n");
        exit(1);
    }
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        fprintf(stderr,"Error opening file.\n");
        exit(2);
    }
    int i;
    for (i=0;1;i++)
    {
        int j;
        
        /* Have we gone over our line allocation? */
        if (i >= lines_allocated)
        {
            int new_size;
            
            /* Double our allocation and re-allocate */
            new_size = lines_allocated*2;
            words = (char **)realloc(words,sizeof(char*)*new_size);
            if (words==NULL)
            {
                fprintf(stderr,"Out of memory.\n");
                exit(3);
            }
            lines_allocated = new_size;
    }
        /* Allocate space for the next line */
        words[i] = malloc(max_line_len);
        if (words[i]==NULL)
        {
            fprintf(stderr,"Out of memory (3).\n");
            exit(4);
        }
        if (fgets(words[i],max_line_len-1,fp)==NULL)
            break;
        
        /* Get rid of CR or LF at end of line */
        for (j=strlen(words[i])-1;j>=0 && (words[i][j]=='\n' || words[i][j]=='\r');j--)
            ;
        words[i][j+1]='\0';
    }
    /* Close file */
    fclose(fp);
    
    
    /* convert buffer into array of floats*/
    V = (double*)malloc(sizeof(double)*i);
    V1 = (double*)malloc(sizeof(double)*i);
    
    int j;
    for(j = 0; j < i; j++)
    {
        V[j]=atof(words[j]);
        V1[j]=V[j];
    }
    n=i;
    
    /* Good practice to free memory */
    for (;i>=0;i--)
        free(words[i]);
    free(words);
    
    
}










