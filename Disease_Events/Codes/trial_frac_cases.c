#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <time.h> //Used for clock

#define year 365 //days
#define T_MAX 3*year //In days
#define deltat 1 //In days
#define AG 15
#define countries 3
#define N_trials 10000
#define Ncontrol 1639
#define Nvaccine 1626


void start_sto_sim(double *AG_C_SA, double *AG_V_SA, double *AG_C_K, double *AG_V_K, double *AG_C_Z, double *AG_V_Z);
int eventDrivenTrial(int a, int i_ctr, int iter, int NL, int NF, double beta, double p, double q, double r, double r_L);
double bern_binom(double p, double n);
double gaussrand();

FILE *f;

void start_sto_sim(double *AG_C_SA, double *AG_V_SA, double *AG_C_K, double *AG_V_K, double *AG_C_Z, double *AG_V_Z){
    int age;
    int scan, i, j, iter;
    int D_c;
    FILE *fp;
    
    double PRC[3][2][AG];
    for(i=0; i<AG; i++){
        PRC[0][0][i] = AG_C_SA[i];
        PRC[0][1][i] = AG_V_SA[i];
        PRC[1][0][i] = AG_C_K[i];
        PRC[1][1][i] = AG_V_K[i];
        PRC[2][0][i] = AG_C_Z[i];
        PRC[2][1][i] = AG_V_Z[i];
    }
    
    int total_population_control = Ncontrol;
    
    double L_frac[3][501][AG]; //[country][Model iteration][age group]; countries are loaded as 0 -> South Africa, 1 -> Kenya, 2 -> Zambia
    
    //Loading F_L distribution in South Africa;
    fp=fopen("Data_countries/Fracs_INF_South_Africa.txt", "rt");
    for(iter=0; iter<501; iter++){
        for(i=0; i<AG; i++){
            scan = fscanf(fp, "%lf ",&L_frac[0][iter][i]);
        }
    }
    fclose(fp);
    
    //Loading F_L distribution in Kenya;
    fp=fopen("Data_countries/Fracs_INF_Kenya.txt", "rt");
    for(iter=0; iter<501; iter++){
        for(i=0; i<AG; i++){
            scan = fscanf(fp, "%lf ",&L_frac[1][iter][i]);
        }
    }
    fclose(fp);
    
    //Loading F_L distribution in Zambia;
    fp=fopen("Data_countries/Fracs_INF_Zambia.txt", "rt");
    for(iter=0; iter<501; iter++){
        for(i=0; i<AG; i++){
            scan = fscanf(fp, "%lf ",&L_frac[2][iter][i]);
        }
    }
    fclose(fp);
    
    double param[3][501][6];
    fp=fopen("Data_countries/param_South_Africa.txt", "rt");
    for(iter=0; iter<501; iter++){
        scan = fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &param[0][iter][0], &param[0][iter][1], &param[0][iter][2], &param[0][iter][3], &param[0][iter][4], &param[0][iter][5]);
    }
    fclose(fp);
    
    fp=fopen("Data_countries/param_Kenya.txt", "rt");
    for(iter=0; iter<501; iter++){
        scan = fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &param[1][iter][0], &param[1][iter][1], &param[1][iter][2], &param[1][iter][3], &param[1][iter][4], &param[1][iter][5]);
    }
    fclose(fp);
    
    fp=fopen("Data_countries/param_Zambia.txt", "rt");
    for(iter=0; iter<501; iter++){
        scan = fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &param[2][iter][0], &param[2][iter][1], &param[2][iter][2], &param[2][iter][3], &param[2][iter][4], &param[2][iter][5]);
    }
    fclose(fp);
    
    double beta_inf[3][501][AG];
    //Loading force of infection South Africa;
    fp=fopen("Data_countries/F_INF_South_Africa.txt", "rt");
    for(iter=0; iter<501; iter++){
        for(i=0; i<AG; i++){
            scan = fscanf(fp, "%lf ",&beta_inf[0][iter][i]);
        }
    }
    fclose(fp);
    //Loading force of infection Kenya;
    fp=fopen("Data_countries/F_INF_Kenya.txt", "rt");
    for(iter=0; iter<501; iter++){
        for(i=0; i<AG; i++){
            scan = fscanf(fp, "%lf ",&beta_inf[1][iter][i]);
        }
    }
    fclose(fp);
    //Loading force of infection Zambia;
    fp=fopen("Data_countries/F_INF_Zambia.txt", "rt");
    for(iter=0; iter<501; iter++){
        for(i=0; i<AG; i++){
            scan = fscanf(fp, "%lf ",&beta_inf[2][iter][i]);
        }
    }
    fclose(fp);
    
    //--------------------------------------------------------------------------//
    
    //Initiate C's own random generator.
    srand(time(NULL));
    
    int iter_select;
    double p, q, r, r_L;
    
    int L_ini = 0, F_ini = 0, calc_aux = 0;
    
    int i_c = 0;
    double D_c_age[3];
    double F_con[3];
    
    for(j = 0; j < N_trials; j++){
        for(i_c=0; i_c<3; i_c++){
            D_c_age[i_c] = 0;
            F_con[i_c] = 0;

                
            iter_select = (double)rand()/RAND_MAX*501;
            
            q = param[i_c][iter_select][3];
            r = param[i_c][iter_select][4];
            r_L = param[i_c][iter_select][5];
            
            int aux_starting_pop;
            
            for(age = 0; age < 15; age++){
                
                //Control
                aux_starting_pop = total_population_control; //ages_cont[age];
                calc_aux = round(L_frac[i_c][iter_select][age]/100 * aux_starting_pop); 
                
                L_ini = calc_aux;
                F_ini = aux_starting_pop - calc_aux;
                
                if(age==0){
                    p = param[i_c][iter_select][0];
                }
                else if(age==1){
                    p = param[i_c][iter_select][1];
                }
                else{
                    p = param[i_c][iter_select][2];
                }
                
                D_c = eventDrivenTrial(age, i_c, j, L_ini, F_ini, beta_inf[i_c][iter_select][age], p, q, r, r_L);

                D_c_age[i_c] += (double)D_c*PRC[i_c][0][age];
                
                F_con[i_c] += (double)F_ini*PRC[i_c][0][age];
            }
        }
    }
}

double bern_binom(double p, double n){
    int i, c=0;
    double aux;
    for(i=0; i<n; i++){
        aux = (double)rand()/RAND_MAX;
        if(aux < p){
            c++;
        }
    }
    return c;
}

double gaussrand(){
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}

int eventDrivenTrial(int a, int i_ctr, int iter, int NL, int NF, double beta, double p, double q, double r, double r_L){
    double t=0, dt, P;
    int i, n_l=0, n_f=0, n_i=0, n_i_total = 0, calc_aux = 0;
    int n_f_reinf=0;
    
    double re[4], R=0, rn;

    //printf("Control cohort");
    n_l=NL;
    n_f=NF;
    
    FILE *out;
    char filename[500];
    //sprintf(filename, "Case_Dist/Fracs_%02d_%d.txt", a, i_ctr);
    sprintf(filename, "Case_Dist/Fracs_%d.txt", i_ctr);
    //if(iter==0){
    if(iter==0 && a==0){
        out = fopen(filename, "w");
        fprintf(out, "Age Iter L_ini F_ini Dfin Reinf ReinfFP  FastProg SlowProg\n");
    }
    else{
        out = fopen(filename, "a");
    }
    
    int reinf=0, slow_prog=0, fast_prog=0, reinf_fp=0;
    
    t=0;
    while(t < T_MAX/(double)year){
          
        re[0] = r_L*n_l; //Reactivation from Latency
        re[1] = beta*q*p*n_l; //Reinfection in Latency
        re[2] = r*n_f_reinf; //Fast progress 2
        re[3] = r*n_f; //Fast progress
        

        R=0;
        for(i=0; i<4; i++)
            R += re[i];
    
        rn = (double)rand()/RAND_MAX; //generador();
        while(rn<1E-24){
            rn = (double)rand()/RAND_MAX;// (double)rand()/RAND_MAX;//generador();
        }
        dt = log(rn)*(-1.0/R);
        P = (double)rand()/RAND_MAX*R; //generador()*R;

        if(P < re[0]){ //Slow progress
            if(n_l>0){
                n_l -= 1;
                n_i += 1;
            }
            slow_prog++;
        }
        else if( P < (re[0] + re[1]) ){ //Reinfection
            if(n_l>0){
                n_l -= 1;
                n_f_reinf += 1;
            }
            reinf++;
        }
        else if( P < (re[0] + re[1] + re[2]) ){ //Reinfection + fast_progress
            if(n_f_reinf>0){ //Fast progress
                n_f_reinf -= 1;
                n_i += 1;
            }
            reinf_fp++;
        }
        else{
            if(n_f>0){ //Fast progress
                n_f -= 1;
                n_i += 1;
            }
            fast_prog++;
        }
        
        t+=dt;
    }
    
    fprintf(out, "%d %d %d %d %d %d %d %d %d\n", a, iter, NL, NF, n_i, reinf, reinf_fp, fast_prog, slow_prog);
    fclose(out);
    n_i_total = n_i;
    
    return n_i_total;
}









