#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "string.h"
#include "time.h"
#include "unistd.h"

//Common functions for runner and fitter
#include "TB_functions.h"

//#define INF_MATRIX_CALC

//Models:
int model_run(char filedem[string_length],int start,int end, double *objeto);

//Loading functions
void load_dynamic_variables(char demo_r[string_length], double s[states][AG][EG],double S[states],int start,int end, int is_file);
void load_initial_conditions(double s[states][AG][EG],double S[states]);
void load_fit_pyramids(double s[states][AG][EG],char pyr[string_length],int inicio,int end);
void load_initial_conditions_from_s_resource(double s[states][AG][EG],double S[states]);
void read_pyr_coef();

//Recording functions
void save_output_micro(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int year, double *vv);
void save_output_micro_file(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int start);
void rewrite_output_object();
void write_output_run_inf_matrix(int uncertainty, int phase,int iter,double inf[AG][AG],  double inf_1[AG][AG], double inf_2[AG][AG], double inf_3[AG][AG]);

//Integrator
int rk4(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG], double S[states], double f_micro[AG][EG][internal_fluxes], double t);
int deriv(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG], double ds[states][AG][EG], double f_micro[AG][EG][internal_fluxes], double t,double S[states]);
int calculate_fluxes(int a, int e,double t,double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double f[fluxes],double S[states]);

double cdf(int a, int av, double tw, double sigma);
void ageing(double ds[states][AG][EG], double s[states][AG][EG],double S[states],double vp[variable_parameters][AG][EG],double t);
double semiempiric_derivative_fitted_population(int a,int e,double t);

char gp_global[string_length], rp_global[string_length], vp_global[string_length], cm_global[string_length], demo_global[string_length];
int phase, phase_temp;
double b[EG][max_records*4*3*H1], mu[AG][EG][max_records*4*3*H1], pyramid_coefficient_aux[AG][EG][fitting_degree+1][2];
double s_stat[states][AG][EG],S_stat[states];


double scan;

int temp_iter=0;

void R2C_runner(char **mf, int *adjust, int *vaciter, int *ittemp, int *con_mod, int *vaccinate, double *objeto1, double *objeto2, int *out_handler){
    /* Running has 4 phases */
	/* 1º Baseline until the year of the vaccine */
	/* 2º Baseline from year of the vaccine till the end */
	/* 3ª Vaccination until the year of the vaccine (same as baseline) */
	/* 4º Vaccination from year of the vaccine till the end */
	/* We always need the baseline to obtain the evolution of the demographic that we will apply on the run of the vaccine */

    temp_iter=*ittemp;
    
    double correction, sick = 0;
    int a, e, a1, i, span, err=0;
    char aux_str[string_length], aux_command[string_length];

    strcpy(master_file, mf[0]);
    master_file_reader(master_file);
    load_window(window_file);
    span = year_max - year_min + 1;

    cm_mode_all = *con_mod;

	printf("____________________________________________________________________\n");
	printf("RUNNING: Run dynamics between %d and %d with vaccine in %d\n",year_min,year_end,year_vac);

	strcpy(gp_global, global_parameters_file);
	strcpy(rp_global, regional_parameters_file);
	strcpy(cm_global, contacts_file);
	strcpy(demo_global, demography_file);

	iter=0;
	uncertainty=0;
	phase=0;
	phase_temp=0;

	load_parameters(gp_global, global_parameters, gp);
	load_parameters(rp_global, regional_parameters, rp);
    load_contact_matrix(cm_global, cm);
	read_pyr_coef();

    if (strcmp(name_country, simulation_name) != 0)
        sprintf(aux_str, "../Outputs/%s/Iter/%s/Fitted_inputs/Parameters_to_fit.txt", name_country, simulation_name);
    else
	    sprintf(aux_str, "../Outputs/%s/Fitted_inputs/Parameters_to_fit.txt", simulation_name);

    double aux_h[5][AG];
    load_fitting_parameters(aux_str, gp, rp, aux_h);

    load_dynamic_variables(demo_global, s, S, year_min, year_end, 1);

    for(a=0;a<AG;a++){
		for(e=0;e<EG;e++){
			for(i=0;i<regional_parameters;i++)
				rp_base[i][a][e]=rp[i][a][e];
		}
	}

    for(a=0; a<AG; a++){
        for(e=0; e<EG; e++){
	        rp[3][a][e] = rp_base[3][a][e]*aux_h[0][a]; //Initial diagnosis rate smear-pos
			for(i=4; i<6; i++)
				rp[i][a][e] = rp_base[i][a][e]*aux_h[0][a]*rp[18][a][e];//Diagnosis rate smear-neg. and extra-pulm, times eta
			rp[6][a][e] = aux_h[1][a];
            rp[7][a][e] = aux_h[2][a];
			rp[8][a][e] = aux_h[3][a];
            rp[2][a][e] = aux_h[4][a];
        }
	}

	for(a=0; a<AG; a++)
	{
		for(e=0; e<EG; e++)
		{
			s_resource[reservories][a][e] = 0.0;
			if(rp[2][a][e]<0)
			{
				correction = 0;
				for(i=1; i<reservories; i++)
				{
					s_resource[i][a][e] = (1 + rp[2][a][e])*s[i][a][e];
					s_resource[reservories][a][e] += s_resource[i][a][e];
					correction += s[i][a][e] - s_resource[i][a][e];
				}
				s_resource[0][a][e] = s[0][a][e] + correction;
				s_resource[reservories][a][e] += s_resource[0][a][e];
				for(i=reservories+1; i<states; i++)
					s_resource[i][a][e] = 0.0;
			}
			else
			{
		    	s_resource[0][a][e] = (1-rp[2][a][e])*s[0][a][e];
				s_resource[reservories][a][e] += s_resource[0][a][e];
				correction = s[0][a][e] - s_resource[0][a][e];
				sick = 0;
				for(i=1; i<reservories; i++)
					sick += s[i][a][e];
				if(sick==0)
					sick = 0.000001;
				for(i=1; i<reservories; i++)
				{
					s_resource[i][a][e] = s[i][a][e]*(1 + correction/sick);
					s_resource[reservories][a][e] += s_resource[i][a][e];
				}
				for(i=reservories+1; i<states; i++)
					s_resource[i][a][e] = 0.0;
			}
		}
	}

    //sprintf(aux_command, "mkdir -p lambda_iters/%s/", name_country);
    //system(aux_command);
    
	phase = 1;
	cont_clone = 0;
    load_variable_fitter(variable_parameters, vp);
	load_dynamic_variables(demo_global, s, S, year_min, year_end, 0);


    int long_save = EG*(1+AG*(reservories + internal_fluxes))*(year_end - year_min + 1);
    double save_output_to_R[long_save];
    for(a1=0; a1<long_save; a1++)
        save_output_to_R[a1] = 0.0;

	phase_temp = 1;
	err = model_run(demo_global, year_min, year_vac, save_output_to_R);
	if(err==1){
        *out_handler=1;
	}

	phase_temp = 2;
	err = model_run(demo_global, year_vac+1, year_end, save_output_to_R);
	if(err==1){
        *out_handler=1;
	} /*else{
        rewrite_output_object();
    }*/

    for(a1=0; a1<long_save; a1++){
        objeto1[a1] = save_output_to_R[a1];
        save_output_to_R[a1] = 0.0;
    }

    //Vaccine run starts here
    if(*vaccinate == 1){
        phase = 2;
        cont_clone = 0;
        load_dynamic_variables(demo_global, s, S, year_min, year_end, 0);
        phase_temp = 1;
        err = model_run(demo_global, year_min, year_vac, save_output_to_R);
        if(err==1){
            *out_handler=1;
        }

        char VacProfile[500];
        if(*adjust==0){
            sprintf(VacProfile, "../Inputs/General_Inputs/perfil.txt");
        }
        else{
            sprintf(VacProfile, "../Inputs_Iter/%s/Perfil/perfil_%03d.txt", name_country, temp_iter);
        }

        //load_vaccine_modifiers(VacProfile, variable_parameters, vp);
        load_vaccine_modifiers_H(VacProfile, variable_parameters, vp);

        phase_temp = 2;
        err = model_run(demo_global, year_vac+1, year_end, save_output_to_R);
        if(err==1){
            *out_handler=1;
        } /*else{
            rewrite_output_object();
        }*/

        for(a1=0; a1<long_save; a1++)
            objeto2[a1] = save_output_to_R[a1];
    }

}

int model_run(char filedem[string_length], int start, int end, double *objeto)
{
	int record, i, a, e, span;
	double t, taux;
	int err = 0;
	span = end-start + 1;
    double f_micro[AG][EG][internal_fluxes];

	t = start - year_min;
	printf("recording from %d to %d\n", start, end);
    int record_2 = start;

	for(record=0; record<span; record++) //Evolución hasta el máximo de años.
	{
		t_auxiliar=record;


        for(i=0; i<internal_fluxes; i++){
		    for(a=0; a<AG; a++){
			    for(e=0; e<EG; e++){
                 f_micro[a][e][i] = 0.0;
			    }
		    }
	    }

		for(taux=0;taux<record_window&&err==0;t+=H) //Evolución de un año al siguiente
		{

			err=rk4(gp, vp, rp, cm, s, S, f_micro, t);
			if (err==1)
				return 1;
			taux+=H;
		}

        //save_output_micro_file(s, f_micro, 0, record_2);
        //save_output_micro_file(s, f_micro, 1, record_2);
        //save_output_micro_file(s, f_micro, 2, record_2);

        save_output_micro(s, f_micro, 0, record_2, objeto);
        save_output_micro(s, f_micro, 1, record_2, objeto);
        save_output_micro(s, f_micro, 2, record_2, objeto);

        record_2++;

	}

	return err;
}

void save_output_micro(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int year, double *vv){
	int i, a, shift, index;
    int max_years = year_end - year_min + 1;
    int max_col_per_year = AG*reservories + 1 + internal_fluxes*AG;

    shift = (year - year_min)*max_col_per_year + e*max_col_per_year*max_years;
    vv[shift] = year;
	for(i=0; i<reservories; i++){
        for(a=0; a<AG; a++){
            index = shift + 1 + a + i*AG;
		    vv[index] = s[i][a][e];
        }
    }
    for(i=0; i<internal_fluxes; i++){
        for(a=0; a<AG; a++){
            index = shift + 1 + a + (i + reservories)*AG;
		    vv[index] = f[a][e][i];
        }
    }
}

void save_output_micro_file(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int start){
    FILE *fp1;
	char com_aux[1000];
    if (strcmp(name_country, simulation_name) != 0)
        sprintf(com_aux, "../Outputs/%s/Iter/%s/", name_country, simulation_name);
    else
	    sprintf(com_aux, "../Outputs/%s/", simulation_name);

    sprintf(com_aux, "%s/Output_object_%d_%d.txt", com_aux, phase, e);

	if(start==year_min){
		fp1=fopen(com_aux, "w");
	}
	else{
	    fp1=fopen(com_aux,"a");
	}

	int i, a;
    fprintf(fp1, "%d ", start);
	for(i=0; i<reservories; i++){
        for(a=0; a<AG; a++)
		    fprintf(fp1, "%.12lf ", s[i][a][e]);
    }
    for(i=0; i<internal_fluxes; i++){
        for(a=0; a<AG; a++)
		    fprintf(fp1, "%.12lf ", f[a][e][i]);
    }

	fprintf(fp1, "\n");
	fclose(fp1);
}

void rewrite_output_object(){
    FILE *fp1;
    FILE *fp2;
	char com_aux[1000], com[1000], sys_com[1000], line[30000], *RetChar;
    int E, i, systemRet;

    if (strcmp(name_country, simulation_name) != 0)
        sprintf(com_aux, "../Outputs/%s/Iter/%s/", name_country, simulation_name);
    else
	    sprintf(com_aux, "../Outputs/%s/", simulation_name);

    E = 0;
    sprintf(com, "%s/Output_object_%d_%d.txt", com_aux, phase, E);
	fp1 = fopen(com, "a");

    for(E=1; E<EG; E++){
        sprintf(com, "%s/Output_object_%d_%d.txt", com_aux, phase, E);
	    fp2 = fopen(com, "r");
        for(i = 0; i<(year_end - year_min + 1); i++){
            RetChar = fgets(line, 30000, fp2);
	        fprintf(fp1, "%s", line);
        }
        fclose(fp2);
        sprintf(sys_com, "rm %s", com);
        systemRet = system(sys_com);
    }
    fclose(fp1);
    sprintf(sys_com, "cp %s/Output_object_%d_%d.txt %s/Output_object_%d.txt", com_aux, phase, 0, com_aux, phase);
    systemRet = system(sys_com);
    sprintf(sys_com, "rm %s/Output_object_%d_%d.txt", com_aux, phase, 0);
    systemRet = system(sys_com);
}

int auxtemp=0;

int calculate_fluxes(int a, int e,double t,double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double f[fluxes],double S[states])
{
	double infect,infect_aux,pop,lambda;
	int i,a_prime,e_prime;
	int err=0;
	double beta,d,P_matrix[AG][AG];

	/* Calculating Scaled Infection Force and Diagnosis Rates */
	beta = funcion_param_t(rp[6][a][e],rp[7][a][e],t,2*rp[6][a][e],0);
	d = funcion_param_t(rp[3][a][e],rp[8][a][e],t,12.166667,0);

	f[9]=vp[24][a][e]*d*s[3][a][e];
	f[10]=vp[25][a][e]*rp[18][a][e]*d*s[4][a][e];
	f[11]=vp[26][a][e]*rp[18][a][e]*d*s[5][a][e];


	lambda=0;
	//Unscaled Infection force
	if(contact_type==1){
        if(cm_mode_all == 0){
            cm_M0(cm, P_matrix);
        }
        else if(cm_mode_all == 1){
            cm_M1(cm, P_matrix);
        }
        else if(cm_mode_all == 2){
            cm_M2(cm, P_matrix);
        }
        else if(cm_mode_all == 3){
            cm_M3(cm, P_matrix);
        }

		for(a_prime=0;a_prime<AG;a_prime++)
		{
			infect=0;
			pop=0;
			for(e_prime=0;e_prime<EG;e_prime++)
			{
				infect+=vp[19][a_prime][e_prime]*
				(
				 s[3][a_prime][e_prime]+
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*s[4][a_prime][e_prime]+
				 vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[9][a_prime][e_prime]+
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[11][a_prime][e_prime]+
				 s[18][a_prime][e_prime]
				);

				pop+=s[reservories][a_prime][e_prime];

			}
			lambda+=beta*P_matrix[a][a_prime]*infect/pop;
		}
	}
	else
	{
		infect=0;
		pop=0;

		for(a_prime=0;a_prime<AG;a_prime++)
			for(e_prime=0;e_prime<EG;e_prime++)
				pop+=s[reservories][a_prime][e_prime];

		for(a_prime=0;a_prime<AG;a_prime++)
		{
			for(e_prime=0;e_prime<EG;e_prime++)
			{
				infect_aux=
				(
				 s[3][a_prime][e_prime]+
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*s[4][a_prime][e_prime]+
				 vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[9][a_prime][e_prime]+
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[11][a_prime][e_prime]+
				 s[18][a_prime][e_prime]
				);

				infect+=infect_aux;
			}
		}
		lambda=beta*2*infect/pop;
	}
	
	/*FILE *out;
    char name_aux[1000];
    
	if(t>15 && t < 18){
        if(e==0){
            sprintf(name_aux, "lambda_iters/%s/force_of_infection_%s_%03d.txt", name_country, name_country, temp_iter);
            if((temp_iter==0)&&(auxtemp==0)){
                out = fopen(name_aux, "wt");
                auxtemp=1;
            }
            else{
                out = fopen(name_aux, "a");
            }
            fprintf(out, "%d %.12lf\n", a, lambda);
            fclose(out);
        }
    }*/
    

	/* Infection */
	f[1]=lambda*vp[6][a][e]*vp[10][a][e]*rp[9][a][e]*s[0][a][e]; //To fast-progression
	f[2]=lambda*vp[6][a][e]*(1-vp[10][a][e]*rp[9][a][e])*s[0][a][e]; //To slow-progression

	/*Fast-Progression*/
	f[3]=vp[12][a][e]*gp[0][a][e]*rp[10][a][e]*s[1][a][e]; //To smear-positive
	f[4]=vp[12][a][e]*gp[0][a][e]*(1-rp[10][a][e]-rp[11][a][e])*s[1][a][e]; //To smear-negative
	f[5]=vp[12][a][e]*gp[0][a][e]*rp[11][a][e]*s[1][a][e]; //To extra-pulmonary

	/*Slow-Progression*/
	f[6]=vp[13][a][e]*gp[1][a][e]*rp[10][a][e]*s[2][a][e]; //To smear-positive
	f[7]=vp[13][a][e]*gp[1][a][e]*(1-rp[10][a][e]-rp[11][a][e])*s[2][a][e]; //To smear-negative
	f[8]=vp[13][a][e]*gp[1][a][e]*rp[11][a][e]*s[2][a][e]; //To extra-pulmonary

	/*Natural Recovery */
	f[12]=vp[18][a][e]*gp[10][a][e]*s[3][a][e]; //In smear-positives
	f[13]=vp[18][a][e]*gp[10][a][e]*s[4][a][e]; //In smear-negatives
	f[14]=vp[18][a][e]*gp[10][a][e]*s[5][a][e]; //In extra-pulmonary

	/*Treatment*/
	f[15]=gp[6][a][e]*rp[13][a][e]*s[6][a][e]; //smear-positive: default
	f[16]=gp[6][a][e]*rp[16][a][e]*s[7][a][e]; //smear-negative: default
	f[17]=gp[6][a][e]*rp[16][a][e]*s[8][a][e]; //extra-pulmonary: default
	f[18]=gp[6][a][e]*(1-rp[12][a][e]-rp[13][a][e]-rp[14][a][e])*s[6][a][e]; //smear-positive: success
	f[19]=gp[6][a][e]*(1-rp[15][a][e]-rp[16][a][e]-rp[17][a][e])*s[7][a][e]; //smear-negative: success
	f[20]=gp[6][a][e]*(1-rp[15][a][e]-rp[16][a][e]-rp[17][a][e])*s[8][a][e]; //extra-pulmonary: default

	/*Smear-progression*/
	f[21]=vp[19][a][e]*gp[7][a][e]*s[4][a][e]; //Undiagnosed
	f[22]=vp[19][a][e]*gp[7][a][e]*s[7][a][e]; //Diagnosed (in treatment)

	/*Relapses*/
	f[23]=vp[20][a][e]*gp[8][a][e]*s[9][a][e]; //Default, smear-positive
	f[24]=vp[20][a][e]*gp[8][a][e]*s[11][a][e]; //Default, smear-negative
	f[25]=vp[20][a][e]*gp[8][a][e]*s[13][a][e]; //Default, extra-pulmonary
	f[26]=vp[21][a][e]*gp[9][a][e]*s[10][a][e]; //Natural, smear-positive
	f[27]=vp[21][a][e]*gp[9][a][e]*s[12][a][e]; //Natural, smear-negative
	f[28]=vp[21][a][e]*gp[9][a][e]*s[14][a][e]; //Natural, smear-extra-pulmonary

	/*Reinfections*/
	for(i=0;i<6;i++){
		f[29+i] = vp[6][a][e]*vp[10][a][e]*rp[9][a][e]*vp[11][a][e]*gp[5][a][e]*lambda*s[9+i][a][e];	//After treatment
        //f[29+i] = vp[6][a][e]*(1.0-vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[9+i][a][e];
    }
	f[35] = vp[6][a][e]*vp[10][a][e]*rp[9][a][e]*vp[11][a][e]*gp[5][a][e]*lambda*s[2][a][e]; //Before treatment (from slow to fast latency)

	/*Deaths by TB*/
	f[36] = vp[7][a][e]*gp[2][a][e]*s[3][a][e]; //Smear-positive
	f[37] = vp[8][a][e]*gp[3][a][e]*s[4][a][e]; //Smear-negative
	f[38] = vp[9][a][e]*gp[4][a][e]*s[5][a][e]; //Extra-pulmonary

	if(dynamics_type < 3)
	{
		int i;
		f[0]=0.0;//births and not TB deaths are not directly calculated, but from the database
		for(i=39;i<56;i++)
			f[i]=0.0;
	}
	else
	{
		printf("algorithm not ready for dyn_type>2\n");
		exit(1);
	}

	/*Treatment (Second part)*/
	f[56]=rp[12][a][e]*gp[6][a][e]*s[6][a][e]; //Failure, smear-positive
	f[57]=gp[6][a][e]*rp[14][a][e]*s[6][a][e]; //Death, smear-positive
	f[58]=gp[6][a][e]*rp[17][a][e]*s[7][a][e]; //Death, smear-negative
	f[59]=gp[6][a][e]*rp[17][a][e]*s[8][a][e]; //Death, extra-pulmonary

	f[60]=gp[2][a][e]*s[18][a][e]; //Death of failure

	/* Relapses from success */
	f[61]=vp[22][a][e]*gp[14][a][e]*s[15][a][e]; //Smear-positive
	f[62]=vp[22][a][e]*gp[14][a][e]*s[16][a][e]; //Smear-negative
	f[63]=vp[22][a][e]*gp[14][a][e]*s[17][a][e]; //Extra-pulmonary

	/*Reinfections from success to fast latency*/
	for(i=0;i<3;i++){
		f[64+i] = vp[6][a][e]*vp[10][a][e]*rp[9][a][e]*vp[11][a][e]*gp[5][a][e]*lambda*s[15+i][a][e];
        //f[64+i] = 0.0;
    }

	/*Reinfections from success to slow latency (not considered)*/
	for(i=0;i<3;i++){
        f[67+i] = 0;
        //f[67+i] = vp[6][a][e]*(1 - vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[15+i][a][e];
	}

	/*Treatment (Third part)*/
	f[70]=rp[15][a][e]*gp[6][a][e]*s[7][a][e]; //Failure, smear-negative
	f[71]=rp[15][a][e]*gp[6][a][e]*s[8][a][e]; //Failure, extra-pulmonary

	/*Reinfections followed by slow-progression: Do not have an effect in the dynamics*/
	for(i=0;i<6;i++)
		f[72+i]=vp[6][a][e]*(1.0-vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[9+i][a][e];	//After treatment
	f[78] = vp[6][a][e]*(1.0-vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[2][a][e]; //Before treatment (from slow to fast latency)
	for(i=0;i<3;i++)
		f[79+i]=vp[6][a][e]*(1.0-vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[15+i][a][e];

	/*Potential infection*/
	/*Infections we'll have if all the population was susceptible*/
	f[82]=0.0; //vp[6][a][e]*lambda*s[reservories][a][e];

	int ii;
	for(ii=0;ii<internal_fluxes;ii++)
	{
		if(f[ii]<0)
		{
			printf("negative flux\n");
			printf("rp[9][a][e] %lf s[0][a][e] %lf,(a %d e %d t %lf) flux[%d]=%.16lf lambda %.16lf (%.16lf %lf %lf)(%.16lf %lf %lf)\n",rp[9][a][e],s[0][a][e],a,e,t,ii,f[ii],lambda,rp[3][a][e],rp[4][a][e],rp[5][a][e],rp[6][a][e],rp[7][a][e],rp[8][a][e]);
			err=1;
		}
	}
	return err;
}

int deriv(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG], double ds[states][AG][EG], double f_micro[AG][EG][internal_fluxes], double t,double S[states])
{
	int i, a, e;
	double f[fluxes];
	int err=0;

	for(a=0;a<AG;a++)
	{
		for(e=0;e<EG;e++)
		{
			for(i=0;i<states;i++)
			{
				//if(s[i][a][e]<0){
				//	printf("t=%lf negative value in deriv: s[%d][%d][%d]=%lf\n", t, i, a, e, s[i][a][e]);
                //}
			}
			err=calculate_fluxes(a, e, t, gp, vp, rp, cm, s, f, S);
			calculate_derivatives(a, e, ds, f);
            for(i=0; i<internal_fluxes; i++){
                f_micro[a][e][i] = f[i];
            }
			if(err==1)
				return(err);
		}
	}
	ageing(ds, s, S, vp, t);
	return err;
}

int rk4(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG], double S[states], double f_micro[AG][EG][internal_fluxes], double t)
{
	int i,a,e;
	double ds[states][AG][EG], dsm[states][AG][EG], dst[states][AG][EG], st[states][AG][EG],St[states];
    double f_m1[AG][EG][internal_fluxes], f_m2[AG][EG][internal_fluxes], f_m3[AG][EG][internal_fluxes], fs[AG][EG][internal_fluxes];
	double hh,h6;
	int err=0;
	hh=H*0.5;
	h6=H/6.0;

	err=deriv(gp, vp, rp, cm, s, ds, fs, t, S);                /*ds=f(x0)*/
	if(err==1)
	{
		printf("Runge-Kutta fails in first call to deriv\n");
		return err;
	}

	for(i=0;i<states;i++)
	{
		St[i]=0;
		for(a=0;a<AG;a++)
		{
			for(e=0;e<EG;e++)
			{
				st[i][a][e] = s[i][a][e] + hh*ds[i][a][e];
				St[i] += st[i][a][e];
			}
		}
	}

	err=deriv(gp, vp, rp, cm, st, dst, f_m1, t+hh, St);           /*st es x1, dst=f(x1)*/
	if(err==1)
		return err;

	for(i=0;i<states;i++)
	{
		St[i]=0;
		for(a=0;a<AG;a++)
		{
			for(e=0;e<EG;e++)
			{
				st[i][a][e] = s[i][a][e] + hh*dst[i][a][e];
				St[i] += st[i][a][e];
			}
		}
	}

	err=deriv(gp, vp, rp, cm, st, dsm, f_m2, t+hh, St);           /*xt es x2, dxm=f(x2)*/
	if(err==1)
		return err;

	for(i=0;i<states;i++)
	{
		St[i]=0;
		for(a=0;a<AG;a++)
		{
			for(e=0;e<EG;e++)
			{
				st[i][a][e] = s[i][a][e] + H*dsm[i][a][e];  /*st es x3, dsm=f(x1)+f(x2)*/
				St[i] += st[i][a][e];
				dsm[i][a][e] += dst[i][a][e];
			}
		}
	}

	err=deriv(gp, vp, rp, cm, st, dst, f_m3, t+H, St);           /*dst=f(x3)*/
	if(err==1)
		return err;

	for(i=0;i<states;i++)
	{
		S[i]=0;
		for(a=0;a<AG;a++)
		{
			for(e=0;e<EG;e++)
			{
				s[i][a][e] += h6*(ds[i][a][e] + dst[i][a][e] + 2.0*dsm[i][a][e]);  /*st es x3, dsm=f(x1)+f(x2)*/
				S[i] += s[i][a][e];
			}
		}
	}
    for(i=0; i<internal_fluxes; i++){
		for(a=0; a<AG; a++){
			for(e=0; e<EG; e++){
                f_micro[a][e][i] += h6*(fs[a][e][i] + f_m3[a][e][i] + 2*f_m1[a][e][i] + 2*f_m2[a][e][i] );
			}
		}
	}

	return err;
}

void load_initial_conditions(double s[states][AG][EG], double S[states])
{
	FILE *fich;
	int entries;
	int i,a,e;
	int amin,amax,emin,emax,index;
	double value,vmin,vmax;
	char name[256];
	char file_name[string_length];

	sprintf(file_name,"%s/Initial_Conditions/initial_conditions.txt",path_fitted_inputs);
	fich=fopen(file_name,"rt");

	scan=fscanf(fich,"entries = %d \n",&entries);
	scan=fscanf(fich,"description \t index \t amin \t amax \t emin \t emax \t value \t interval \n");

	for(i=0;i<entries;i++)
	{
		scan=fscanf(fich,"%s \t %d \t %d \t %d \t %d \t %d \t %lf \t %lf \t %lf \n",name,&index,&amin,&amax,&emin,&emax,&value,&vmin,&vmax);
		if(index==reservories) continue;//We dont load the populations here, that is done from the demo. pyramid:
		//Reason is that we might want to run intervals considering, for example, demographic uncertainty, but not uncertainty from initial conditions
		//So we load S[reservories] (which is total population) from the demo. pyr., and the rest of states from here (ini. cond.)
		//Always in that order, first demo. pyr. then initial conditions.
		amin=(int)(amin/age_window);
		amax=(int)(amax/age_window);

		if(amax>AG&&amin<=AG) amax=AG;
		if(amin<0||amin>AG||amax<0||amax>AG||emin<0||emin>EG||emax<0||emax>EG||amin>amax||emin>emax)
		{
			printf("error de rango en la lectura de las condiciones iniciales\n");
			exit(1);
		}
		for (a=amin;a<amax;a++)
		{
			for(e=emin;e<emax;e++)
			{
				s[index][a][e]=value;
				if(index>reservories)
					s[index][a][e]=0.0;//If for any reason you pass me obsevables in the file (you shouldnt), I read them, but immediately put them at zero
			}
		}
	}

	double auxiliar;
	for(i=0;i<states;i++)//I have to calculate S[reservories] (age-aggregated population) here, it's not done when I load demography
		S[i]=0;
	for (a=0;a<AG;a++)
	{
		for(e=0;e<EG;e++)
		{
			auxiliar=0;
			for(i=1;i<states;i++)
			{
				if(i<reservories)
				{
					s[i][a][e]*=s[reservories][a][e];
					auxiliar+=s[i][a][e];
					S[i]+=s[i][a][e];
				}
				else
				{
					if(i>reservories)
						s[i][a][e]=0.0;
					S[i]+=s[i][a][e];
				}
			}
			s[0][a][e]=s[reservories][a][e]-auxiliar;
			if(s[0][a][e]<0)
			{
				printf("error loading initial conditions: negative class (a=%d e=%d) s[0][a][e]=%f\n",a,e,s[0][a][e]);
				exit(1);
			}
			S[0]+=s[0][a][e];
		}
	}
	fclose(fich);
}

void ageing(double ds[states][AG][EG], double s[states][AG][EG],double S[states],double vp[variable_parameters][AG][EG],double t){
	int i, a, e;
	double sick_density, numerator, denominator, vac_flux[reservories+1][EG-1], flux, flux_aux;
    
    double flux_AON = 0.0;

	flux = ds[0][0][1];
	flux_aux = ds[0][1][1];

	//Newborns
	for(i=0; i<=reservories; i++)
	{
		for(e=0; e<EG; e++)
			ds[i][0][e] -= s[i][0][e]/(double)age_window;
	}

    //Vaccination. Rutine over ageing.
	for(a=1; a<AG; a++){
        flux_AON = 0.0;
		for(e=0; e<EG-1; e++){
            vac_flux[0][e] = vp[1][a][e]*s[0][a-1][e];
			vac_flux[1][e] = vp[2][a][e]*s[1][a-1][e];
			vac_flux[2][e] = vp[2][a][e]*s[2][a-1][e];
            
            if(eps_AON==1 && e==0 && a==age_AON){
                flux_AON = cov_AON*(s[1][a-1][0])/(double)age_window; 
                //printf("%lf\n", flux_AON);
            }

			for(i=3; i<reservories; i++){
				vac_flux[i][e] = 0.0;
            }
			for(i=9; i<18; i++){
                vac_flux[i][e] = vp[28][a][e]*s[i][a-1][e];
            }
			vac_flux[reservories][e] = 0.0;
			for(i=0; i<reservories; i++){
				vac_flux[reservories][e] += vac_flux[i][e];
            }
		}

		for(i=0; i<=reservories; i++){
            if(eps_AON==1){
                if(i==1){
                    ds[1][a][0] += (s[i][a-1][0] - s[i][a][0] - vac_flux[1][0])/(double)age_window; //- flux_AON;
                    ds[1][a][1] += (s[i][a-1][1] - s[i][a][1])/(double)age_window;
                    ds[1][a][2] += (s[i][a-1][2] - s[i][a][2])/(double)age_window;
                }
                else if(i==2){
                    ds[2][a][0] += (s[i][a-1][0] - s[i][a][0] - vac_flux[2][0])/(double)age_window;// + flux_AON;
                    ds[2][a][1] += (s[i][a-1][1] - s[i][a][1] + vac_flux[2][0] + vac_flux[1][0])/(double)age_window;
                    ds[2][a][2] += (s[i][a-1][2] - s[i][a][2])/(double)age_window;
                }
                else{
                    ds[i][a][0] += (s[i][a-1][0] - s[i][a][0] - vac_flux[i][0])/(double)age_window;
                    ds[i][a][1] += (s[i][a-1][1] - s[i][a][1] + vac_flux[i][0])/(double)age_window;
                    ds[i][a][2] += (s[i][a-1][2] - s[i][a][2])/(double)age_window;
                }
            }
            else{
                ds[i][a][0] += (s[i][a-1][0] - s[i][a][0] - vac_flux[i][0])/(double)age_window;
                ds[i][a][1] += (s[i][a-1][1] - s[i][a][1] + vac_flux[i][0])/(double)age_window;
                ds[i][a][2] += (s[i][a-1][2] - s[i][a][2])/(double)age_window;
            }

    	}
		for(i=0; i<reservories; i++)
			ds[reservories+4][a][0] += vac_flux[i][0]/(double)age_window;
	}

	double delta;

	if(dynamics_type==0)
	{
		if((phase==1)||(phase==3))//phase=1 es run sin intervención (sin vacunacion)
		{
			numerator = 0;
			denominator = 0;
			for(a=start_maternity; a<=stop_maternity; a++)
			{
				for(e=0; e<EG; e++)
				{
					numerator += (s[3][a][e] + s[4][a][e] + s[5][a][e] + s[9][a][e] + s[11][a][e] + s[13][a][e] + s[18][a][e])*gp[13][a][e]*vp[23][a][e];
					for(i=0; i<reservories; i++)
						denominator += s[i][a][e];
				}
			}

			if(denominator!=0)
				sick_density = numerator/denominator;
			else
				sick_density = 0;
			delta = 0;
			for(e=0; e<EG; e++)
			{
				for(a=0; a<AG; a++)
				{
					for(i=0; i<reservories; i++)
					{
						delta -= ds[i][a][e];//Total increment of population
					}
				}
			}

			for(e=0; e<EG; e++)
			{
				ds[0][0][e] += (vp[0][0][e])*delta*(1 - sick_density);
				ds[1][0][e] += (vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
				ds[2][0][e] += (vp[0][0][e])*delta*(sick_density)*(1 - rp[9][0][e]*vp[10][a][e]);
				if(mother_child_infections_accounted==1)
					ds[reservories+5][0][e] += (vp[0][0][e])*delta*(sick_density);
				ds[reservories][0][e] += delta*vp[0][0][e];
				b[e][cont_clone] = delta/S[reservories];
			}
			ds[reservories+4][0][0]+=(vp[0][0][1])*delta;
			cont_clone++;

		}
		else
		{
			if((phase==2)||(phase==4))
			{
				numerator = 0;
				denominator = 0;
				for(a=start_maternity; a<=stop_maternity; a++)
				{
					for(e=0; e<EG; e++)
					{
						numerator += (s[3][a][e] + s[4][a][e] + s[5][a][e] + s[9][a][e] + s[11][a][e] + s[13][a][e] + s[18][a][e])*gp[13][a][e]*vp[23][a][e];
						for(i=0; i<reservories; i++)
							denominator += s[i][a][e];
					}
				}

				if(denominator!=0)
					sick_density = numerator/denominator;
				else
					sick_density = 0;
				for(e=0; e<EG; e++)
				{
					delta = b[e][cont_clone]*S[reservories];
					ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
					ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
					ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
					ds[reservories][0][e]+=delta*vp[0][0][e];
					if(mother_child_infections_accounted==1)
						ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
				}
				cont_clone++;
				ds[reservories+4][0][0]+=(vp[0][0][1])*delta;

			}
			else
			{
				if((phase==0)||(phase==5))
				{
					numerator=0;
					denominator=0;
					for(a=start_maternity;a<=stop_maternity;a++)
					{
						for(e=0;e<EG;e++)
						{
							numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[23][a][e];
							for(i=0;i<reservories;i++)
								denominator+=s[i][a][e];
						}
					}

					if(denominator!=0)
						sick_density=numerator/denominator;
					else
						sick_density=0;
					delta=0;
					for(e=0;e<EG;e++)
					{
						for(a=0;a<AG;a++)
						{
							for(i=0;i<reservories;i++)
							{
								delta-=ds[i][a][e];
							}
						}
					}

					for(e=0;e<EG;e++)
					{
						ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
						ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
						ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
						if(mother_child_infections_accounted==1)
							ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
						ds[reservories][0][e]+=delta*vp[0][0][e];
					}
					cont_clone++;
					ds[reservories+4][0][0]+=(vp[0][0][1])*delta;
				}
				else
				{
					printf("error with phase\n");
					exit(1);
				}
			}
		}
	}
	else
	{
		if (dynamics_type==1)
		{
			if((phase==1)||(phase==3))
			{

				numerator=0;
				denominator=0;
				for(a=start_maternity;a<=stop_maternity;a++)
				{
					for(e=0;e<EG;e++)
					{
						numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[23][a][e];
						for(i=0;i<reservories;i++)
							denominator+=s[i][a][e];
					}
				}

				if(denominator!=0)
					sick_density=numerator/denominator;
				else
					sick_density=0;
				delta=0;
				for(e=0;e<EG;e++)
				{
					for(i=0;i<reservories;i++)
						delta-=ds[i][0][e];
				}

				for(e=0;e<EG;e++)
				{
					ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
					ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
					ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
					ds[reservories][0][e]+=delta*vp[0][0][e];
					if(mother_child_infections_accounted==1)
						ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
					b[e][cont_clone]=delta/S[reservories];
				}
				ds[reservories+4][0][0]+=(vp[0][0][1])*delta;

				for(a=1;a<AG;a++)
				{
					delta=0;
					denominator=0;

					for(e=0;e<EG;e++)
					{
						for(i=0;i<reservories;i++)
						{
							delta-=ds[i][a][e];
							denominator+=s[i][a][e];
						}
					}

					for(e=0;e<EG;e++)
					{
						for(i=0;i<=reservories;i++)
							ds[i][a][e]+=delta*(s[i][a][e]/denominator);
						mu[a][e][cont_clone]=delta/denominator;
					}
				}
				cont_clone++;
			}
			else
			{
				if((phase==2)||(phase==4))
				{

					numerator=0;
					denominator=0;
					for(a=start_maternity;a<=stop_maternity;a++)
					{
						for(e=0;e<EG;e++)
						{
							numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[23][a][e];
							for(i=0;i<reservories;i++)
								denominator+=s[i][a][e];
						}
					}

					if(denominator!=0)
						sick_density=numerator/denominator;
					else
						sick_density=0;


					for(e=0;e<EG;e++)
					{
						delta=b[e][cont_clone]*S[reservories];
						ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
						ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
						ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
						if(mother_child_infections_accounted==1)
							ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
						ds[reservories][0][e]+=delta*vp[0][0][e];
					}
					ds[reservories+4][0][0]+=(vp[0][0][1])*delta;
					for(a=1;a<AG;a++)
					{
						delta=0;
						denominator=0;

						for(e=0;e<EG;e++)
						{
							for(i=0;i<reservories;i++)
							{
								delta-=ds[i][a][e];
								denominator+=s[i][a][e];
							}
						}

						for(e=0;e<EG;e++)
						{
							delta=mu[a][e][cont_clone]*denominator;
							for(i=0;i<=reservories;i++)
								ds[i][a][e]+=delta*(s[i][a][e]/denominator);
						}
					}
					cont_clone++;
				}
				else
				{
					if((phase==0)||(phase==5))
					{

						numerator=0;
						denominator=0;
						for(a=start_maternity;a<=stop_maternity;a++)
						{
							for(e=0;e<EG;e++)
							{
								numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[23][a][e];
								for(i=0;i<reservories;i++)
									denominator+=s[i][a][e];
							}
						}

						if(denominator!=0)
							sick_density=numerator/denominator;
						else
							sick_density=0;
						delta=0;
						for(e=0;e<EG;e++)
						{
							for(i=0;i<reservories;i++)
								delta-=ds[i][0][e];
						}

						for(e=0;e<EG;e++)
						{
							ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
							ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
							ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
							if(mother_child_infections_accounted==1)
								ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
							ds[reservories][0][e]+=delta*vp[0][0][e];
						}
						ds[reservories+4][0][0]+=(vp[0][0][1])*delta;
						for(a=1;a<AG;a++)
						{
							delta=0;
							denominator=0;
							for(e=0;e<EG;e++)
							{
								for(i=0;i<reservories;i++)
								{
									delta-=ds[i][a][e];
									denominator+=s[i][a][e];
								}
							}
							for(e=0;e<EG;e++)
							{
								for(i=0;i<=reservories;i++)
									ds[i][a][e]+=delta*(s[i][a][e]/denominator);
							}
						}
					}
					else
					{
						printf("error with phase bis\n");
						exit(1);
					}
				}
			}
		}
		else
		{
			if (dynamics_type==2)
			{
				if((phase==1)||(phase==3))
				{
					numerator=0;
					denominator=0;
					for(a=start_maternity;a<=stop_maternity;a++)
					{
						for(e=0;e<EG;e++)
						{
							numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[23][a][e];
							for(i=0;i<reservories;i++)
								denominator+=s[i][a][e];
						}
					}

					if(denominator!=0)
						sick_density=numerator/denominator;
					else
						sick_density=0;
					delta=0;
					for(e=0;e<EG;e++)
					{
						for(i=0;i<reservories;i++)
							delta-=ds[i][0][e];
					}
					delta+=semiempiric_derivative_fitted_population(0,0,t);
					for(e=0;e<EG;e++)
					{
						ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
						ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
						ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
						if(mother_child_infections_accounted==1)
							ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
						ds[reservories][0][e]+=delta*vp[0][0][e];
						b[e][cont_clone]=delta/S[reservories];
					}
					ds[reservories+4][0][0]+=(vp[0][0][1])*delta;
					for(a=1;a<AG;a++)
					{
						delta=0;
						denominator=0;
						for(e=0;e<EG;e++)
						{
							for(i=0;i<reservories;i++)
							{
								delta-=ds[i][a][e];
								denominator+=s[i][a][e];
							}
						}
						delta+=semiempiric_derivative_fitted_population(a,0,t);

						for(e=0;e<EG;e++)
						{
							for(i=0;i<=reservories;i++)
								ds[i][a][e]+=delta*(s[i][a][e]/denominator);
							mu[a][e][cont_clone]=delta/denominator;
						}
					}
					cont_clone++;
				}
				else
				{
					if((phase==2)||(phase==4))
					{
						numerator=0;
						denominator=0;
						for(a=start_maternity;a<=stop_maternity;a++)
						{
							for(e=0;e<EG;e++)
							{
								numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[23][a][e];
								for(i=0;i<reservories;i++)
									denominator+=s[i][a][e];
							}
						}

						if(denominator!=0)
							sick_density=numerator/denominator;
						else
							sick_density=0;


						for(e=0;e<EG;e++)
						{
							delta=b[e][cont_clone]*S[reservories];
							ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
							ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
							ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
							if(mother_child_infections_accounted==1)
								ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
							ds[reservories][0][e]+=delta*vp[0][0][e];
						}
						ds[reservories+4][0][0]+=(vp[0][0][1])*delta;
						for(a=1;a<AG;a++)
						{
							delta=0;
							denominator=0;

							for(e=0;e<EG;e++)
							{
								for(i=0;i<reservories;i++)
								{
									delta-=ds[i][a][e];
									denominator+=s[i][a][e];
								}
							}

							for(e=0;e<EG;e++)
							{
								delta=mu[a][e][cont_clone]*denominator;
								for(i=0;i<=reservories;i++)
									ds[i][a][e]+=delta*(s[i][a][e]/denominator);
							}
						}
						cont_clone++;
					}
					else
					{
						if((phase==0)||(phase==5))
						{

							numerator=0;
							denominator=0;
							for(a=start_maternity;a<=stop_maternity;a++)
							{
								for(e=0;e<EG;e++)
								{
									numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[23][a][e];
									for(i=0;i<reservories;i++)
										denominator+=s[i][a][e];
								}
							}
							if(denominator!=0)
								sick_density=numerator/denominator;
							else
								sick_density=0;
							delta=0;
							for(e=0;e<EG;e++)
							{
								for(i=0;i<reservories;i++)
									delta-=ds[i][0][e];
							}
							delta+=semiempiric_derivative_fitted_population(0,0,t);
							for(e=0;e<EG;e++)
							{
								ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
								ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
								ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
								ds[reservories][0][e]+=delta*vp[0][0][e];
								if(mother_child_infections_accounted==1)
									ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
							}
							ds[reservories+4][0][0]+=(vp[0][0][1])*delta;
							for(a=1;a<AG;a++)
							{
								delta=0;
								denominator=0;
								for(e=0;e<EG;e++)
								{
									for(i=0;i<reservories;i++)
									{
										delta-=ds[i][a][e];
										denominator+=s[i][a][e];
									}
								}
								delta+=semiempiric_derivative_fitted_population(a,0,t);

								for(e=0;e<EG;e++)
								{
									for(i=0;i<=reservories;i++)
										ds[i][a][e]+=delta*(s[i][a][e]/denominator);
								}
							}
						}
						else
						{
							printf("error with phase bis\n");
							exit(1);
						}
					}
				}
			}
			else
			{
				printf("error with dynamic type\n");
				exit(1);
			}
		}
	}
}

double semiempiric_derivative_fitted_population(int a,int e,double t)
{
	int i;
	double result;
	result=0;
	for(i=1;i<=fitting_degree;i++)
		result+=pyramid_coefficient[a][e][i]*i*pow(t,(i-1));
	return(result);
}

void load_fit_pyramids(double s[states][AG][EG],char pyr[string_length],int start,int end){

	FILE *fich;
	double pyramid[AG][EG][demo_span_max];
	int entries,span;
	int i,a,e,t,j,k;
	int amin,amax,emin,emax,tmin,tmax;
	double value,vmin,vmax;
	span=end-start+1;

	if(span>demo_span_max)
	{
		printf("Too many years selected when loading pyramids\n");
		exit(1);
	}

	fich=fopen(pyr,"rt");

	scan=fscanf(fich,"entries=%d \n",&entries);
	scan=fscanf(fich,"amin \t amax \t emin \t emax \t tmin \t tmax \t value \t interval\n");

	for(i=0;i<entries;i++)
	{
		scan=fscanf(fich,"%d \t %d \t %d \t %d \t %d \t %d \t %lf \t %lf \t %lf\n",&amin,&amax,&emin,&emax,&tmin,&tmax,&value,&vmin,&vmax);

		amin=(int)(amin/age_window);
		amax=(int)(amax/age_window);

		if(amax>AG&&amin<=AG) amax=AG;

		if(tmin<year_min) continue;

		if(dynamics_type!=2)
		{
			//If dynamic type is not 2, only read the inicial pyramid, and then consider it constant
			if(tmin==start)
				tmax=end+1;
			else
				continue;
		}
		for (a=amin;a<amax;a++)
		{
			for(e=emin;e<emax;e++)
			{
				for(t=tmin;t<tmax;t++)
					pyramid[a][e][t-year_min]=value;
			}
		}
	}
	fclose(fich);

	if(phase==0)
	{
		for(j=0;j<AG;j++)
		{
			for(k=0;k<EG;k++)
			{
				for(i=0;i<fitting_degree+1;i++)
					pyramid_coefficient[j][k][i]=pyramid_coefficient_aux[j][k][i][0];
			}
		}

	}
	else
	{
		for(j=0;j<AG;j++)
		{
			for(k=0;k<EG;k++)
			{
				for(i=0;i<fitting_degree+1;i++)
					pyramid_coefficient[j][k][i]=pyramid_coefficient_aux[j][k][i][1];
			}
		}
	}

    for(a=0;a<AG;a++){
		s[reservories][a][0] = pyramid[a][0][start-year_min];//Initial Populations
        s[reservories][a][1] = 0.0;
        s[reservories][a][2] = 0.0;
    }
}

void load_dynamic_variables(char demo_r[100], double s[states][AG][EG], double S[states], int start, int end, int is_file)
{
	load_fit_pyramids(s, demo_r, start, end);
	if(is_file==1)
		load_initial_conditions(s, S);
	else
		load_initial_conditions_from_s_resource(s, S);
}

void load_initial_conditions_from_s_resource(double s[states][AG][EG], double S[states])
{
	int i, a, e;

	for(i=0; i<states; i++)
	{
		S[i] = 0.0;
		for(a=0; a<AG; a++)
		{
			for(e=0; e<EG; e++)
			{
				s[i][a][e] = s_resource[i][a][e];
				S[i] += s[i][a][e];
			}
		}
	}
}

void read_pyr_coef()
{
	FILE *fp_local;
	char file_name[string_length];
	int i, j, j_aux, k_aux;

	sprintf(file_name, "%s/Pyr_coef_fit.txt", path_fitted_inputs);
	fp_local=fopen(file_name, "rt");
	for(j=0; j<AG; j++)
	{
		//for(k=0; k<EG; k++)
		//{
			scan = fscanf(fp_local, "%d %d ", &j_aux, &k_aux);
			for(i=0; i<fitting_degree+1; i++){
				scan=fscanf(fp_local,"%lf ", &pyramid_coefficient_aux[j][0][i][0]);
                pyramid_coefficient_aux[j][1][i][0] = 0.0;
                pyramid_coefficient_aux[j][2][i][0] = 0.0;
            }
			scan=fscanf(fp_local,"\n");
		//}
	}
	fclose(fp_local);

	sprintf(file_name, "%s/Pyr_coef_run.txt", path_fitted_inputs);
	fp_local = fopen(file_name, "rt");
	for(j=0; j<AG; j++)
	{
		//for(k=0; k<EG; k++)
		//{
			scan=fscanf(fp_local, "%d %d ", &j_aux, &k_aux);
			for(i=0; i<fitting_degree+1; i++){
				scan = fscanf(fp_local, "%lf ", &pyramid_coefficient_aux[j][0][i][1]);
                pyramid_coefficient_aux[j][1][i][1] = 0.0; //k itera donde ahora hay [j][k][i][r], es el segundo índice
                pyramid_coefficient_aux[j][2][i][1] = 0.0;
            }
			scan = fscanf(fp_local, "\n");
		//}
	}
	fclose(fp_local);

}

void write_output_run_inf_matrix(int uncertainty, int phase,int iter,double inf[AG][AG], double inf_1[AG][AG], double inf_2[AG][AG],double inf_3[AG][AG])
{
	FILE *fp1;
	int a1,a2;
	char com_aux[1000],path_def[1000],out[1000];
    if (strcmp(name_country, simulation_name) != 0)
        sprintf(com_aux, "../Outputs/%s/Iter/%s/", name_country, simulation_name);
    else
	    sprintf(com_aux, "../Outputs/%s/", simulation_name);
	sprintf(path_def,"%s/macro",com_aux);

	sprintf(com_aux,"mkdir -p %s",path_def);
	scan=system(com_aux);

	sprintf(out,"%s/output_run_potential_infections_matrix_%d_phase_%d.txt",path_def,uncertainty>0,phase);
	if(iter==0)
		fp1=fopen(out,"w");
	else
		fp1=fopen(out,"a");
	for(a1=0;a1<AG;a1++)
	{
		for(a2=0;a2<AG;a2++)
			fprintf(fp1,"%lf ",inf[a1][a2]);
		fprintf(fp1,"\n");
	}
	fprintf(fp1,"\n");
	fclose(fp1);

	sprintf(out,"%s/output_run_primary_infections_matrix_%d_phase_%d.txt",path_def,uncertainty>0,phase);
	if(iter==0)
		fp1=fopen(out,"w");
	else
		fp1=fopen(out,"a");
	for(a1=0;a1<AG;a1++)
	{
		for(a2=0;a2<AG;a2++)
			fprintf(fp1,"%lf ",inf_1[a1][a2]);
		fprintf(fp1,"\n");
	}
	fprintf(fp1,"\n");
	fclose(fp1);

	sprintf(out,"%s/output_run_re_infections+fast_matrix_%d_phase_%d.txt",path_def,uncertainty>0,phase);
	if(iter==0)
		fp1=fopen(out,"w");
	else
		fp1=fopen(out,"a");
	for(a1=0;a1<AG;a1++)
	{
		for(a2=0;a2<AG;a2++)
			fprintf(fp1,"%lf ",inf_2[a1][a2]);
		fprintf(fp1,"\n");
	}
	fprintf(fp1,"\n");
	fclose(fp1);

	sprintf(out,"%s/output_run_re_infections+slow_matrix_%d_phase_%d.txt",path_def,uncertainty>0,phase);
	if(iter==0)
		fp1=fopen(out,"w");
	else
		fp1=fopen(out,"a");
	for(a1=0;a1<AG;a1++)
	{
		for(a2=0;a2<AG;a2++)
			fprintf(fp1,"%lf ",inf_3[a1][a2]);
		fprintf(fp1,"\n");
	}
	fprintf(fp1,"\n");
	fclose(fp1);
}
