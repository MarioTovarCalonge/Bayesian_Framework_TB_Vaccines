#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "string.h"
#include "time.h"
#include "unistd.h"

//Common functions, macros and variables for runner and fitter
#include "TB_functions.h"

//Models:
//void model_fit(double *c1, double *c2, double *c3, double *c4, double *c5, char **mf, int *writ);
void model_fit(char **mf, int *writ, int *con_mod, double *objeto);

//Loading functions
void load_dynamic_variables(char demo_r[100],double s[states][AG][EG],double S[states],int inicio,int final);
void load_initial_conditions(double s[states][AG][EG],double S[states]);
void load_fit_pyramids(double s[states][AG][EG], char pyr[string_length], int inicio, int final);

//Recording functions
void write_in_cond();
void save_output_micro(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int year, double *vv);
void read_fast(char **nation, int *iter, int *ind, int *len, int *len_row, int *len_col, double *p);

void save_output_micro_file(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int start);
void rewrite_output_object();

//Integrator
int rk4(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double S[states], double f_micro[AG][EG][internal_fluxes], double t);
int deriv(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double ds[states][AG][EG],double t,double S[states], double f_micro[AG][EG][internal_fluxes]);
int calculate_fluxes(int a, int e,double t,double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double f[fluxes],double S[states]);
void ageing(double ds[states][AG][EG], double s[states][AG][EG],double S[states],double vp[variable_parameters][AG][EG],double t);
double semiempiric_derivative_fitted_population(int a,int e,double t);

int write_flag;
double scan;
char master_file_arg[string_length];

void model_fit(char **mf, int *writ, int *con_mod, double *objeto){
	double correction, t, taux, sick = 0;
	int dynamics_type_aux, records, i, a, e, stationarity_reached = 0, err = 0;
	double slope_calcs[2][2];
    char aux_str[string_length];

    cm_mode_all = *con_mod;

    //######################################--LOAD-OF-INPUTS--############################################
    strcpy(master_file_arg, mf[0]);
    master_file_reader(master_file_arg);
    load_window(window_file);
    
    load_dictionary(dictionary_file); //Reading states names
    load_parameters(global_parameters_file, global_parameters, gp);
	load_parameters(regional_parameters_file, regional_parameters, rp);
    load_variable_fitter(variable_parameters, vp);
    load_contact_matrix(contacts_file, cm);
    //load_cm_modified(contacts_file, cm);
	load_dynamic_variables(demography_file, s, S, year_min, year_max);
    
    if (strcmp(name_country, simulation_name) != 0)
        sprintf(aux_str, "../Outputs/%s/Iter/%s/Fitted_inputs/Parameters_to_fit.txt", name_country, simulation_name);
    else
	    sprintf(aux_str, "../Outputs/%s/Fitted_inputs/Parameters_to_fit.txt", simulation_name);

    double aux_h[5][AG];
    load_fitting_parameters(aux_str, gp, rp, aux_h);

    //--------------------------------------------------------------------------------------------------//

    write_flag = *writ;
    int span = year_max - year_min + 1;

    //#############################--Alocación de memoria--####################################
    double macro[observables+annuals];
	double **temporal_measure_macro;//*inc,*mort,*pop,*vac;
    double f_micro[AG][EG][internal_fluxes];

	temporal_measure_macro = (double **)malloc((observables+annuals)*sizeof(double *));
	for(i=0; i<observables+annuals; i++)
		temporal_measure_macro[i] = (double *)malloc((span)*sizeof(double));


    //--------------------------------------------------------------------------------------------------//

	for(a=0; a<AG; a++){
		for(e=0; e<EG; e++){
			for(i=0; i<regional_parameters; i++)
				rp_base[i][a][e] = rp[i][a][e];
	    }
	}
	
    for(a=0; a<AG; a++){
        for(e=0; e<EG; e++){
	        rp[3][a][e] = rp_base[3][a][e]*aux_h[0][a]; //Initial diagnosis rate smear-pos
		    for(i=4; i<6; i++)
			    rp[i][a][e] = rp_base[i][a][e]*aux_h[0][a]*rp[18][a][e];//Diagnosis rate smear-neg. and extra-pulm, times eta
		    rp[6][a][e] = aux_h[1][a];
            rp[7][a][e] = 0;
			rp[8][a][e] = 0; //For termalization, time variation is zero
            rp[2][a][e] = aux_h[4][a];
        }
	}

	dynamics_type_aux = dynamics_type;
	if(dynamics_type==2)
		dynamics_type = 1;

    //######################################--Term--############################################
	t = 0;
	//printf("termalizing\n");
	for(records=0; (stationarity_reached==0); records++){
		start_record(macro, s, S);

		for(taux=0; (taux<record_window && err==0); t+=H){
			err = rk4(gp, vp, rp, cm, s, S, f_micro, t);
			taux += H;
			if(err==1) break;
		}
		if(err==1){
			forced_fitting_flag = 2;
			printf("Problems termalizing\n");
			break;
		}
		stationarity_reached = end_record(macro, s, S, slope_calcs, records);
	}
	//printf("termalized \n");
    
    //---------------------------------------------------------------------------------------------//

	if(err==0){
		for(a=0; a<AG; a++){
			for(e=0; e<EG; e++){
				rp[7][a][e] = aux_h[2][a];//Loading time variation now
				rp[8][a][e] = aux_h[3][a];
			}
		}
		
		if(write_flag==1){
			write_in_cond();
			//return;
		}
        
		for(i=0; i<states; i++)
			S[i] = 0;

		for(a=0; a<AG; a++){
			for(e=0; e<EG; e++)
			{
				s[reservories][a][e] = 0;
				if(rp[2][a][e] < 0)
				{
					correction = 0;
					for(i=1; i<reservories; i++)
					{
						correction -= s[i][a][e]*rp[2][a][e];
						s[i][a][e] = (1+rp[2][a][e])*s[i][a][e];
						s[reservories][a][e] += s[i][a][e];
						S[i] += s[i][a][e];
					}
					s[0][a][e] = s[0][a][e] + correction;
					s[reservories][a][e] += s[0][a][e];
					S[0] += s[0][a][e];
					for(i=reservories+1; i<states; i++)
						s[i][a][e] = 0;
				}
				else
				{
					correction = rp[2][a][e]*s[0][a][e];
					s[0][a][e] = (1-rp[2][a][e])*s[0][a][e];
					s[reservories][a][e] += s[0][a][e];
					S[0] += s[0][a][e];
					sick=0;
					for(i=1; i<reservories; i++)
						sick += s[i][a][e];
					if(sick == 0)
						sick = 0.000001;
					for(i=1; i<reservories; i++)
					{
						s[i][a][e] += correction*s[i][a][e]/sick;
						s[reservories][a][e] += s[i][a][e];
						S[i] += s[i][a][e];
					}
					for(i=reservories+1; i<states; i++)
						s[i][a][e] = 0;
				}
			}
		}
        
		for(i=0;i<reservories;i++)
			S[reservories] += S[i];

		t=0;
		dynamics_type = dynamics_type_aux;

		for(records=0; records<=(year_max-year_min); records++){

            for(i=0; i<internal_fluxes; i++){
		        for(a=0; a<AG; a++){
			        for(e=0; e<EG; e++){
                        f_micro[a][e][i] = 0.0;
			        }   
		        }
	        }

			start_record(macro, s, S);
			for(taux=0; (taux<record_window && err==0); t+=H)
			{
				err = rk4(gp, vp, rp, cm, s, S, f_micro, t);
				taux += H;
				if(err==1) break;
			}
			if(err==1)
			{
				forced_fitting_flag = 2;
				break;
			}
			end_record(macro, s, S, slope_calcs, records);

            //save_output_micro_file(s, f_micro, 0, records);
            //save_output_micro_file(s, f_micro, 1, records);
            //save_output_micro_file(s, f_micro, 2, records);

            save_output_micro(s, f_micro, 0, records, objeto);
            
			for(i=0;i<observables+annuals;i++){
				temporal_measure_macro[i][records] = macro[i];
			}
		}
	}

    //rewrite_output_object();
    
	/*if(err==0){		
		printf("Inc %d: %lf\n", year_min, temporal_measure_macro[1][0]);
        printf("Inc %d: %lf\n", year_min+span-1, temporal_measure_macro[1][span-1]);
		printf("Mort %d: %lf\n", year_min, temporal_measure_macro[3][0]);
		printf("Mort %d: %lf\n\n", year_min+span-1, temporal_measure_macro[3][span-1]);
	}*/
}

int calculate_fluxes(int a, int e,double t,double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double f[fluxes],double S[states])
{
	double infect, pop, lambda;
	int i, a_prime, e_prime;
	int err;
	double beta, d, P_matrix[AG][AG];

	lambda=0;

	//Unscaled Infection force
	if(contact_type==1)
	{
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

		for(a_prime=0; a_prime<AG; a_prime++)
		{
			for(e_prime=0; e_prime<EG; e_prime++)
			{
				infect = P_matrix[a][a_prime]*vp[19][a_prime][e_prime]*
				(
				 s[3][a_prime][e_prime]+
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*s[4][a_prime][e_prime]+
				 vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[9][a_prime][e_prime]+
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[11][a_prime][e_prime]+
				 s[18][a_prime][e_prime]
				);

				if(infect>0)
					lambda += infect/s[reservories][a_prime][e_prime];
			}
		}
	}
	else
	{
		infect = 0;
		pop = 0;
		for(a_prime=0; a_prime<AG; a_prime++)
		{
			for(e_prime=0; e_prime<EG; e_prime++)
			{
				infect+=
				(
				 s[3][a_prime][e_prime] +
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*s[4][a_prime][e_prime] +
				 vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[9][a_prime][e_prime] +
				 vp[16][a_prime][e_prime]*gp[11][a_prime][e_prime]*vp[17][a_prime][e_prime]*gp[12][a_prime][e_prime]*s[11][a_prime][e_prime] +
				 s[18][a_prime][e_prime]
				);

				pop += s[reservories][a_prime][e_prime];
			}
		}
		lambda = 2*infect/pop;
	}

	//Calculate beta
	beta = funcion_param_t(rp[6][a][e], rp[7][a][e], t, 2*rp[6][a][e], 0);
	//Calculate d
	d = funcion_param_t(rp[3][a][e], rp[8][a][e], t, 12.166667, 0);

	/* Calculating Scaled Infection Force and Diagnosis Rates */
	lambda*=beta;
	f[9] = vp[24][a][e]*d*s[3][a][e];
	f[10] = vp[25][a][e]*rp[18][a][e]*d*s[4][a][e];
	f[11] = vp[26][a][e]*rp[18][a][e]*d*s[5][a][e];

	/* Infection */
	f[1] = lambda*vp[6][a][e]*vp[10][a][e]*rp[9][a][e]*s[0][a][e]; //To fast-progression
	f[2] = lambda*vp[6][a][e]*(1-vp[10][a][e]*rp[9][a][e])*s[0][a][e]; //To slow-progression

    //printf("%lf %lf\n", rp[10][a][e], rp[11][a][e]);
	/*Fast-Progression*/
	f[3] = vp[12][a][e]*gp[0][a][e]*rp[10][a][e]*s[1][a][e]; //To smear-positive
	f[4] = vp[12][a][e]*gp[0][a][e]*(1-rp[10][a][e]-rp[11][a][e])*s[1][a][e]; //To smear-negative
	f[5] = vp[12][a][e]*gp[0][a][e]*rp[11][a][e]*s[1][a][e]; //To extra-pulmonary

	/*Slow-Progression*/
	f[6] = vp[13][a][e]*gp[1][a][e]*rp[10][a][e]*s[2][a][e]; //To smear-positive
	f[7] = vp[13][a][e]*gp[1][a][e]*(1-rp[10][a][e]-rp[11][a][e])*s[2][a][e]; //To smear-negative
	f[8] = vp[13][a][e]*gp[1][a][e]*rp[11][a][e]*s[2][a][e]; //To extra-pulmonary

	/*Natural Recovery */
	f[12] = vp[18][a][e]*gp[10][a][e]*s[3][a][e]; //In smear-positives
	f[13] = vp[18][a][e]*gp[10][a][e]*s[4][a][e]; //In smear-negatives
	f[14] = vp[18][a][e]*gp[10][a][e]*s[5][a][e]; //In extra-pulmonary

	/*Treatment*/
	f[15] = gp[6][a][e]*rp[13][a][e]*s[6][a][e]; //smear-positive: default
	f[16] = gp[6][a][e]*rp[16][a][e]*s[7][a][e]; //smear-negative: default
	f[17] = gp[6][a][e]*rp[16][a][e]*s[8][a][e]; //extra-pulmonary: default
	f[18] = gp[6][a][e]*(1-rp[12][a][e]-rp[13][a][e]-rp[14][a][e])*s[6][a][e]; //smear-positive: success
	f[19] = gp[6][a][e]*(1-rp[15][a][e]-rp[16][a][e]-rp[17][a][e])*s[7][a][e]; //smear-negative: success
	f[20] = gp[6][a][e]*(1-rp[15][a][e]-rp[16][a][e]-rp[17][a][e])*s[8][a][e]; //extra-pulmonary: success

	/*Smear-progression*/
	f[21] = vp[19][a][e]*gp[7][a][e]*s[4][a][e]; //Undiagnosed
	f[22] = vp[19][a][e]*gp[7][a][e]*s[7][a][e]; //Diagnosed (in treatment)

	/*Relapses*/
	f[23] = vp[20][a][e]*gp[8][a][e]*s[9][a][e]; //Default, smear-positive
	f[24] = vp[20][a][e]*gp[8][a][e]*s[11][a][e]; //Default, smear-negative
	f[25] = vp[20][a][e]*gp[8][a][e]*s[13][a][e]; //Default, extra-pulmonary
	f[26] = vp[21][a][e]*gp[9][a][e]*s[10][a][e]; //Natural, smear-positive
	f[27] = vp[21][a][e]*gp[9][a][e]*s[12][a][e]; //Natural, smear-negative
	f[28] = vp[21][a][e]*gp[9][a][e]*s[14][a][e]; //Natural, smear-extra-pulmonary

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
		f[0]=0;//births and not TB deaths are not directly calculated, but from the database
		for(i=39; i<56; i++)
			f[i] = 0;
	}
	else
	{
		printf("algorithm not ready for dyn_type>2\n");
		exit(1);
	}

	/*Treatment (Second part)*/
	f[56] = rp[12][a][e]*gp[6][a][e]*s[6][a][e]; //Failure, smear-positive
	f[57] = gp[6][a][e]*rp[14][a][e]*s[6][a][e]; //Death, smear-positive
	f[58] = gp[6][a][e]*rp[17][a][e]*s[7][a][e]; //Death, smear-negative
	f[59] = gp[6][a][e]*rp[17][a][e]*s[8][a][e]; //Death, extra-pulmonary

	f[60] = gp[2][a][e]*s[18][a][e]; //Death of failure

	/* Relapses from success */
	f[61] = vp[22][a][e]*gp[14][a][e]*s[15][a][e]; //Smear-positive
	f[62] = vp[22][a][e]*gp[14][a][e]*s[16][a][e]; //Smear-negative
	f[63] = vp[22][a][e]*gp[14][a][e]*s[17][a][e]; //Extra-pulmonary

	/*Reinfections from success to fast latency*/
	for(i=0; i<3; i++){
		f[64+i] = vp[6][a][e]*vp[10][a][e]*rp[9][a][e]*vp[11][a][e]*gp[5][a][e]*lambda*s[15+i][a][e];
        //f[64+i] = 0.0;
    }

	/*Reinfections from success to slow latency (not considered)*/
	for(i=0; i<3; i++){
        f[67+i] = 0;
        //f[67+i] = vp[6][a][e]*(1 - vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[15+i][a][e];
	}

	/*Treatment (Third part)*/
	f[70] = rp[15][a][e]*gp[6][a][e]*s[7][a][e]; //Failure, smear-negative
	f[71] = rp[15][a][e]*gp[6][a][e]*s[8][a][e]; //Failure, extra-pulmonary

	/*Reinfections followed by slow-progression*/
	/*Do not have an effect in the dynamics, but it is a good observable*/
	for(i=0; i<6; i++)
		f[72+i] = vp[6][a][e]*(1.0-vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[9+i][a][e];	//After treatment
	f[78] = vp[6][a][e]*(1.0-vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[2][a][e]; //Before treatment (from slow to fast latency)
	for(i=0; i<3; i++)
		f[79+i] = vp[6][a][e]*(1.0-vp[10][a][e]*rp[9][a][e])*vp[11][a][e]*gp[5][a][e]*lambda*s[15+i][a][e];

	/*Potential infection*/
	/*Infections that we will have if all the population was susceptible*/
	f[82] = vp[6][a][e]*lambda*s[reservories][a][e];

	int ii;
	for(ii=0; ii<internal_fluxes; ii++)
	{
		if(f[ii]<0)
		{
			//printf("negative flux\n");
			//printf("rp[9][a][e] %lf s[0][a][e] %lf,(a %d e %d t %lf) flux[%d]=%.16lf lambda %.16lf (%.16lf %lf %lf)(%.16lf %lf %lf)\n",rp[9][a][e],s[0][a][e],a,e,t,ii,f[ii],lambda,rp[3][a][e],rp[4][a][e],rp[5][a][e],rp[6][a][e],rp[7][a][e],rp[8][a][e]);
			err = 0;
            //err = 1;
			//break;
		}
		else
			err = 0;
	}
	return err;
}

int deriv(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double ds[states][AG][EG],double t,double S[states], double f_micro[AG][EG][internal_fluxes])
{
	int i, a, e;
	double f[fluxes];
	int err = 0;
	for(a=0; a<AG; a++)
	{
		for(e=0; e<EG; e++)
		{
			/*for(i=0; i<=reservories; i++)
			{
				if(s[i][a][e]<0)
					printf("t=%lf negative value in deriv: s[%d][%d][%d]=%lf\n", t, i, a, e, s[i][a][e]);
			}*/
			err = calculate_fluxes(a, e, t, gp, vp, rp, cm, s, f, S);
			calculate_derivatives(a, e, ds, f);
            for(i=0; i<internal_fluxes; i++){
                f_micro[a][e][i] = f[i];
            }
			if(err == 1)
				return(err);
		}

		if(err == 1)
			return(err);
	}
	ageing(ds, s, S, vp, t);
	return err;
}

int rk4(double gp[global_parameters][AG][EG],double vp[variable_parameters][AG][EG],double rp[regional_parameters][AG][EG],double cm[AG][AG][EG][EG], double s[states][AG][EG],double S[states], double f_micro[AG][EG][internal_fluxes], double t)
{
	int i,a,e;
	double ds[states][AG][EG], dsm[states][AG][EG], dst[states][AG][EG], st[states][AG][EG], St[states];
    double f_m1[AG][EG][internal_fluxes], f_m2[AG][EG][internal_fluxes], f_m3[AG][EG][internal_fluxes], fs[AG][EG][internal_fluxes];
	double hh, h6;
	int err = 0;
	hh = H*0.5;
	h6 = H/6.0;

	err = deriv(gp, vp, rp, cm, s, ds, t, S, fs);                /*ds=f(x0)*/
	if(err==1)
	{
		printf("Runge-Kutta fails in first call to deriv\n");
		return err;
	}

	for(i=0; i<states; i++)
	{
		St[i] = 0;
		for(a=0; a<AG; a++)
		{
			for(e=0; e<EG; e++)
			{
				st[i][a][e] = s[i][a][e] + hh*ds[i][a][e];
				St[i] += st[i][a][e];
			}
		}
	}
	err = deriv(gp, vp, rp, cm, st, dst, t+hh, St, f_m1);           /*st es x1, dst=f(x1)*/
	if(err==1)
		return err;

	for(i=0; i<states; i++)
	{
		St[i] = 0;
		for(a=0; a<AG; a++)
		{
			for(e=0; e<EG; e++)
			{
				st[i][a][e] = s[i][a][e]+hh*dst[i][a][e];
				St[i] += st[i][a][e];
			}
		}
	}
	err = deriv(gp, vp, rp, cm, st, dsm, t+hh, St, f_m2);           /*xt es x2, dxm=f(x2)*/
	if(err==1)
		return err;

	for(i=0; i<states; i++)
	{
		St[i] = 0;
		for(a=0; a<AG; a++)
		{
			for(e=0; e<EG; e++)
			{
				st[i][a][e] = s[i][a][e] + H*dsm[i][a][e];  /*st es x3, dsm=f(x1)+f(x2)*/
				St[i] += st[i][a][e];
				dsm[i][a][e] += dst[i][a][e];
			}
		}
	}
	err = deriv(gp, vp, rp, cm, st, dst, t+H, St, f_m3);           /*dst=f(x3)*/
	if(err==1)
		return err;

	for(i=0; i<states; i++)
	{
		S[i] = 0;
		for(a=0; a<AG; a++)
		{
			for(e=0; e<EG; e++)
			{
				s[i][a][e] += h6*(ds[i][a][e]+dst[i][a][e] + 2.0*dsm[i][a][e]);  /*st es x3, dsm=f(x1)+f(x2)*/
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

void load_initial_conditions(double s[states][AG][EG],double S[states])
{
	FILE *fich;
	int entries;
	int i, a, e;
	int amin, amax, emin, emax, index;
	double value, vmin, vmax;
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
		//Reason is that we might want to run intervals considering, for example, demographic uncertainty,
		//but not uncertainty from initial conditions
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
					s[index][a][e]=0;//If for any reason you pass me obsevables in the file (you shouldnt), I read them, but immediately put them at zero
			}
		}
	}

	double auxiliar;
	for(i=0;i<states;i++)//I calculate S[reservories] (age-aggregated population) here because it is not done when I load demography
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
					s[i][a][e] *= s[reservories][a][e];
					auxiliar += s[i][a][e];
					S[i] += s[i][a][e];
				}
				else
				{
					if(i>reservories)
						s[i][a][e]=0;
					S[i] += s[i][a][e];
				}
			}
			s[0][a][e] = s[reservories][a][e] - auxiliar;
			if(s[0][a][e]<0)
			{
				printf("error loading initial conditions: negative class (a=%d e=%d) s[0][a][e]=%f\n",a,e,s[0][a][e]);
                s[0][a][e] = -s[0][a][e];
				//exit(1);
			}
			S[0] += s[0][a][e];
		}
	}
	//write_in_cond();
	fclose(fich);
}

void ageing(double ds[states][AG][EG], double s[states][AG][EG],double S[states],double vp[variable_parameters][AG][EG],double t)
{
	int i,a,e;
	double sick_density,numerator,denominator;
	for(e=0;e<EG;e++)
	{
		for(i=0;i<=reservories;i++)
		{
			ds[i][0][e]-=s[i][0][e]/(double)age_window;
			for(a=1;a<AG;a++)
				ds[i][a][e] += (s[i][a-1][e]-s[i][a][e])/(double)age_window;
		}
	}
	double delta;

	if (dynamics_type==0)
	{
		numerator=0;
		denominator=0;
		for(a=start_maternity;a<=stop_maternity;a++)
		{
			for(e=0;e<EG;e++)
			{
				numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[26][a][e];
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
					delta-=ds[i][a][e];//Total increment in the population (changing sign)
			}
		}

		for(e=0;e<EG;e++)
		{
			ds[0][0][e]+=(vp[0][0][e])*delta*(1-sick_density);
			ds[1][0][e]+=(vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
			ds[2][0][e]+=(vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
			if(mother_child_infections_accounted==1)
				ds[reservories+5][0][e]+=(vp[0][0][e])*delta*(sick_density);
			ds[reservories][0][e]+=delta*vp[0][0][e]; //Thus I cancel total increment, introducing newborns
		}
		cont_clone++;
		ds[reservories+4][0][0]+=(vp[0][0][1])*delta;//When summing vaccinated individual, we have to add for newborn vaccination
	}
	else
	{
		if (dynamics_type==1)
		{
			numerator=0;
			denominator=0;
			for(a=start_maternity;a<=stop_maternity;a++)
			{
				for(e=0;e<EG;e++)
				{
					numerator+=(s[3][a][e]+s[4][a][e]+s[5][a][e]+s[9][a][e]+s[11][a][e]+s[13][a][e]+s[18][a][e])*gp[13][a][e]*vp[26][a][e];
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
					delta-=ds[i][0][e];//Population lost in agegroup 0
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
			ds[reservories+4][0][0]+=(vp[0][0][1])*delta;//Correcting vaccination with newborns vacs
			for(a=1;a<AG;a++)
			{
				delta=0;
				denominator=0;

				for(e=0;e<EG;e++)
				{
					for(i=0;i<reservories;i++)
					{
						delta-=ds[i][a][e];//Population lost in agegroup a
						denominator+=s[i][a][e];//Population aggregated for agegroup a
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
			if (dynamics_type==2)
			{
				numerator = 0;
				denominator = 0;
				for(a=start_maternity; a<=stop_maternity; a++)
				{
					for(e=0; e<EG; e++)
					{
						numerator += (s[3][a][e] + s[4][a][e] + s[5][a][e] + s[9][a][e] + s[11][a][e] + s[13][a][e] + s[18][a][e])*gp[13][a][e]*vp[26][a][e];
						for(i=0; i<reservories; i++)
							denominator += s[i][a][e];
					}
				}

				if(denominator != 0)
					sick_density = numerator/denominator;
				else
					sick_density = 0;
				delta = 0;
				for(e=0; e<EG; e++)
				{
					for(i=0; i<reservories; i++)
						delta -= ds[i][0][e];//Population lost agegroup 0
				}
				delta += semiempiric_derivative_fitted_population(0,0,t);//Only once, this pyramid is the sum for all e-groups
				for(e=0; e<EG; e++)
				{
					ds[0][0][e] += (vp[0][0][e])*delta*(1-sick_density);
					ds[1][0][e] += (vp[0][0][e])*delta*(sick_density)*rp[9][0][e]*vp[10][a][e];
					ds[2][0][e] += (vp[0][0][e])*delta*(sick_density)*(1-rp[9][0][e]*vp[10][a][e]);
					ds[reservories][0][e] += delta*vp[0][0][e];
					if(mother_child_infections_accounted == 1)
						ds[reservories+5][0][e] += (vp[0][0][e])*delta*(sick_density);
				}
				ds[reservories+4][0][0] += (vp[0][0][1])*delta;
				for(a=1;a<AG;a++)
				{
					delta=0;
					denominator=0;

					for(e=0;e<EG;e++)
					{
						for(i=0;i<reservories;i++)
						{
							delta -= ds[i][a][e];//Population lost in agegroup a
							denominator += s[i][a][e];//Population aggregated for agegroup a
						}
					}
					delta += semiempiric_derivative_fitted_population(a,0,t);
					for(e=0; e<EG; e++)
					{
						for(i=0; i<=reservories; i++)
							ds[i][a][e] += delta*(s[i][a][e]/denominator);
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

double semiempiric_derivative_fitted_population(int a, int e, double t)
{
	int i;
	double result;
	result=0;
	for(i=1; i<=fitting_degree; i++)
		result += pyramid_coefficient[a][e][i]*i*pow(t,(i-1));
	return(result);
}

void load_fit_pyramids(double s[states][AG][EG], char pyr[string_length], int start, int end)
{
	FILE *fich;
	double pyramid[AG][EG][demo_span_max];
	int entries, span;
	int i, a, e, t;
	int amin, amax, emin, emax, tmin, tmax;
	double value, vmin, vmax;
	span = end - start + 1;

	if(span>demo_span_max)
	{
		printf("too many years when loading pyramids\n");
		exit(1);
	}

    //#############-Lectura-pirámides-##################
	fich = fopen(pyr, "rt");
	fflush(stdout);

	scan=fscanf(fich,"entries=%d \n",&entries);
	scan=fscanf(fich,"amin \t amax \t emin \t emax \t tmin \t tmax \t value \t interval\n");

	for(i=0;i<entries;i++)
	{
		scan = fscanf(fich,"%d \t %d \t %d \t %d \t %d \t %d \t %lf \t %lf \t %lf\n", &amin, &amax, &emin, &emax, &tmin, &tmax, &value, &vmin, &vmax);

		amin = (int)(amin/age_window);
		amax = (int)(amax/age_window);

		if(amax>AG && amin<=AG) amax=AG; //Only case that could fail

		if(tmin<year_min) continue;

		if(dynamics_type!=2)
		{
			//If dynamic type is not 2, only the inicial pyramid is read and considered constant
			if(tmin==start)
				tmax = end+1;
			else
				continue;
		}

		for (a=amin; a<amax; a++)
		{
			for(e=emin; e<emax; e++)
			{
				for(t=tmin; t<tmax; t++)
					pyramid[a][e][t-year_min] = value;
			}
		}
	}
	fclose(fich);
    //###################-FIN-##########################

    FILE *fp_local;
	char file_name[string_length];
	int j, j_aux, k_aux;

    
	sprintf(file_name, "%s/Pyr_coef_fit.txt", path_fitted_inputs);
	fp_local=fopen(file_name, "rt");
	for(j=0; j<AG; j++)
	{
	    //for(k=0; k<EG; k++)
	    //{
		    scan = fscanf(fp_local, "%d %d ", &j_aux, &k_aux);
		    for(i=0; i<fitting_degree+1; i++){
			    scan = fscanf(fp_local, "%lf ", &pyramid_coefficient[j][0][i]);
                pyramid_coefficient[j][1][i] = 0.0;
                pyramid_coefficient[j][2][i] = 0.0;
            }
		    scan = fscanf(fp_local, "\n");
	    //}
	   }
	fclose(fp_local);

	for(a=0; a<AG; a++)
	    s[reservories][a][0] = pyramid[a][0][start-year_min];//Initial Populations

	for (e=1; e<EG; e++){
	    for(a=0; a<AG; a++)
		    s[reservories][a][e] = 0;//pyramid[a][e][start-year_min];
	}
    
}

void load_dynamic_variables(char demo_r[100],double s[states][AG][EG],double S[states],int start,int end)
{
	load_fit_pyramids(s, demo_r, start, end);
	load_initial_conditions(s, S);
}

void write_in_cond()
{
	FILE *fp_local;
	int a, i;
	double frac[reservories][AG], total;
	char file_name[string_length];


	sprintf(file_name,"%s/Initial_Conditions/initial_conditions.txt",path_fitted_inputs);
	fp_local = fopen(file_name,"wt");

	for(a=0; a<AG; a++)
	{
		total = 0;
		for(i=0; i<reservories; i++)
			total += s[i][a][0];
		for(i=0; i<reservories; i++)
			frac[i][a] = s[i][a][0]/total;
	}

	fprintf(fp_local, "entries = %d \n", reservories*(AG+1));
	fprintf(fp_local, "description \t index \t amin \t amax \t emin \t emax \t value \t interval \n");
	for(i=0; i<reservories; i++)
	{
		for(a=0; a<AG; a++)
			fprintf(fp_local, "%s \t %d \t %d \t %d \t %d \t %d \t %.12lf \t %.12lf \t %.12lf \n", states_dictionary[i],i,a*5,(a+1)*5,0,1,frac[i][a],frac[i][a],frac[i][a]);
		fprintf(fp_local, "%s \t %d \t %d \t %d \t %d \t %d \t %.12lf \t %.12lf \t %.12lf \n", states_dictionary[i],i,0,120,1,3,0.0,0.0,0.0);
	}
	fclose(fp_local);
}

void save_output_micro(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int year, double *vv){
	int i, a, shift, index;
    int max_fit_years = year_max - year_min + 1;
    int max_col_per_year = AG*reservories + 1 + internal_fluxes*AG;

    //shift = year*max_col_per_year;
    shift = year*max_col_per_year + e*max_col_per_year*max_fit_years;
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

void read_fast(char **nation, int *iter, int *ind, int *len, int *len_row, int *len_col, double *p){
    FILE *fp1;

	char com_aux[1000];
    int i, j, k, l, scan, ag=15, year, ind_aux;
    double val, fl[*len_row], prev[*len_col][ag];

    for(i=0; i<(*len_col); i++){
        for(j=0; j<ag; j++){
            prev[i][j] = 0.0;
        }
    }

    if (*iter==0){
        sprintf(com_aux, "../Outputs/%s/Fitter_object.txt", nation[0]);
    }
    else{
	    sprintf(com_aux, "../Outputs/%s/Iter/%s_%03d/Fitter_object.txt", nation[0], nation[0], *iter);
    }

    fp1 = fopen(com_aux, "rt");

    for(i=0; i<*len_col; i++){
        scan=fscanf(fp1, "%d ", &year);
        for(j=1; j<*len_row; j++){
            scan=fscanf(fp1, "%lf ", &val);
            fl[j] = val;
        }
        scan=fscanf(fp1, "\n");

        for(k=0; k<*len; k++){
            for(l=0; l<(ag); l++){
                ind_aux = ag*ind[k] + 1 + l;
                prev[i][l] += fl[ind_aux] ;
            }     
        }
    }

    for(k=0; k<(*len_col); k++){
        for(l=0; l<ag; l++){
            p[k + l*(*len_col)] = prev[k][l];
        }
    }
    fclose(fp1);
}

void save_output_micro_file(double s[states][AG][EG], double f[AG][EG][internal_fluxes], int e, int start){
    FILE *fp1;
	char com_aux[1000];
    if (strcmp(name_country, simulation_name) != 0)
        sprintf(com_aux, "../Outputs/%s/Iter/%s/", name_country, simulation_name);
    else
	    sprintf(com_aux, "../Outputs/%s/", simulation_name);

    sprintf(com_aux, "%s/Fitter_object_%d.txt", com_aux, e);

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
    sprintf(com, "%s/Fitter_object_%d.txt", com_aux, E);
	fp1 = fopen(com, "a");

    for(E=1; E<EG; E++){
        sprintf(com, "%s/Fitter_object_%d.txt", com_aux, E);
	    fp2 = fopen(com, "r");
        for(i = 0; i<(year_max - year_min + 1); i++){
            RetChar = fgets(line,30000,fp2);
	        fprintf(fp1,"%s",line);
        }
        fclose(fp2);
        sprintf(sys_com, "rm %s", com);
        systemRet = system(sys_com);
    }
    fclose(fp1);
    sprintf(sys_com, "cp %s/Fitter_object_%d.txt %s/Fitter_object.txt", com_aux, 0, com_aux);
    systemRet = system(sys_com);
    sprintf(sys_com, "rm %s/Fitter_object_%d.txt", com_aux, 0);
    systemRet = system(sys_com);
}


