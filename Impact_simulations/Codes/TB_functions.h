//COMMON MACROS
#define string_length_short 100
#define string_length 500 
#define annuals 7 //Variables that we record annually (Susceptibles, Infectious, Disease Prevalence, Smear-positive disease prevalence and Pulmonary disease prevalence (both smear-pos and neg)), Infection Force, Susceptibles_weighted
#define observables 9 //Varibles that we follow continuosly: Population, Incidence over T states, Incidence over D states (real incidence), Mortality, Vaccination, Primary Infections, re-infections+fast-prog., re-infections+slow-prog., potential_infections
#define reservories 19 //Microscopic states of the disease (Compartments)
#define states observables+reservories //Sum of both
#define AG 15 //Age Groups
#define age_window 5 //Duration of each age group (years)
#define EG 3  //Vaccination branches
#define global_parameters 15 //Disease paremeters considered the same in every region
#define regional_parameters 19 //Disease parameters that vary among countries
#define variable_parameters 29  //Define the effect of vaccination
#define fitting_degree 9 //Degree of the polynomial p(t) to which we fit the evolution of each age segment of the demographic pyramid. WARNING: Cant be more than the years to fit
#define WHO 2 //0: Treat incidence and mortality separatedly, 1: Group Inc and Mort, 2: Group Inc, Mort and p_treatment (default)
#define max_records 200 //Maximum number of records in model run (*time window has to be enough to reach stationary)
#define demo_span_max 200 //Maximum temporal span covered by demographic pyramids
#define H 0.003 //Integration time step (in years, so it is about 9 hours)
#define H1 1000 // 1/H (so we can spare some divisions in the code)
#define record_window 1 //time between records (in years)
#define internal_fluxes 83 //fluxes inside the same vaccination branch
#define external_fluxes 22 //fluxes that occur between vac branches
#define fluxes external_fluxes+internal_fluxes //total fluxes
#define scaling 1000000.0 //Rates are calculated every "scaling" individuals
#define incidence_to_fit 0 //=0:fit active cases. =1: fit diagnosed cases
#define threshold_slope 0.01 //Precision in the determination of the steady state
#define min_reg_fit 50 //Minimum number of iterations in the steady state before accepting it
#define start_maternity 3 //Age group in which maternity begins (15 years)
#define stop_maternity 7 //Age group in which maternity ends (35 years)
#define maxspan 120 //Maximum number of records (years of simulation)

//FILES AND PATHS
char simulation_name[string_length];
char name_country[string_length];
char master_file[string_length];
char dictionary_file[string_length];
char global_parameters_file[string_length];
char regional_parameters_file[string_length];
char vaccine_parameters_file[string_length];
char contacts_file[string_length];
char demography_file[string_length];
char path_fitted_inputs[string_length];
char initial_conditions_file[string_length];
char fitting_parameters_file[string_length];
char window_file[string_length];

//COMMON VARIABLES
char states_dictionary[states][string_length_short];
double gp[global_parameters][AG][EG], vp[variable_parameters][AG][EG], rp[regional_parameters][AG][EG]; //Disease parameters
double cm[AG][AG][EG][EG];  //Contact matrix
int year_min, year_max, year_end, year_vac, time_stabilization;
int dynamics_type, contact_type;
int parameters_to_fit;
int start_global, cont_clone, iter, forced_fitting_flag;
int flag_demo, uncertainty;
int mother_child_infections_accounted=1;
double pyramid_coefficient[AG][EG][fitting_degree+1],rp_base[regional_parameters][AG][EG];
double s_resource[states][AG][EG], s[states][AG][EG], S[states];
double temporal_measure_macro[observables+annuals][maxspan];
double t_auxiliar;
int cm_mode_all;

//FUNCTIONS
void master_file_reader(char *master_file);
void load_dictionary(char file[string_length]);
void load_parameters(char file[string_length],int parameters, double p[parameters][AG][EG]);
void load_variable_fitter(int parameters, double p[parameters][AG][EG]);
void load_vaccine_modifiers(char file[string_length], int parameters, double varp[parameters][AG][EG]);
void load_vaccine_modifiers_H(char file[string_length], int parameters, double varp[parameters][AG][EG]);
void modify(int parameters, char S_L_R[string_length], double p[parameters][AG][EG], int age, double cov, double blocking, double waning, char tran[string_length], double eff);
void sort(int *v, int n, int *index);

void load_contact_matrix(char file[string_length], double k[AG][AG][EG][EG]);
void load_fitting_parameters(char fn[string_length], double gp[global_parameters][AG][EG], double rp[regional_parameters][AG][EG], double ax[2][AG]);
void load_fitting_parameters_from_R(int *ind, int *cosa, int *a_min, int *a_max, double *valu, double *lb, double *ub, int entries, double gp[global_parameters][AG][EG], double rp[regional_parameters][AG][EG], double ax[5][AG]);
void load_window(char window[string_length]);
double funcion_param_t(double l0, double r0, double t, double lsup, double linf);
void start_record(double macro[observables+annuals], double s[states][AG][EG], double S[states]);
int end_record(double macro[observables+annuals], double s[states][AG][EG],double S[states],double slope_calcs[2][2],int reg);
void calculate_derivatives(int a, int e, double ds[states][AG][EG],double f[fluxes]);
void correct_contact_matrix(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]);
void load_cm_modified(char file[string_length], double k[AG][AG][EG][EG]);
void cm_M0(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]);
void cm_M1(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]);
void cm_M2(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]);
void cm_M3(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]);

void master_file_reader(char *master_file){	
	FILE *fp;
	int scan;
	char command[string_length];
	
	fp=fopen(master_file, "rt");
	
	scan=fscanf(fp, "SIMULATION_NAME %s\n", simulation_name);
    scan=fscanf(fp, "COUNTRY %s\n", name_country);
	
	scan=fscanf(fp, "EXECUTION_MODES\n");
	scan=fscanf(fp, "dynamic_type %d\n", &dynamics_type);
	scan=fscanf(fp, "contact_type %d\n", &contact_type);
	parameters_to_fit = 5;
	
	scan=fscanf(fp, "INPUTS: \n");
	scan=fscanf(fp, "states_dictionary %s\n", dictionary_file);
	scan=fscanf(fp, "global_parameters %s\n", global_parameters_file);
	scan=fscanf(fp, "regional_parameters %s\n", regional_parameters_file);
	scan=fscanf(fp, "contacts %s \n", contacts_file);
	scan=fscanf(fp, "demography %s\n", demography_file);
	scan=fscanf(fp, "window %s\n", window_file);
	
    if (strcmp(name_country, simulation_name) != 0)
        sprintf(path_fitted_inputs, "../Outputs/%s/Iter/%s/Fitted_inputs", name_country, simulation_name);
    else
        sprintf(path_fitted_inputs, "../Outputs/%s/Fitted_inputs", simulation_name);	
	
	fclose(fp);
}

void load_dictionary(char file[string_length]){
	FILE *fp;
	int i,scan;
	fp=fopen(file,"rt");
	for(i=0;i<states;i++)
		scan=fscanf(fp,"%s\n",states_dictionary[i]);
	fclose(fp);
}

void load_parameters(char file[string_length],int parameters,double p[parameters][AG][EG]){
	FILE *fich;
	int entries;
	int i,a,e,scan;
	int amin,amax,emin,emax,index;
	double value,vmin,vmax;
	char name[100];
	
	fich=fopen(file,"rt");
	
	scan=fscanf(fich,"entries = %d \n",&entries);
	scan=fscanf(fich,"description \t index \t amin \t amax \t emin \t emax \t value \t interval\n");
	for(i=0;i<entries;i++){
		scan=fscanf(fich,"%s \t %d \t %d \t %d \t %d \t %d \t %lf \t %lf \t %lf \n",name,&index,&amin,&amax,&emin,&emax,&value,&vmin,&vmax);
		amin=(int)(amin/age_window);
		amax=(int)(amax/age_window);
		if(amax>AG&&amin<=AG) amax=AG; //Único caso factible en el cual las cosas pueden salirse de orden
		if(amin<0||amin>AG||amax<0||amax>AG||emin<0||emin>EG||emax<0||emax>EG||amin>amax||emin>emax) 
		{
			printf("Range error in the reading of the parameters\n");
			exit(1);
		}
		for (a=amin;a<amax;a++){
			for(e=emin;e<(emax);e++)
				p[index][a][e]=value;			
		}		
	}

	fclose(fich);
}

void load_variable_fitter(int parameters, double p[parameters][AG][EG])
{
	int i, a, e;

    for (a=0;a<AG;a++){
        p[0][a][0] = 1.0;
        p[0][a][1] = 0.0;
        p[0][a][2] = 0.0;
    }		

    for(i=1; i<3; i++){
		for(a=0; a<AG; a++){
            for(e=0; e<EG; e++)
                p[i][a][e] = 0.0;
        }
	}
	
	for(i=3; i<5; i++){
		for(a=0; a<AG; a++){
            for(e=0; e<EG; e++)
                p[i][a][e] = 1.0;
        }
	}

    for(a=0; a<AG; a++){
        for(e=0; e<EG; e++)
            p[5][a][e] = 0.0;
    }
	
	for(i=6; i<27; i++){
		for(a=0; a<AG; a++){
            for(e=0; e<EG; e++)
                p[i][a][e] = 1.0;
        }
	}

    for(i=27; i<29; i++){
		for(a=0; a<AG; a++){
            for(e=0; e<EG; e++)
                p[i][a][e] = 0.0;
        }
	}
	
}

void load_vaccine_modifiers(char file[string_length], int parameters, double varp[parameters][AG][EG])
{
    FILE *fich;
	double wan, waning_low, waning_hi, block, blocking_low, blocking_hi, coverage, cov_low, cov_hi;
	int scan, age;
    double efficacy, eff_low, eff_hi;
    char tt[string_length], slr[string_length];

    int max_num = 0;
	fich=fopen(file,"rt");	
    scan=fscanf(fich, "transition age S_L_R efficacy e_low e_hi coverage c_low c_hi waning w_low w_hi blocking b_low b_hi\n");
    while(fscanf(fich, "%s %d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", tt, &age, slr, &efficacy, &eff_low, &eff_hi, &coverage, &cov_low, &cov_hi, &wan, &waning_low, &waning_hi, &block, &blocking_low, &blocking_hi)==15){
        max_num++;
    }
    fclose(fich);
    printf("%d\n", max_num);

    int sent = 0;
    int ages[max_num], index[max_num];
    double things[max_num][12];
    char transit[max_num][string_length];
    char status_target[max_num][string_length];
    int i, j;
    for(i = 0; i< max_num; i++){
        index[i] = i;
        ages[i] = 0;
        for(j = 0; j<12; j++)
            things[i][j] = 0.0;
    }

    fich=fopen(file,"rt");	
    scan=fscanf(fich, "transition age S_L_R efficacy e_low e_hi coverage c_low c_hi waning w_low w_hi blocking b_low b_hi\n");

    while(fscanf(fich, "%s %d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", transit[sent], &ages[sent], status_target[sent], &things[sent][0], &things[sent][1], &things[sent][2], &things[sent][3], &things[sent][4], &things[sent][5], &things[sent][6], &things[sent][7], &things[sent][8], &things[sent][9], &things[sent][10], &things[sent][11])==15){

        sent++;
    }

    printf("%d\n", sent);
    sort(ages, max_num, index);

    for(i=0; i< max_num; i++){
        modify(parameters, status_target[index[i]], varp, ages[i], things[index[i]][3], things[index[i]][9], things[index[i]][6], transit[index[i]], things[index[i]][0]);
    }
    fclose(fich);		
}

int target=0;
double waning_h=0;

void load_vaccine_modifiers_H(char file[string_length], int parameters, double varp[parameters][AG][EG])
{
    FILE *fich;
	double wan, waning_low, waning_hi, block, blocking_low, blocking_hi, coverage, cov_low, cov_hi, trash;
	int scan, age, inttrash;
    double efficacy, eff_low, eff_hi;
    char tt[string_length], slr[string_length];

    int max_num = 0;
	fich=fopen(file,"rt");	
    scan=fscanf(fich,"target = %d \n",&inttrash);
    scan=fscanf(fich,"waning_h = %lf \n",&trash);
    scan=fscanf(fich, "transition age S_L_R efficacy e_low e_hi coverage c_low c_hi waning w_low w_hi blocking b_low b_hi\n");
    while(fscanf(fich, "%s %d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", tt, &age, slr, &efficacy, &eff_low, &eff_hi, &coverage, &cov_low, &cov_hi, &wan, &waning_low, &waning_hi, &block, &blocking_low, &blocking_hi)==15){
        max_num++;
    }
    fclose(fich);
    printf("%d\n", max_num);

    int sent = 0;
    int ages[max_num], index[max_num];
    double things[max_num][12];
    char transit[max_num][string_length];
    char status_target[max_num][string_length];
    int i, j;
    for(i = 0; i< max_num; i++){
        index[i] = i;
        ages[i] = 0;
        for(j = 0; j<12; j++)
            things[i][j] = 0.0;
    }

    fich=fopen(file,"rt");	
    scan=fscanf(fich,"target = %d \n",&target);
    scan=fscanf(fich,"waning_h = %lf \n",&waning_h);
    scan=fscanf(fich, "transition age S_L_R efficacy e_low e_hi coverage c_low c_hi waning w_low w_hi blocking b_low b_hi\n");

    while(fscanf(fich, "%s %d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", transit[sent], &ages[sent], status_target[sent], &things[sent][0], &things[sent][1], &things[sent][2], &things[sent][3], &things[sent][4], &things[sent][5], &things[sent][6], &things[sent][7], &things[sent][8], &things[sent][9], &things[sent][10], &things[sent][11])==15){

        sent++;
    }

    printf("%d\n", sent);
    sort(ages, max_num, index);
    
    for(i=0; i< max_num; i++){
        modify(parameters, status_target[index[i]], varp, ages[i], things[index[i]][3], things[index[i]][9], things[index[i]][6], transit[index[i]], things[index[i]][0]);
    }
  
    fclose(fich);		
}

int eps_AON = 0;
int age_AON = 0;
double cov_AON = 0;

void modify(int parameters, char S_L_R[string_length], double p[parameters][AG][EG], int age, double cov, double blocking, double waning, char tran[string_length], double eff)
{
    int a;
    double shift = 0.0;
    //Index 0: Newborn coverage gives the fraction of newborns added to each branch f_vn
	if(age==0){		
        for (a=0;a<AG;a++){
	        p[0][a][0] = 1.0 - cov;
            p[0][a][1] = cov*(1-blocking);
            p[0][a][2] = cov*blocking;
        }
	}
	
	//Index 1-2-28: Coverage. Applies to the transfer from e=0 (Non-vac branch) to e=1 (Vac branch) and e=2 (Blocked branch)
	if(age!=0){		
        if(strcmp(S_L_R, "S")==0){
            p[1][age][0] = cov*(1-blocking);
            p[1][age][1] = cov*blocking;
            p[1][age][2] = 0;
        }
        else if(strcmp(S_L_R, "L")==0){
            p[2][age][0] = cov*(1-blocking);
            p[2][age][1] = cov*blocking;
            p[2][age][2] = 0;
        }
        else if(strcmp(S_L_R, "R")==0){
            p[28][age][0] = cov*(1-blocking);
            p[28][age][1] = cov*blocking;
            p[28][age][2] = 0;
        }
        else if(strcmp(S_L_R, "SL")==0 || strcmp(S_L_R, "LS")==0){
            p[1][age][0] = cov*(1-blocking);
            p[1][age][1] = cov*blocking;
            p[1][age][2] = 0;
            p[2][age][0] = cov*(1-blocking);
            p[2][age][1] = cov*blocking;
            p[2][age][2] = 0;
        }
        else if(strcmp(S_L_R, "SR")==0 || strcmp(S_L_R, "RS")==0){
            p[1][age][0] = cov*(1-blocking);
            p[1][age][1] = cov*blocking;
            p[1][age][2] = 0;
            p[28][age][0] = cov*(1-blocking);
            p[28][age][1] = cov*blocking;
            p[28][age][2] = 0;
        }
        else if(strcmp(S_L_R, "LR")==0 || strcmp(S_L_R, "RL")==0){
            p[2][age][0] = cov*(1-blocking);
            p[2][age][1] = cov*blocking;
            p[2][age][2] = 0;
            p[28][age][0] = cov*(1-blocking);
            p[28][age][1] = cov*blocking;
            p[28][age][2] = 0;
        }
        else if(strcmp(S_L_R, "SLR")==0 || strcmp(S_L_R, "LRS")==0 || strcmp(S_L_R, "RSL")==0 || strcmp(S_L_R, "LSR")==0 || strcmp(S_L_R, "RLS")==0){
            p[1][age][0] = cov*(1-blocking);
            p[1][age][1] = cov*blocking;
            p[1][age][2] = 0;
            p[2][age][0] = cov*(1-blocking);
            p[2][age][1] = cov*blocking;
            p[2][age][2] = 0;
            p[28][age][0] = cov*(1-blocking);
            p[28][age][1] = cov*blocking;
            p[28][age][2] = 0;
        }
        else{
            printf("No vaccinaton target selected. Rerun with S, L, R, SL, SR, LR or SLR\n");
        }
	}

	//Index 6: Protection against infection, E_beta
    if(strcmp(tran, "E_beta")==0){
        for(a=age; a<AG; a++){
			p[6][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
		}
    }

    //Index 7 Protection against Mortality by pulmonary smear positive TB, E_mu_i+ 
    if(strcmp(tran, "E_mu_i+")==0){ 
        for(a=age; a<AG; a++)
			p[7][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 8 Protection against Mortality by pulmonary smear negative TB, E_mu_i-
    if(strcmp(tran, "E_mu_i-")==0){ 
        for(a=age; a<AG; a++)
			p[8][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 9 Protection against Mortality by non pulmonary TB, E_mu_np
    if(strcmp(tran, "E_mu_np")==0){ 
        for(a=age; a<AG; a++)
			p[9][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }
	
	//Index 10 Protection against fast-progression probability, E_p
	if(strcmp(tran, "E_p")==0){
        for(a=age; a<AG; a++){
			p[10][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
            //printf("%lf \n", p[10][a][1]);
		}
    }

	//Index 11 Protection against reinfection, E_q
    if(strcmp(tran, "E_q")==0){
        for(a=age; a<AG; a++){
			p[11][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
		}
		if(strcmp(S_L_R, "SLR")==0 || strcmp(S_L_R, "LRS")==0 || strcmp(S_L_R, "RSL")==0 || strcmp(S_L_R, "LSR")==0 || strcmp(S_L_R, "RLS")==0){
            for(a=age; a<AG; a++){
                p[10][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
                //printf("%lf \n", p[10][a][1]);
            }
        }
    }
	
	//Index 12 Protection against fast progression rate, E_r
    if(strcmp(tran, "E_r")==0){
        for(a=age; a<AG; a++){
			p[12][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
		}
    }
	
    //Index 13 Protection against slow progression rate, E_rl
    if(strcmp(tran, "E_rl")==0){ 
        for(a=age; a<AG; a++){
			p[13][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
		}
    }

    //Index 16 Protection against Infectiousness reduction coefficient of Dp− with respect to Dp+, E_phi_i-
    if(strcmp(tran, "E_phi_i-")==0){ 
        for(a=age; a<AG; a++)
			p[16][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }
     
    //Index 17 Protection against Infectiousness reduction coefficient of Rp+D with respect to Dp+, E_phi_f
    if(strcmp(tran, "E_phi_f")==0){ 
        for(a=age; a<AG; a++)
			p[17][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 18 Protection against Natural recovery rate, E_nu
    if(strcmp(tran, "E_nu")==0){ 
        for(a=age; a<AG; a++)
			p[18][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 19 Protection against Smear progression rate, E_theta
    if(strcmp(tran, "E_theta")==0){ 
        for(a=age; a<AG; a++)
			p[19][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 20 Protection against Relapse rate for individuals who defaulted treatment, E_rdef
    if(strcmp(tran, "E_rdef")==0){ 
        for(a=age; a<AG; a++)
			p[20][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 21 Protection against Relapse rate for naturally recovered individuals, E_rn
    if(strcmp(tran, "E_rn")==0){ 
        for(a=age; a<AG; a++)
			p[21][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 22 Protection against Relapse rate for individuals who succeeded treatment, E_rs
    if(strcmp(tran, "E_rs")==0){ 
        for(a=age; a<AG; a++)
			p[22][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 23 Protection against Proportion of mothers that infect their newborn children, E_mc
    if(strcmp(tran, "E_mc")==0){ 
        for(a=age; a<AG; a++)
			p[23][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 24 Protection against Diagnosis rate for pulmonary smear positive TB, E_cdr_i+ 
    if(strcmp(tran, "E_cdr_i+")==0){ 
        for(a=age; a<AG; a++)
			p[24][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 25 Protection against Diagnosis rate for pulmonary smear negative TB, E_cdr_i-
    if(strcmp(tran, "E_cdr_i-")==0){ 
        for(a=age; a<AG; a++)
			p[25][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }

    //Index 26 Protection against Diagnosis rate for non pulmonary TB, E_cdr_np
    if(strcmp(tran, "E_cdr_np")==0){ 
        for(a=age; a<AG; a++)
			p[26][a][1] = 1.0 - eff*exp(-waning*( (a - (age - shift))*5) );
    }
    
    
    //No index Protection against fast-progression AON, E_p_AON
	if(strcmp(tran, "E_p_AON")==0){
		eps_AON = 1;
        age_AON = age;
        cov_AON = cov;
    }

}

void sort(int *v, int n, int *index){
    int i,j;
    int tmp;
	for (i = 0; i < n; i++)                     //Loop for ascending ordering
	{
		for (j = 0; j < n; j++)             //Loop for comparing other values
		{
			if (v[j] > v[i])                //Comparing other array elements
			{
				tmp = v[i];         //Using temporary variable for storing last value
				v[i] = v[j];            //replacing value
				v[j] = tmp; 
                tmp = index[i];
                index[i] = index[j];
                index[j] = tmp;  
			}  
		}
	}
}

void load_contact_matrix(char file[string_length], double k[AG][AG][EG][EG])
{
	FILE *fich;
	double value;
	fich = fopen(file,"rt");
	double kaux[3][AG][AG][EG][EG];
	
	int i, j, l, m, n, scan;
	
	for(n=0;n<3;n++)//Cada matriz de las 3 que contiene el fichero.
	{
		for(i=0;i<AG;i++)//Grupos de edad=15
		{
			for(j=0;j<AG;j++)//Grupos de edad porque la matriz es cuadrada.
			{
				scan=fscanf(fich,"%lf ",&value);
				for(l=0;l<EG;l++)//Ramas de vacunación. Ni idea todavía.
				{
					for(m=0;m<EG;m++)//Repite porque cuadrada, supongo.
						kaux[n][i][j][l][m]=value;
				}
			}
			scan=fscanf(fich,"\n");
		}
		scan=fscanf(fich,"\n");
	}	
	fclose(fich);
	
	for(i=0;i<AG;i++)
	{
		for(j=0;j<AG;j++)
		{
			for(l=0;l<EG;l++)
			{
				for(m=0;m<EG;m++)
					k[i][j][l][m]=kaux[0][i][j][l][m];				
			}
		}
	}
	
}

void load_fitting_parameters(char fn[string_length], double gp[global_parameters][AG][EG], double rp[regional_parameters][AG][EG], double ax[5][AG]){
    FILE *fp;
    int scan, entries, i, j, k, a, e, index, amin, amax;
    double v1, v2, v3;
    char aux[500], type[500];
    fp = fopen(fn, "rt");
    scan = fscanf(fp, "entries = %d\n", &entries);
    scan = fscanf(fp, "Parameter Index Type Amin Amax Value LB UB\n");
    double phi_mu = 2.857143;

    for(j=0; j<5; j++){
        for(k=0; k<AG; k++){
            ax[j][k] = 0.0;
        }
    }
    for(i=0; i < entries; i++){
        scan = fscanf(fp, "%s %d %s %d %d %lf %lf %lf\n", aux, &index, type, &amin, &amax, &v1, &v2, &v3);
        amin=(int)(amin/age_window);
		amax=(int)(amax/age_window);
        if(amax>=AG)
            amax = AG;
        if(strcmp(type, "Regional")==0){
            for (a=amin; a<amax; a++){
			    for(e=0; e<EG; e++){
				    rp[index][a][e] = v1*v3 + v2*(1 - v1);			
		        }
            }
        }
        else if(strcmp(type, "Global")==0){
            if(index==2){ //Ratio constante para mutb
                for(a=amin; a<amax; a++){
			        for(e=0; e<EG; e++){
                        gp[index][a][e] = v1*v3 + v2*(1 - v1);
			            gp[index+1][a][e] = gp[index][a][e] / (phi_mu);			
		            }
                }
            } 
            else{
                for(a=amin; a<amax; a++){
			        for(e=0; e<EG; e++){
			            gp[index][a][e] = v1*v3 + v2*(1 - v1);			
		            }
                }
            }
        }
        else if(strcmp(type, "Diag_force")==0){ //d0 == 0, beta0 == 1, beta1 == 2, d1 == 3, dist_est == 4
            for (a=amin; a<amax; a++){
		        ax[index][a] = v1*v3 + v2*(1 - v1);			
            }
        }
        else{
            printf("Error with parameters to fit\n");
            break;
        }
        
    } 
    fclose(fp);
}

void load_fitting_parameters_from_R(int *ind, int *cosa, int *a_min, int *a_max, double *valu, double *lb, double *ub, int entries, double gp[global_parameters][AG][EG], double rp[regional_parameters][AG][EG], double ax[5][AG]){
    int i, j, k, a, e, index, amin, amax;
    double phi_mu = 2.857143;
    int tipo;
    double v1, v2, v3;

    for(j=0; j<5; j++){
        for(k=0; k<AG; k++){
            ax[j][k] = 0.0;
        }
    }

    for(i=0; i<entries; i++){
        index = (int)ind[i];
        tipo = (int)cosa[i];
        amin = (int)(a_min[i]/age_window);
		amax = (int)(a_max[i]/age_window);
        v1 = valu[i];
        v2 = lb[i];
        v3 = ub[i];
        if(amax>=AG)
            amax = AG;
        if(tipo==1){ //REGIONALES
            for (a=amin; a<amax; a++){
			    for(e=0; e<EG; e++){
				    rp[index][a][e] = v1*v3 + v2*(1 - v1);			
		        }
            }
        }
        else if(tipo==0){ // GLOBALES
            if(index==2){ //Ratio constante para mutb
                for(a=amin; a<amax; a++){
			        for(e=0; e<EG; e++){
                        gp[index][a][e] = v1*v3 + v2*(1 - v1);
			            gp[index+1][a][e] = gp[index][a][e] / (phi_mu);			
		            }
                }
            } 
            else{
                for(a=amin; a<amax; a++){
			        for(e=0; e<EG; e++){
			            gp[index][a][e] = v1*v3 + v2*(1 - v1);			
		            }
                }
            }
        }
        else if(tipo==2){ //d0 == 0, beta0 == 1, beta1 == 2, d1 == 3, dist_est == 4
            for(a=amin; a<amax; a++){
		        ax[index][a] = v1*v3 + v2*(1 - v1);	
            }
        }
        else{
            printf("Error with parameters to fit\n");
            break;
        }
        
    } 
}

void load_window(char window[string_length])
{
	FILE *fp1;
	int scan;
	fp1=fopen(window,"r");
	scan=fscanf(fp1,"%d %d %d %d %d\n", &year_min, &year_max, &year_end, &time_stabilization, &year_vac);
	fclose(fp1);
}

double funcion_param_t(double l0, double r0, double t, double lsup, double linf){
	double lt;
    if(r0>1E-15){
		lt = l0 + (lsup - l0)*t/(t + 1.0/r0);
    }
	else
	{
		if(r0<-1E-15)
			lt = l0 + (linf - l0)*t/(t - 1.0/r0);
		else
			lt = l0;
	}
	
	return lt;
}

void start_record(double macro[observables+annuals],double s[states][AG][EG],double S[states]){
	int i;
	for(i=0; i<observables; i++){
		macro[i] = S[reservories+i];
	}
}

int end_record(double macro[observables+annuals],double s[states][AG][EG],double S[states],double slope_calcs[2][2],int reg)
{
	double slope[2];
	int i;
	for(i=1; i<observables; i++){//Incrementales: incidencia, mortandad, vacunacion e infecciones
		macro[i]=(-macro[i]+S[reservories+i])/(((macro[0]+S[reservories])/2.0)/scaling);//inc variable 1
    }
	
	slope_calcs[1][0] = slope_calcs[0][0];
	slope_calcs[1][1] = slope_calcs[0][1];
#if incidence_to_fit==0
	slope_calcs[0][0] = macro[1];
#else	
	slope_calcs[0][0] = macro[2];
#endif
	slope_calcs[0][1] = macro[3];	
	if(reg==0)
		return(0);
	else
	{
		slope[0] = (slope_calcs[0][0] - slope_calcs[1][0])/(double)record_window;
		slope[1] = (slope_calcs[0][1] - slope_calcs[1][1])/(double)record_window;
		
		
		if(fabs(slope[0])<threshold_slope && fabs(slope[1])<threshold_slope && reg>min_reg_fit)
			return(1);		
		else
			return(0);		
	}
}

void calculate_derivatives(int a, int e, double ds[states][AG][EG], double f[fluxes])
{
	int i;
	
	ds[0][a][e] = f[0]-f[1]-f[2]-f[41]; //Susceptibles
	ds[1][a][e] = f[1]+f[35]+f[39]-f[42]-f[3]-f[4]-f[5]+f[29]+f[30]+f[31]+f[32]+f[33]+f[34]+f[64]+f[65]+f[66]; //Fast-latency
    //ds[1][a][e] = f[1]+f[35]+f[39]-f[42]-f[3]-f[4]-f[5]+f[64]+f[65]+f[66]; //Fast-latency
	//ds[2][a][e] = f[2]-f[35]+f[40]-f[43]-f[6]-f[7]-f[8]+f[72]+f[73]+f[74]+f[75]+f[76]+f[77]+f[67]+f[68]+f[69]; //Slow-Latency

    ds[2][a][e] = f[2]-f[35]+f[40]-f[43]-f[6]-f[7]-f[8]+f[67]+f[68]+f[69]; //Slow-Latency

	ds[3][a][e] = f[3]+f[6]-f[36]-f[44]-f[9]+f[23]+f[26]+f[21]-f[12]+f[61]; //Active disease smear-positive
	ds[4][a][e] = f[4]+f[7]-f[37]-f[45]-f[10]+f[24]+f[27]-f[21]-f[13]+f[62]; //Active disease smear-negative
	ds[5][a][e] = f[5]+f[8]-f[38]-f[46]-f[11]+f[25]+f[28]-f[14]+f[63]; //Active disease extra-pulmonary	
	ds[6][a][e] = f[9]-f[15]-f[18]-f[47]+f[22]-f[56]-f[57];  //Under treatment smear-positive
	ds[7][a][e] = f[10]-f[16]-f[19]-f[48]-f[22]-f[58]-f[70]; //Under treatment smear-negative
	ds[8][a][e] = f[11]-f[17]-f[20]-f[49]-f[59]-f[71]; //Under treatment extra-pulmonary
	ds[9][a][e] = f[15]-f[23]-f[29]-f[50];//   -f[72]  ;  //Default (smear-positive)
	ds[10][a][e] = f[12]-f[26]-f[30]-f[51];//  -f[73]  ; //Natural recovered (smear-positive)
	ds[11][a][e] = f[16]-f[24]-f[31]-f[52];//  -f[74]  ; //Default (smear-negative)
	ds[12][a][e] = f[13]-f[27]-f[32]-f[53];//  -f[75]  ; //Natural recovered (smear-negative)
	ds[13][a][e] = f[17]-f[25]-f[33]-f[54];//  -f[76]  ; //Default (non-pulmonary)
	ds[14][a][e] = f[14]-f[28]-f[34]-f[55];//  -f[77]  ; //Natural recovered (non-pulmonary)	
	ds[15][a][e] = f[18]-f[61]-f[64]-f[67]; //Success (smear-positive)
	ds[16][a][e] = f[19]-f[62]-f[65]-f[68]; //Success (smear-negative)
	ds[17][a][e] = f[20]-f[63]-f[66]-f[69]; //Success (non-pulmonary)
	ds[18][a][e] = f[56]-f[60]+f[70]+f[71]; //Failure	
	
	ds[reservories][a][e] = f[0]+f[39]+f[40]-f[36]-f[37]-f[38];
	for(i=41;i<56;i++)
		ds[reservories][a][e] -= f[i];
	ds[reservories+1][a][e] = f[3]+f[4]+f[5]+f[6]+f[7]+f[8]+f[23]+f[24]+f[25]+f[26]+f[27]+f[28]+f[61]+f[62]+f[63];
	ds[reservories+2][a][e] = f[9]+f[10]+f[11];
	ds[reservories+3][a][e] = f[36]+f[37]+f[38]+f[57]+f[58]+f[59]+f[60];
    ds[reservories+5][a][e] = 0.0;
    ds[reservories+4][a][e] = 0.0;
    ds[reservories+6][a][e] = 0.0;
    ds[reservories+7][a][e] = 0.0;
    ds[reservories+8][a][e] = 0.0;
}

void correct_contact_matrix(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]) //M3 density + normalization method.
{
	int a1, a2, e1;
	double total_N, N_t[AG];
	double sum_dem = 0;
	
	total_N = 0;
	for(a1=0; a1<AG; a1++){
		N_t[a1] = 0;
		for(e1=0; e1<EG; e1++){
			total_N += s[reservories][a1][e1];
			N_t[a1] += s[reservories][a1][e1];
		}		
	}
	
	for(a1=0; a1<AG; a1++){
		for(a2=0; a2<AG; a2++){
			sum_dem += cm[a1][a2][0][0]*N_t[a2]*N_t[a1];
		}
	}
	
	for(a1=0; a1<AG; a1++){
		for(a2=0; a2<AG; a2++){
			P_matrix[a1][a2] = 2*cm[a1][a2][0][0]*N_t[a2]*total_N/sum_dem;
		}
	}	
}

void load_cm_modified(char file[string_length], double k[AG][AG][EG][EG])
{
	FILE *fich;
	double value;
    fich = fopen(file,"rt");
	//fich = fopen("Ref_matrix.txt", "rt");
	int i, j, l, m, scan;
	
    for(i=0;i<AG;i++)//AG=15
	{
		for(j=0;j<AG;j++)
		{
			scan=fscanf(fich,"%lf ",&value);
			for(l=0;l<EG;l++)
			{
				for(m=0;m<EG;m++)
					k[i][j][l][m]=value;
			}
		}
		scan=fscanf(fich,"\n");
	}

	fclose(fich);	
}	

void cm_M0(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]){ //M0 and M1 work only with intensive contact matrices (M)
	int a1, a2;
	for(a1=0; a1<AG; a1++){
		for(a2=0; a2<AG; a2++){
			P_matrix[a1][a2] = cm[a1][a2][0][0];
		}
	}	
}

void cm_M1(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]){
	int a1,a2,e1;
	double N_t[AG];
	
	for(a1=0;a1<AG;a1++){
		N_t[a1]=0;
		for(e1=0;e1<EG;e1++){
			N_t[a1]+=s[reservories][a1][e1];
		}		
	}
	
	for(a1=0;a1<AG;a1++){
		for(a2=0;a2<AG;a2++){
			P_matrix[a1][a2] = (cm[a1][a2][0][0]*N_t[a1] + cm[a2][a1][0][0]*N_t[a2])/(2*N_t[a1]);
		}
	}	
}

void cm_M2(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]){ //cm_M2 and cm_M3 work only with Gamma matrices (_pi).
	int a1,a2,e1;
	double N_t[AG];
    double total_N;
	
	total_N=0;
	for(a1=0;a1<AG;a1++){
		N_t[a1]=0;
		for(e1=0;e1<EG;e1++){
			total_N += s[reservories][a1][e1];
			N_t[a1] += s[reservories][a1][e1];
		}		
	}
	
	for(a1=0;a1<AG;a1++){
		for(a2=0;a2<AG;a2++){
			P_matrix[a1][a2] = cm[a1][a2][0][0]*N_t[a2]/total_N;
		}
	}	
}

void cm_M3(double cm[AG][AG][EG][EG], double P_matrix[AG][AG]){
	int a1, a2, e1;
	double total_N, N_t[AG];
	double sum_dem = 0;
	
	total_N = 0;
	for(a1=0; a1<AG; a1++){
		N_t[a1] = 0;
		for(e1=0; e1<EG; e1++){
			total_N += s[reservories][a1][e1];
			N_t[a1] += s[reservories][a1][e1];
		}		
	}
	
	for(a1=0; a1<AG; a1++){
		for(a2=0; a2<AG; a2++){
			sum_dem += cm[a1][a2][0][0]*N_t[a2]*N_t[a1];
		}
	}
	
	for(a1=0; a1<AG; a1++){
		for(a2=0; a2<AG; a2++){
			P_matrix[a1][a2] = 2*cm[a1][a2][0][0]*N_t[a2]*total_N/sum_dem;
		}
	}	
}


