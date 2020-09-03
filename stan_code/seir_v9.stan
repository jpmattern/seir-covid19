functions {
    vector compute_b_spline(int nx, vector x, real[] knots_padded, int p, int i);
    vector compute_b_spline(int nx, vector x, real[] knots_padded, int p, int i){
        //////////////////////////////////////////
        // function for computing b-splines recursively
        vector[nx] alpha1;
        vector[nx] alpha2;

        if (p == 0){
            for (ix in 1:nx){
                if (x[ix] < knots_padded[i] || x[ix] >= knots_padded[i+1]){
                    alpha1[ix] = 0.0;
                } else {
                    alpha1[ix] = 1.0;
                }
            }
            return alpha1;
        }
        if (knots_padded[i+p] == knots_padded[i]){
            alpha1 = rep_vector(0.0, nx);
        } else {
            for (ix in 1:nx){
                alpha1[ix] = (x[ix]-knots_padded[i])/(knots_padded[i+p]-knots_padded[i]);
            }
        }
        if (knots_padded[i+p+1] == knots_padded[i+1]){
            alpha2 = rep_vector(0.0, nx);
        } else {
            for (ix in 1:nx){
                alpha2[ix] = (knots_padded[i+p+1]-x[ix])/(knots_padded[i+p+1]-knots_padded[i+1]);
            }
        }
        return alpha1 .* compute_b_spline(nx, x, knots_padded, p-1, i) + alpha2 .* compute_b_spline(nx, x, knots_padded, p-1, i+1);
    }
    matrix compute_b_splines(int nx, vector x, int n, real[] knots, int p){
        //////////////////////////////////////////
        // function for pre-computing b-splines for given knots  
        matrix[nx,n+p-1] result;
        real knots_padded[n+2*p];

        int i0 = 1;
        int i1 = 1;
        for (i in 1:p){
            knots_padded[i] = knots[1];
            knots_padded[p+n+i] = knots[n];
        }
        for (i in 1:n){
            knots_padded[p+i] = knots[i];
        }
        for (i in 1:n+p-1){
            result[:,i] = compute_b_spline(nx, x, knots_padded, p, i);
        }
        result[nx,n+p-1] = 1.0;
        return result;
    }
}
data {
    //////////////////////////////////////////
    // data required to run model  
    int<lower=0> nobs_I;                // number of timesteps with observations
    int tobs_I[nobs_I];                 // associated obs times
    real<lower=0> obs_I[nobs_I];        // observations cumulative number of infected

    int<lower=0> nobs_Rmort;            // number of timesteps with number of deaths observations
    int tobs_Rmort[nobs_Rmort];         // associated obs times
    real<lower=0> obs_Rmort[nobs_Rmort];// observations cumulative number of deaths

    int<lower=0> nobs_Hmod;             // number of timesteps with non-ICU hospitalization observations
    int tobs_Hmod[nobs_Hmod];           // associated obs times
    real<lower=0> obs_Hmod[nobs_Hmod];  // observations current number of non-ICU hospitalized

    int<lower=0> nobs_Hicu;             // number of timesteps with ICU hospitalization observations
    int tobs_Hicu[nobs_Hicu];           // associated obs times
    real<lower=0> obs_Hicu[nobs_Hicu];  // observations current number of ICU hospitalized
    
    int<lower=0> nt;                    // number of time steps
    int<lower=1> nage;                  // number of age groups in the model
    vector[nage] npop;                  // total population for each age group
    vector[nage] frac_hosp;             // ICU + non-ICU
    vector[nage] frac_icu;              // ICU
    vector[nage] frac_mort;             // mortality for each age group
    vector[nage] frac_asym;             // fraction of asymptomatic
    
    //////////////////////////////////////////
    // prior parameter distributions
    real<lower=1.0> mu_duration_lat;        // mean duration in "exposed" stage
    real<lower=0.0> sigma_duration_lat;     // sd duration in "exposed" stage
    real<lower=1.0> mu_duration_rec_asym;   // mean duration in "infectious" stage for asymptomatic cases
    real<lower=0.0> sigma_duration_rec_asym;// sd duration in "infectious" stage for asymptomatic cases
    real<lower=1.0> mu_duration_rec_mild;   // mean duration in "infectious" stage for mild cases
    real<lower=0.0> sigma_duration_rec_mild;// sd duration in "infectious" stage for mild cases
    real<lower=1.0> mu_duration_pre_hosp;   // mean duration in "infectious" stage for hospitalized cases
    real<lower=0.0> sigma_duration_pre_hosp;// sd duration in "infectious" stage for hospitalized cases
    real<lower=1.0> mu_duration_hosp_mod;   // mean duration in hospital for non-ICU cases
    real<lower=0.0> sigma_duration_hosp_mod;// sd duration in hospital for non-ICU cases
    real<lower=1.0> mu_duration_hosp_icu;   // mean duration in hospital for ICU cases
    real<lower=0.0> sigma_duration_hosp_icu;// sd duration in hospital for ICU cases
    
    real<lower=0.0> mu_beta1;       // mean initial beta estimate
    real<lower=0.0> sigma_beta1;    // sd initial beta estimate
    
    real<lower=0.0> lambda_ini_exposed; // parameter for initial conditions of "exposed"
   
    // hierarchical priors for standard deviation of model observation misfit for each observation type
    real<lower=0.0> lambda_Iobs;    // parameter for Iobs misfit
    real<lower=0.0> lambda_Hmod;    // parameter for Hmod misfit
    real<lower=0.0> lambda_Hicu;    // parameter for Hicu misfit
    real<lower=0.0> lambda_Rmort;   // parameter for Rmort misfit

    real<lower=0.0> alpha_multiplier; 
    
    //////////////////////////////////////////
    // spline-related + future interventions
    int<lower=0> dknot;                       // distance in time between the knots used to construct the splines 
    int<lower=0> itoday;                      // time index for todays date at which to stop estimating beta
    int<lower=0> ninter;                      // number of interventions
    int<lower=0> t_inter[ninter];             // start time of each interventions
    int<lower=0> len_inter[ninter];           // length of each intervention
    real<lower=0.0> mu_beta_inter[ninter];    // mean change in beta through intervention 
    real<lower=0.0> sigma_beta_inter[ninter]; // sd change in beta through intervention
    // splinemode = 1: estimate up to itoday but enforce derivative of 0
    // splinemode = 2: constant value for beta in [itoday-dknot,itoday]
    int<lower=1,upper=2> splinemode;

    //////////////////////////////////////////
    // spline of order p_fractest for fraction tested
    int nknots_fractest;                                // number of knots for fraction tested splines
    int p_fractest;                                     // degree of fraction tested splines
    real knots_fractest[nknots_fractest];               // knots for fraction tested splines
    real mu_fractest[nknots_fractest+p_fractest-1];     // mean of fraction tested a spline control points
    real sigma_fractest[nknots_fractest+p_fractest-1];  // sd of fraction tested a spline control points
}
transformed data {
    //assigning indices for state matrix x
    int S = 1;
    int E = 2;
    int Iasym = 3;
    int Imild = 4;
    int Ipreh = 5; 
    int Hmod  = 6;
    int Hicu  = 7;
    int Rlive = 8;
    int Rmort = 9;

    int p = 3;

    matrix<lower=0.0>[5,nage] alpha_frac_I;
    
    real npoptotal = sum(npop);
    vector[nt] t;
    
    // without splinemode
    //int nbeta_est = (itoday-1)/dknot + 2;
    //int nbeta = nbeta_est + 1 + ninter*4;

    // with splinemode
    int nbeta_est = (itoday-1)/dknot + 3 - splinemode;
    int nbeta = nbeta_est + ninter*4 + splinemode;

    real knots_beta[nbeta];
    matrix[nt,nbeta+p-1] bsplines;

    matrix[nt,nknots_fractest+p_fractest-1] bsplines_fractest;

    // the alpha parameter for the Dirichlet distribution  
    for (iage in 1:nage){
        // 1: asym, 2: mild, 3: hosp (non-ICU), 4: non-mort ICU, 5: mort ICU 
        alpha_frac_I[1,iage] = alpha_multiplier * frac_asym[iage];
        alpha_frac_I[3,iage] = alpha_multiplier * (frac_hosp[iage]-frac_icu[iage]);
        alpha_frac_I[4,iage] = alpha_multiplier * (frac_icu[iage]-frac_mort[iage]);
        alpha_frac_I[5,iage] = alpha_multiplier * frac_mort[iage];
        // rest defaults to mild
        alpha_frac_I[2,iage] = alpha_multiplier - alpha_frac_I[1,iage] - alpha_frac_I[3,iage] - alpha_frac_I[4,iage] - alpha_frac_I[5,iage];
    }
    
    //////////////////////////////////////////
    // compute knots for beta b-splines
    
    print("itoday/dknot = ",itoday/dknot,", nbeta = ", nbeta)
    knots_beta[1] = 1;
    {
        int it = dknot;
        int ibeta = 2;
        real alpha_shape = 0.1; // must be between 0.0 and 0.5
        while (it < itoday){
            knots_beta[ibeta] = it;
            ibeta += 1;
            it += dknot;
        }
        knots_beta[ibeta-1] = itoday - dknot;
        knots_beta[ibeta] = itoday;
        ibeta += 1;
        
        // set knots for each intervention
        for (iinter in 1:ninter){
            knots_beta[ibeta] = t_inter[iinter];
            if (knots_beta[ibeta] <  knots_beta[ibeta-1]){
                if (iinter == 1){
                    reject("Interventions are not allowed to be in the past (problem with intervention ",iinter,").")
                }
                reject("Interventions are not allowed to overlap (problem with intervention ",iinter,").")
            }
            knots_beta[ibeta+1] = t_inter[iinter]+alpha_shape*len_inter[iinter];
            knots_beta[ibeta+2] = t_inter[iinter]+(1.0-alpha_shape)*len_inter[iinter];
            knots_beta[ibeta+3] = t_inter[iinter]+len_inter[iinter];
            ibeta += 4;
        }
    }
    
    t[1] = 1.0;
    for (i in 2:nt){
        t[i] = t[i-1] + 1.0;
    }
    
    knots_beta[nbeta] = nt;
    print("knots_beta = ", knots_beta)
    
    // compute b-splines, now that knots are known
    bsplines = compute_b_splines(nt, t, nbeta, knots_beta, p);
    //print("bsplines = ", bsplines)

    //////////////////////////////////////////
    // compute b-splines for fraction tested
    if (knots_fractest[1] != 1){
        reject("First entry in knots_fractest must be 1.")
    }
    if (knots_fractest[nknots_fractest] != nt){
        reject("Last entry in knots_fractest must be nt.")
    }
    bsplines_fractest = compute_b_splines(nt, t, nknots_fractest, knots_fractest, p_fractest);
}
parameters {
    //specifying parameters bounds (not always necessary) for the parameters that are being estimated
    real<lower=0.0, upper=10.0> beta_est[nbeta_est];
    real<lower=0.0> sigma_beta;
    real<lower=0.0> beta_change[ninter];
    
    real<lower=1.0> duration_lat; // duration is a minimum of 1 which is the stepsize of this model
    real<lower=1.0> duration_rec_asym;
    real<lower=1.0> duration_rec_mild;
    real<lower=1.0> duration_pre_hosp;
    real<lower=1.0> duration_hosp_mod;
    real<lower=1.0> duration_hosp_icu;

    simplex[5] frac_I[nage];

    vector<lower=0.0, upper=1.0>[nknots_fractest+p_fractest-1] control_fractest;

    real<lower=0> ini_exposed;
    
    real<lower=0> sigma_Iobs;
    real<lower=0> sigma_Rmort;
    real<lower=0> sigma_Hmod;
    real<lower=0> sigma_Hicu;
}
transformed parameters {
    matrix<lower=0.0>[9,nt] x;
    vector[nt] beta;
    vector[nbeta+p-1] control_beta;
    real mu_Iobs[nt];
    real Rt[nt];
    vector[nt] fractest = bsplines_fractest * control_fractest;
    {
        // variables in curly brackets will not have output, they are local variables
        int ibetaoffset;

        real newE;
        real newI;
        real newI_sum;
        real newrec_asym;
        real newrec_mild;
        real newrec_mod;
        real newrec_icu;
        real newhosp_mod;
        real newhosp_icu;
        real newhosp_mort;
        real newmort;
        
        real Ipreh_mod_new;
        real Ipreh_icu_new;
        real Ipreh_mort_new;
        real Hicu_live_new;
        real Hicu_mort_new;
        
        real I_cur;
        vector[nage] S_new;
        vector[nage] E_new;
        vector[nage] S_cur;
        vector[nage] E_cur;
        
        // initial cond
        real Ipreh_mod_cur = 0.0;
        real Ipreh_icu_cur = 0.0;
        real Ipreh_mort_cur = 0.0;
        real Hicu_live_cur = 0.0;
        real Hicu_mort_cur = 0.0;
        
        Ipreh_mod_new = Ipreh_mod_cur;
        Ipreh_icu_new = Ipreh_icu_cur;
        Ipreh_mort_new = Ipreh_mort_cur;
        
        x[:,1] = rep_vector(0.0, 9); //: means all entries. puts a zero in x1-9 for initial entries except for SEI.

        x[S,1] = npoptotal-ini_exposed;
        x[E,1] = ini_exposed;
        
        I_cur = 0.0;
        E_cur = npop*ini_exposed/npoptotal;
        S_cur = npop - E_cur;
        
        mu_Iobs[1] = I_cur * fractest[1]; //defining first value of mu_Iobs for first time steps based on I true * fraction tested
       
        //////////////////////////////////////////
        // fill control_beta
        for (ibeta in 1:nbeta_est){
            control_beta[ibeta] = beta_est[ibeta];
        }
        if (splinemode == 2){
            control_beta[nbeta_est+1] = beta_est[nbeta_est];
            ibetaoffset = nbeta_est + 1;
        } else {
            ibetaoffset = nbeta_est;
        }
        for (iinter in 1:ninter){
            control_beta[ibetaoffset+4*iinter-3] = control_beta[ibetaoffset+4*iinter-4];
            control_beta[ibetaoffset+4*iinter-2] = control_beta[ibetaoffset+4*iinter-4];
            control_beta[ibetaoffset+4*iinter-1] = control_beta[ibetaoffset+4*iinter-4];
            control_beta[ibetaoffset+4*iinter] = control_beta[ibetaoffset+4*iinter-1] * beta_change[iinter];
        }
        for (ibeta in nbeta:nbeta+p-1){
            control_beta[ibeta] = control_beta[nbeta-1];
        }
        beta = bsplines * control_beta;
        
        //////////////////////////////////////////
        // the SEIR model

        for (it in 1:nt-1){
            //////////////////////////////////////////
            // set transition variables 
            newrec_asym = 1.0/duration_rec_asym * x[Iasym,it];
            newrec_mild = 1.0/duration_rec_mild * x[Imild,it];
            newrec_mod = 1.0/duration_hosp_mod * x[Hmod,it];
            newrec_icu = 1.0/duration_hosp_icu * Hicu_live_cur;
            newmort = 1.0/duration_hosp_icu * Hicu_mort_cur;
            
            newhosp_mod = 1.0/duration_pre_hosp * Ipreh_mod_cur;
            newhosp_icu = 1.0/duration_pre_hosp * Ipreh_icu_cur;
            newhosp_mort = 1.0/duration_pre_hosp * Ipreh_mort_cur;
            
            x[Iasym,it+1] = x[Iasym,it];
            x[Imild,it+1] = x[Imild,it];
            Ipreh_mod_new = Ipreh_mod_cur;    // step not necessary, they are already equal
            Ipreh_icu_new = Ipreh_icu_cur;    // step not necessary, they are already equal
            Ipreh_mort_new = Ipreh_mort_cur;  // step not necessary, they are already equal

            //////////////////////////////////////////
            // S -> E -> I
           
            newI_sum = 0.0;
            Rt[it] = 0.0;
            for (iage in 1:nage){
                newE = fmin(S_cur[iage], beta[it]/npop[iage] * I_cur * S_cur[iage]);
                newI = 1/duration_lat * E_cur[iage]; // duration_lat > 1
                
                S_new[iage] = S_cur[iage] - newE;
                E_new[iage] = E_cur[iage] + newE - newI;
                
                // 1: asym, 2: mild, 3: hosp (non-ICU), 4: non-mort ICU, 5: mort ICU 
                x[Iasym,it+1]  = x[Iasym,it+1]  + frac_I[iage][1] * newI;
                x[Imild,it+1]  = x[Imild,it+1]  + frac_I[iage][2] * newI;
                Ipreh_mod_new  = Ipreh_mod_new  + frac_I[iage][3] * newI;
                Ipreh_icu_new  = Ipreh_icu_new  + frac_I[iage][4] * newI;
                Ipreh_mort_new = Ipreh_mort_new + frac_I[iage][5] * newI;
                
                newI_sum = newI_sum + newI;

                // compute Rt
                // NOTE requires age structure data for S and time-dependent beta 
                // and is hence not easy to calculate in generated quantities block
                Rt[it] = Rt[it] + S_cur[iage]/npoptotal * beta[it] * (frac_I[iage][1]*duration_rec_asym + frac_I[iage][2]*duration_rec_mild + (1.0-frac_I[iage][1]-frac_I[iage][2])*duration_pre_hosp);
            }

            //////////////////////////////////////////
            // I -> recovery or hospital
           
            x[Iasym,it+1]  = x[Iasym,it+1]  - newrec_asym;
            x[Imild,it+1]  = x[Imild,it+1]  - newrec_mild;
            Ipreh_mod_new  = Ipreh_mod_new  - newhosp_mod;
            Ipreh_icu_new  = Ipreh_icu_new  - newhosp_icu;
            Ipreh_mort_new = Ipreh_mort_new - newhosp_mort;
            
            x[Ipreh,it+1] = Ipreh_mod_new + Ipreh_icu_new + Ipreh_mort_new;
            
            //////////////////////////////////////////
            // in hospital
            
            x[Hmod,it+1]  = x[Hmod,it]    + newhosp_mod  - newrec_mod;
            Hicu_live_new = Hicu_live_cur + newhosp_icu  - newrec_icu;
            Hicu_mort_new = Hicu_mort_cur + newhosp_mort - newmort;
            
            x[Hicu,it+1] = Hicu_live_new + Hicu_mort_new;

            //////////////////////////////////////////
            // recovery/death

            x[Rlive,it+1] = x[Rlive,it] + newrec_asym + newrec_mild + newrec_mod + newrec_icu;
            x[Rmort,it+1] = x[Rmort,it] + newmort;
            
            //////////////////////////////////////////
            // observed infected
            
            mu_Iobs[it+1] = mu_Iobs[it] + newI_sum * fractest[it+1];
            
            x[S,it+1] = sum(S_new);
            x[E,it+1] = sum(E_new);
            
            //////////////////////////////////////////
            // test

            if (fabs(sum(x[:,it+1])-npoptotal) > 1e-1){
                reject("Model is leaking, net gain: ", sum(x[:,it+1])-npoptotal) 
            }

            Ipreh_mod_cur = Ipreh_mod_new;
            Ipreh_icu_cur = Ipreh_icu_new;
            Ipreh_mort_cur = Ipreh_mort_new;
            Hicu_live_cur = Hicu_live_new;
            Hicu_mort_cur = Hicu_mort_new;

            S_cur = S_new;
            E_cur = E_new;
            I_cur = x[Iasym,it+1] + x[Imild,it+1] + x[Ipreh,it+1];
        }
        // compute Rt for last timestep
        Rt[nt] = 0.0;
        for (iage in 1:nage){
            Rt[nt] = Rt[nt] + S_cur[iage]/npoptotal * beta[nt] * (frac_I[iage][1]*duration_rec_asym + frac_I[iage][2]*duration_rec_mild + (1.0-frac_I[iage][1]-frac_I[iage][2])*duration_pre_hosp);
        }
    }
}
model {
    //////////////////////////////////////////
    // prior distributions
    sigma_beta ~ exponential(30.0);
    beta_est[1] ~ normal(mu_beta1, sigma_beta1);
    for (ibeta in 2:nbeta_est){
        beta_est[ibeta] ~ normal(beta_est[ibeta-1], sigma_beta);
    }
    beta_change ~ normal(mu_beta_inter, sigma_beta_inter);
    
    duration_lat ~ normal(mu_duration_lat, sigma_duration_lat);
    duration_rec_asym ~ normal(mu_duration_rec_asym, sigma_duration_rec_asym);
    duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild); 
    duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
    duration_hosp_mod ~ normal(mu_duration_hosp_mod, sigma_duration_hosp_mod);
    duration_hosp_icu ~ normal(mu_duration_hosp_icu, sigma_duration_hosp_icu);

    ini_exposed ~ exponential(lambda_ini_exposed);

    for (iage in 1:nage){
        frac_I[iage] ~ dirichlet(alpha_frac_I[:,iage]);
    }

    control_fractest ~ normal(mu_fractest, sigma_fractest);

    //////////////////////////////////////////
    // fitting observations
    
    sigma_Iobs ~ exponential(lambda_Iobs);
    sigma_Rmort ~ exponential(lambda_Rmort);
    
    sigma_Hmod ~ exponential(lambda_Hmod);
    sigma_Hicu ~ exponential(lambda_Hicu);
    {
        real tmp;
        for (iobs in 1:nobs_I){
            tmp = (obs_I[iobs]-mu_Iobs[tobs_I[iobs]])/sigma_Iobs;
            tmp ~ std_normal();
        }    
        for (iobs in 1:nobs_Rmort){
            tmp = (obs_Rmort[iobs]-x[Rmort,tobs_Rmort[iobs]])/sigma_Rmort;
            tmp ~ std_normal();
        }
        for (iobs in 1:nobs_Hmod){
            tmp = (obs_Hmod[iobs]-x[Hmod,tobs_Hmod[iobs]])/sigma_Hmod;
            tmp ~ std_normal();
        }
        for (iobs in 1:nobs_Hicu){
            tmp = (obs_Hicu[iobs]-x[Hicu,tobs_Hicu[iobs]])/sigma_Hicu;
            tmp ~ std_normal();
        }
    }
}
generated quantities{
    real fractested[nobs_I];
    real hospitalized[nt];
    real obs_I_sim[nt];
    for (iobs in 1:nobs_I){
        // dividing here by the cumulative number as well 
        // mu_Iobs[tobs_I[iobs]]/fractest[tobs_I[iobs]] is cumulative number of infectious
        // resulting quantity is the cumulative number of observed infected divided by the cumulative number of infectious in the model
        fractested[iobs] = obs_I[iobs]/(mu_Iobs[tobs_I[iobs]]/fractest[tobs_I[iobs]]);
    }
    for (it in 1:nt){
        hospitalized[it] = x[Hmod,it] + x[Hicu,it];
        obs_I_sim[it] = normal_rng(mu_Iobs[it], sigma_Iobs);
    }
}
