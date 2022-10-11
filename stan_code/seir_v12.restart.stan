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
    int<lower=0> nt;                    // number of time steps
    int<lower=1> nage;                  // number of age groups in the model

    int<lower=0> nobs_I;                // number of timesteps with observations
    int tobs_I[nobs_I];                 // associated obs times
    int<lower=0> obs_I[nage,nobs_I];    // observations cumulative number of infected

    int<lower=0> nobs_Rmort;            // number of timesteps with number of deaths observations
    int tobs_Rmort[nobs_Rmort];         // associated obs times
    real<lower=0> obs_Rmort[nobs_Rmort];// observations cumulative number of deaths

    int<lower=0> nobs_H;                // number of timesteps with non-ICU hospitalization observations
    int tobs_H[nobs_H];                 // associated obs times
    real<lower=0> obs_H[nobs_H];        // observations current number of non-ICU hospitalized
    real<lower=0> obs_Hicu[nobs_H];     // observations current number of ICU hospitalized

    int<lower=0> nvac;                  // number of timesteps with immunizations due vaccination
    int<lower=1> tvac[nvac];            // associated times
    int<lower=0> vac[nage,nvac];        // number of vaccination immunizationsfor each age group

    vector[nage] npop;                  // total population for each age group
    vector[nage] age;                   // a representative age for each age group
    vector[nage] frac_asym;             // fraction of asymptomatic

    //////////////////////////////////////////
    // initial conditions
    matrix[nage, 9] mu_x_ini_2to10;     // mean
    matrix[nage, 9] sigma_x_ini_2to10;  // sd

    vector<lower=0.0>[nage] V_ini;      // total number immunized through vaccination at start

    vector<lower=0.0>[nage] mu_Iobs_ini;

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

    real<upper=0.0> mu_frac_hospmild_i45;
    real<lower=0.0> sigma_frac_hospmild_i45;
    real<upper=0.0> mu_frac_iculive_i45;
    real<lower=0.0> sigma_frac_iculive_i45;
    real<upper=0.0> mu_frac_mort_i45;
    real<lower=0.0> sigma_frac_mort_i45;
    real<lower=0.0> mu_frac_hospmild_slope;
    real<lower=0.0> sigma_frac_hospmild_slope;
    real<lower=0.0> mu_frac_iculive_slope;
    real<lower=0.0> sigma_frac_iculive_slope;
    real<lower=0.0> mu_frac_mort_slope;
    real<lower=0.0> sigma_frac_mort_slope;

    vector[nage] mu_beta1;       // mean initial beta estimate
    real<lower=0.0> sigma_beta1;    // sd initial beta estimate

    // hierarchical priors for standard deviation of model observation misfit for each observation type
    real<lower=0.0> lambda_Iobs;    // parameter for Iobs misfit
    real<lower=0.0> lambda_H;       // parameter for H misfit
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

    //////////////////////////////////////////
    // spline of order p_imports
    // NOTE that imports is using the same knots as fractest
    int p_imports;                                      // degree of imports splines
    real mu_imports[nknots_fractest+p_imports-1];       // mean of imports
    real sigma_imports[nknots_fractest+p_imports-1];    // sd of imports

    //////////////////////////////////////////
    // spline of order p_age for better age interpolation
    int nknots_age;
    real knots_age[nknots_age];

    // dev only
    int<lower=0, upper=1> prior_only[4];
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
    int Vsus = 10;
    int Vrec = 11;

    int nstate = 11;

    int p = 3;

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

    matrix[nt,nknots_fractest+p_imports-1] bsplines_imports;

    int p_age = 3;
    matrix[nage,nknots_age+p_age-1] bsplines_age;

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

    bsplines_imports = compute_b_splines(nt, t, nknots_fractest, knots_fractest, p_imports);

    //////////////////////////////////////////
    // compute b-splines for age-dependence
    bsplines_age = compute_b_splines(nage, age, nknots_age, knots_age, p_age);
}
parameters {
    //specifying parameters bounds (not always necessary) for the parameters that are being estimated
    matrix<lower=0.0>[nage, 9] x_ini_2to10;

    real<lower=0.0, upper=1.0> beta_est[nknots_age,nbeta_est];
    real<lower=0.0> sigma_beta;
    real<lower=0.0> beta_change[ninter];

    real<lower=1.0> duration_lat; // duration is a minimum of 1 which is the stepsize of this model
    real<lower=1.0> duration_rec_asym;
    real<lower=1.0> duration_rec_mild;
    real<lower=1.0> duration_pre_hosp;
    real<lower=1.0> duration_hosp_mod;
    real<lower=1.0> duration_hosp_icu;

    // added a lower limit > 0 on fraction tested
    vector<lower=0.03, upper=1.0>[nknots_fractest+p_fractest-1] control_fractest;

    vector<lower=0.01>[nknots_fractest+p_imports-1] control_imports;

    real<lower=0, upper=100> sigma_Iobs;
    real<lower=0, upper=100> sigma_Rmort;
    real<lower=0, upper=100> sigma_H;
    real<lower=0, upper=100> sigma_Hicu;

    real<upper=0> frac_hospmild_i45;
    real<lower=0, upper=0.1> frac_hospmild_slope;
    real<upper=0> frac_iculive_i45;
    real<lower=0, upper=0.1> frac_iculive_slope;
    real<upper=0> frac_mort_i45;
    real<lower=0, upper=0.1> frac_mort_slope;

    simplex[nage] theta[nobs_I];
    real<lower=1e-10> sigma_alpha;
}
transformed parameters {
    real<lower=0.0> x[nage,nstate,nt];
    matrix[nage,nt] beta;
    matrix[nknots_age+p_age-1,nbeta+p-1] control_beta;
    matrix[nage,nt] mu_Iobs;
    real Rt[nt];
    vector[nt] fractest = bsplines_fractest * control_fractest;

    vector[nt] imports = bsplines_imports * control_imports;

    matrix[nage,5] frac_I;

    {
        // variables in curly brackets will not have output, they are local variables
        int ibetaoffset;

        real newE;
        real newI;
        real newrec_asym;
        real newrec_mild;
        real newrec_mod;
        real newrec_icu;
        real newhosp_mod;
        real newhosp_icu;
        real newhosp_mort;
        real newmort;

        // for error detection
        real x_sum;

        vector[nage] Ipreh_mod_cur;
        vector[nage] Ipreh_icu_cur;
        vector[nage] Ipreh_mort_cur;
        vector[nage] Hicu_live_cur;
        vector[nage] Hicu_mort_cur;

        vector[nage] Ipreh_mod_new;
        vector[nage] Ipreh_icu_new;
        vector[nage] Ipreh_mort_new;
        vector[nage] Hicu_live_new;
        vector[nage] Hicu_mort_new;

        real I_cur;

        real newrec;
        int ivac = 1;

        real tmp;

        //////////////////////////////////////////
        // compute frac_I

        {
            real frac_hospmild;
            real frac_iculive;
            real frac_mort;
            real frac_mild;

            for (iage in 1:nage){
                frac_hospmild = 10^(frac_hospmild_i45 + (age[iage]-45.0) * frac_hospmild_slope);
                frac_iculive =  10^(frac_iculive_i45 +  (age[iage]-45.0) * frac_iculive_slope);
                frac_mort =     10^(frac_mort_i45 +     (age[iage]-45.0) * frac_mort_slope);
                // rest defaults to mild
                frac_mild = 1.0-(frac_asym[iage] + frac_hospmild + frac_iculive + frac_mort);
                if (frac_mild < 0.0){
                    // should only happen during initialization
                    tmp = frac_asym[iage] + frac_hospmild + frac_iculive + frac_mort;
                    frac_mild = 0.0;
                    frac_hospmild = frac_hospmild / tmp;
                    frac_iculive = frac_iculive / tmp;
                    frac_mort = frac_mort / tmp;
                }

                // 1: asym, 2: mild, 3: hosp (non-ICU), 4: non-mort ICU, 5: mort ICU
                frac_I[iage,1] = frac_asym[iage];
                frac_I[iage,2] = frac_mild;
                frac_I[iage,3] = frac_hospmild;
                frac_I[iage,4] = frac_iculive;
                frac_I[iage,5] = frac_mort;
            }
        }

        //////////////////////////////////////////
        // initial conditions
        //
        // S = 1
        // E = 2
        // Iasym = 3
        // Imild = 4
        // Ipreh = 5
        // Hmod  = 6
        // Hicu  = 7
        // Rlive = 8
        // Rmort = 9
        // Vsus = 10
        // Vrec = 11

        for (iage in 1:nage){

            // S is filled further below
            for (ifate in 2:10){
                x[iage,ifate,1] = x_ini_2to10[iage,ifate-1];
            }
            x[iage,Vrec,1] = V_ini[iage] - x[iage,Vsus,1];
            // this can occur but should not happen too often
            if (x[iage,Vrec,1] < 0){
                x[iage,Vsus,1] = V_ini[iage];
                x[iage,Vrec,1] = 0.0;
            }

            // Ipreh is the sum of Ipreh_*_cur
            tmp = frac_I[iage,3] + frac_I[iage,4] + frac_I[iage,5];
            Ipreh_mod_cur[iage]  = x[iage,Ipreh,1] * frac_I[iage,3] / tmp;
            Ipreh_icu_cur[iage]  = x[iage,Ipreh,1] * frac_I[iage,4] / tmp;
            Ipreh_mort_cur[iage] = x[iage,Ipreh,1] - Ipreh_icu_cur[iage] - Ipreh_mod_cur[iage];

            // Hicu is the sum of Hicu_*_cur
            Hicu_live_cur[iage] = x[iage,Hicu,1] * frac_I[iage,4] / (frac_I[iage,4] + frac_I[iage,5]);
            Hicu_mort_cur[iage] = x[iage,Hicu,1] - Hicu_live_cur[iage];

            x[iage,S,1] = npop[iage];
            for (ifate in 2:nstate){
                x[iage,S,1] -= x[iage,ifate,1];
            }
            if (x[iage,S,1] < 0){
                x[iage,Rlive,1] += x[iage,S,1];
                x[iage,S,1] = 0.0;
            }
        }

        for (iage in 1:nage){
            Ipreh_mod_new[iage] = Ipreh_mod_cur[iage];
            Ipreh_icu_new[iage] = Ipreh_icu_cur[iage];
            Ipreh_mort_new[iage] = Ipreh_mort_cur[iage];
        }

        for (iage in 1:nage){
            mu_Iobs[iage,1] = mu_Iobs_ini[iage];
        }

        //////////////////////////////////////////
        // test

        x_sum = 0.0;
        for (iage in 1:nage){
            x_sum = x_sum + sum(x[iage,:,1]);
        }
        //print("debug:    mu Rlive: ", mu_x_ini_2to10[:,Rlive-1])
        //print("debug: sigma Rlive: ", sigma_x_ini_2to10[:,Rlive-1])
        //print("debug:   ini Rlive: ", x_ini_2to10[:,Rlive-1])
        if (fabs(x_sum-npoptotal) > 1e-1){
            print("debug: mu_x_ini_2to10:")
            for (ifate in 1:8){
                print(mu_x_ini_2to10[:,ifate])
            }
            print("debug: sigma_x_ini_2to10:")
            for (ifate in 1:8){
                print(sigma_x_ini_2to10[:,ifate])
            }
            for (ifate in 1:nstate){
                print("ifate = ", ifate, ": sum = ", sum(x[:,ifate,1]))
            }
            for (iage in 1:nage){
                print("age group ",iage,": sum(x) = ",sum(x[iage,:,1]), ", npop = ", npop[iage])
            }
            reject("Initialization failure, net surplus: ", x_sum-npoptotal)
        }

        //////////////////////////////////////////
        // fill control_beta

        //for (iage in 1:nknots_age+p_age-1){
        for (iage in 2:nknots_age+p_age-2){
            for (ibeta in 1:nbeta_est){
                control_beta[iage,ibeta] = beta_est[iage-1,ibeta]; // added -1
            }
            if (splinemode == 2){
                control_beta[iage,nbeta_est+1] = beta_est[iage-1,nbeta_est]; // added -1
                ibetaoffset = nbeta_est + 1;
            } else {
                ibetaoffset = nbeta_est;
            }
            for (iinter in 1:ninter){
                control_beta[iage,ibetaoffset+4*iinter-3] = control_beta[iage,ibetaoffset+4*iinter-4];
                control_beta[iage,ibetaoffset+4*iinter-2] = control_beta[iage,ibetaoffset+4*iinter-4];
                control_beta[iage,ibetaoffset+4*iinter-1] = control_beta[iage,ibetaoffset+4*iinter-4];
                control_beta[iage,ibetaoffset+4*iinter] = control_beta[iage,ibetaoffset+4*iinter-1] * beta_change[iinter];
            }
            for (ibeta in nbeta:nbeta+p-1){
                control_beta[iage,ibeta] = control_beta[iage,nbeta-1];
            }
        }
        control_beta[1,:] = control_beta[2,:];
        control_beta[nknots_age+p_age-1,:] = control_beta[nknots_age+p_age-2,:];

        beta = bsplines_age * control_beta * transpose(bsplines);
        //////////////////////////////////////////
        // the SEIR model

        for (it in 1:nt-1){

            Rt[it] = 0.0;
            I_cur = sum(x[:,Iasym,it]) + sum(x[:,Imild,it]) + sum(x[:,Ipreh,it]) + imports[it];

            //////////////////////////////////////////
            // vaccination immunization effect

            for (iage in 1:nage){
                x[iage,Vsus,it+1] = x[iage,Vsus,it];
                x[iage,Vrec,it+1] = x[iage,Vrec,it];
            }
            if (ivac <= nvac && it == tvac[ivac]){
                for (iage in 1:nage){
                    // fraction S/(S+Rlive)
                    newrec = (x[iage,S,it] / (x[iage,S,it] + x[iage,Rlive,it])) * vac[iage,ivac];
                    x[iage,Vsus,it+1] = x[iage,Vsus,it+1] + newrec;
                    x[iage,S,it] = x[iage,S,it] - newrec;
                    newrec = (vac[iage,ivac] - newrec);
                    x[iage,Vrec,it+1] = x[iage,Vrec,it+1] + newrec;
                    x[iage,Rlive,it] = x[iage,Rlive,it] - newrec;
                }
                ivac = ivac + 1;
            }

            for (iage in 1:nage){

                //////////////////////////////////////////
                // set transition variables
                newrec_asym = 1.0/duration_rec_asym * x[iage,Iasym,it];
                newrec_mild = 1.0/duration_rec_mild * x[iage,Imild,it];
                newrec_mod = 1.0/duration_hosp_mod * x[iage,Hmod,it];
                newrec_icu = 1.0/duration_hosp_icu * Hicu_live_cur[iage];
                newmort = 1.0/duration_hosp_icu * Hicu_mort_cur[iage];

                newhosp_mod = 1.0/duration_pre_hosp * Ipreh_mod_cur[iage];
                newhosp_icu = 1.0/duration_pre_hosp * Ipreh_icu_cur[iage];
                newhosp_mort = 1.0/duration_pre_hosp * Ipreh_mort_cur[iage];

                //Ipreh_mod_new[iage] = Ipreh_mod_cur[iage];    // step not necessary, they are already equal
                //Ipreh_icu_new[iage] = Ipreh_icu_cur[iage];    // step not necessary, they are already equal
                //Ipreh_mort_new[iage] = Ipreh_mort_cur[iage];  // step not necessary, they are already equal

                //////////////////////////////////////////
                // S -> E -> I
                newE = fmin(x[iage,S,it], beta[iage,it]/npop[iage] * I_cur * x[iage,S,it]);
                newI = 1/duration_lat * x[iage,E,it]; // duration_lat > 1

                x[iage,S,it+1] = x[iage,S,it] - newE;
                x[iage,E,it+1] = x[iage,E,it] + newE - newI;

                // frac_I index: 1: asym, 2: mild, 3: hosp (non-ICU), 4: non-mort ICU, 5: mort ICU
                x[iage,Iasym,it+1] = x[iage,Iasym,it] + frac_I[iage,1] * newI;
                x[iage,Imild,it+1] = x[iage,Imild,it] + frac_I[iage,2] * newI;

                Ipreh_mod_new[iage]  = Ipreh_mod_new[iage]  + frac_I[iage,3] * newI;
                Ipreh_icu_new[iage]  = Ipreh_icu_new[iage]  + frac_I[iage,4] * newI;
                Ipreh_mort_new[iage] = Ipreh_mort_new[iage] + frac_I[iage,5] * newI;

                // observed infected
                mu_Iobs[iage,it+1] = mu_Iobs[iage,it] + newI * fractest[it+1];

                // compute Rt
                // NOTE: requires age structure data for S and time-dependent beta
                // and is hence not easy to calculate in generated quantities block
                Rt[it] = Rt[it] + x[iage,S,it]/npop[iage] * beta[iage,it] * (frac_I[iage,1]*duration_rec_asym + frac_I[iage,2]*duration_rec_mild + (1.0-frac_I[iage,1]-frac_I[iage,2])*duration_pre_hosp);

                //////////////////////////////////////////
                // I -> recovery or hospital

                x[iage,Iasym,it+1] = x[iage,Iasym,it+1] - newrec_asym;
                x[iage,Imild,it+1] = x[iage,Imild,it+1] - newrec_mild;
                Ipreh_mod_new[iage]  = Ipreh_mod_new[iage]  - newhosp_mod;
                Ipreh_icu_new[iage]  = Ipreh_icu_new[iage]  - newhosp_icu;
                Ipreh_mort_new[iage] = Ipreh_mort_new[iage] - newhosp_mort;

                x[iage,Ipreh,it+1] = Ipreh_mod_new[iage] + Ipreh_icu_new[iage] + Ipreh_mort_new[iage];

                //////////////////////////////////////////
                // in hospital

                x[iage,Hmod,it+1]  = x[iage,Hmod,it]    + newhosp_mod  - newrec_mod;
                Hicu_live_new[iage] = Hicu_live_cur[iage] + newhosp_icu  - newrec_icu;
                Hicu_mort_new[iage] = Hicu_mort_cur[iage] + newhosp_mort - newmort;

                x[iage,Hicu,it+1] = Hicu_live_new[iage] + Hicu_mort_new[iage];

                //////////////////////////////////////////
                // recovery/death

                x[iage,Rlive,it+1] = x[iage,Rlive,it] + newrec_asym + newrec_mild + newrec_mod + newrec_icu;
                x[iage,Rmort,it+1] = x[iage,Rmort,it] + newmort;

                //////////////////////////////////////////
                // set cur to new

                Ipreh_mod_cur[iage] = Ipreh_mod_new[iage];
                Ipreh_icu_cur[iage] = Ipreh_icu_new[iage];
                Ipreh_mort_cur[iage] = Ipreh_mort_new[iage];
                Hicu_live_cur[iage] = Hicu_live_new[iage];
                Hicu_mort_cur[iage] = Hicu_mort_new[iage];
            }

            //////////////////////////////////////////
            // test

            x_sum = 0.0;
            for (iage in 1:nage){
                x_sum = x_sum + sum(x[iage,:,it+1]);
            }
            if (fabs(x_sum-npoptotal) > 1e-1){
                for(ifate in 1:nstate){
                    print("ifate = ", ifate, ": sum = ", sum(x[:,ifate,it+1]))
                }
                reject("Model is leaking, net gain: ", x_sum-npoptotal, " (it = ", it, ")")
            }
        }
        // compute Rt for last timestep
        Rt[nt] = 0.0;
        for (iage in 1:nage){
            Rt[nt] = Rt[nt] + x[iage,S,nt]/npop[iage] * beta[iage,nt] * (frac_I[iage,1]*duration_rec_asym + frac_I[iage,2]*duration_rec_mild + (1.0-frac_I[iage,1]-frac_I[iage,2])*duration_pre_hosp);
        }
    }
}
model {
    real tmp;
    //////////////////////////////////////////
    // prior distributions

    for (ifate in 1:1){
        x_ini_2to10[:,ifate] ~ lognormal(mu_x_ini_2to10[:,ifate], sigma_x_ini_2to10[:,ifate]);
    }
    for (ifate in 2:9){
        for (iage in 1:nage){
            x_ini_2to10[iage,ifate] ~ normal(mu_x_ini_2to10[iage,ifate], sigma_x_ini_2to10[iage,ifate]) T[0,];
        }
    }

    frac_hospmild_i45 ~ normal(mu_frac_hospmild_i45, sigma_frac_hospmild_i45);
    frac_iculive_i45 ~ normal(mu_frac_iculive_i45, sigma_frac_iculive_i45);
    frac_mort_i45 ~ normal(mu_frac_mort_i45, sigma_frac_mort_i45);

    frac_hospmild_slope ~ normal(mu_frac_hospmild_slope, sigma_frac_hospmild_slope);
    frac_iculive_slope ~ normal(mu_frac_iculive_slope, sigma_frac_iculive_slope);
    frac_mort_slope ~ normal(mu_frac_mort_slope, sigma_frac_mort_slope);

    sigma_beta ~ exponential(30.0);
    beta_est[:,1] ~ normal(mu_beta1, sigma_beta1);
    for (ibeta in 2:nbeta_est){
        beta_est[:,ibeta] ~ normal(beta_est[:,ibeta-1], sigma_beta);
    }
    beta_change ~ normal(mu_beta_inter, sigma_beta_inter);

    duration_lat ~ normal(mu_duration_lat, sigma_duration_lat);
    duration_rec_asym ~ normal(mu_duration_rec_asym, sigma_duration_rec_asym);
    duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);
    duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
    duration_hosp_mod ~ normal(mu_duration_hosp_mod, sigma_duration_hosp_mod);
    duration_hosp_icu ~ normal(mu_duration_hosp_icu, sigma_duration_hosp_icu);

    control_fractest ~ normal(mu_fractest, sigma_fractest);

    control_imports ~ normal(mu_imports, sigma_imports);

    //////////////////////////////////////////
    // fitting observations

    sigma_Iobs ~ exponential(lambda_Iobs);
    sigma_Rmort ~ exponential(lambda_Rmort);

    sigma_H ~ exponential(lambda_H);
    sigma_Hicu ~ exponential(lambda_Hicu);

    sigma_alpha ~ exponential(0.03);
    if (prior_only[1] == 0){
        vector[nage] alpha;
        // fit sum
        for (iobs in 1:nobs_I){
            tmp = (sum(obs_I[:,iobs])-sum(mu_Iobs[:,tobs_I[iobs]]))/sigma_Iobs;
            tmp ~ std_normal();
        }
        // fit age distribution
        for (iobs in 1:nobs_I){
            alpha = mu_Iobs[:,tobs_I[iobs]]/sum(mu_Iobs[:,tobs_I[iobs]]) * sigma_alpha + 1;
            theta[iobs] ~ dirichlet(alpha);
            obs_I[:,iobs] ~ multinomial(theta[iobs]);
        }
    }
    if (prior_only[2] == 0){
        // fit sum
        for (iobs in 1:nobs_Rmort){
            tmp = (obs_Rmort[iobs]-sum(x[:,Rmort,tobs_Rmort[iobs]]))/sigma_Rmort;
            tmp ~ std_normal();
        }
    }
    if (prior_only[3] == 0){
        // fit sum
        for (iobs in 1:nobs_H){
            tmp = (obs_H[iobs]-(sum(x[:,Hmod,tobs_H[iobs]])+sum(x[:,Hicu,tobs_H[iobs]])))/sigma_H;
            tmp ~ std_normal();
        }
    }
    if (prior_only[4] == 0){
        // fit sum
        for (iobs in 1:nobs_H){
            tmp = (obs_Hicu[iobs]-sum(x[:,Hicu,tobs_H[iobs]]))/sigma_Hicu;
            tmp ~ std_normal();
        }
    }
}
generated quantities{
    real fractested[nobs_I];
    real hospitalized[nt];
    real obs_I_sim[nage,nt];
    for (iobs in 1:nobs_I){
        // dividing here by the cumulative number as well
        // mu_Iobs[tobs_I[iobs]]/fractest[tobs_I[iobs]] is cumulative number of infectious
        // resulting quantity is the cumulative number of observed infected divided by the cumulative number of infectious in the model
        fractested[iobs] = sum(obs_I[:,iobs])/(sum(mu_Iobs[:,tobs_I[iobs]])/fractest[tobs_I[iobs]]);
    }
    for (it in 1:nt){
        hospitalized[it] = normal_rng(sum(x[:,Hmod,it]) + sum(x[:,Hicu,it]), sigma_H);
        for (iage in 1:nage){
            obs_I_sim[iage,it] = normal_rng(mu_Iobs[iage,it], sigma_Iobs);
        }
    }
}
