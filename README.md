# The Santa Cruz County COVID-19 model

The Santa Cruz County (SCZ) COVID-19 model is a time-discrete, stochastic SEIR model that uses Bayesian statistical methods, such as Hamiltonian Markov Chain Monte Carlo (MCMC) simulations, to forecast the COVID-19 pandemic in Santa Cruz County, California. The model requires a set of parameters, equations, and local data to help inform its simulations. The model is set to run 4,000 simulations and fine-tune the inputted parameters using the local data (confirmed COVID-19 hospitalizations, confirmed COVID-19 cases, and deaths). The model projects a range of different scenarios that fit the inputs provided and are displayed in the exported plots.

The model contains 9 compartments to divide COVID-19 cases into the asymptomatic, mild, and moderate to severe illness, which better informs hospitalization and death projections (see diagram below). 

![model diagram](plots/seir_diagram.png)

The Jupyter notebook using data from Santa Cruz County can be found [here](seir_santa_cruz.ipynb), a template notebook can be found [here](seir_template.ipynb).

**Note:** If you have issues, questions or find a bug please create an issue in GitHub (above).

## Additional Assumptions of the SCZ COVID-19 Model
- The model's contact rate adjusts every 5 days using spline interpolation.
- COVID-19 cases who recover gain short-term immunity.
- COVID-19 cases can be infectious 2 to 3 days prior to symptom onset.
- COVID-19 hospitalization, ICU, and death rates are calculated based on the overall age demographics of Santa Cruz County.
- Not everyone who tests positive for COVID-19 goes to the hospital.
- COVID-19 cases only die within the hospital.
- Hospitalized COVID-19 patients have a shorter duration of their infectious period because they are less likely to expose others. However, they likely will shed live virus longer, especially if immuno-compromised.
- The model does not account for spatial or network patterns.

## Model updates

### 2020-08-09:
 - Added temporary disclaimer to the plots stating that due to a problem with the State of California’s CalREDIE reporting system, cases have been underreported and projected results likely represent underestimates.
 - This problem is affecting the data being fed into the recent projections, not the model code.

### 2020-07-21:
 - Priors for model misfit parameters (`lambda_Iobs`, `lambda_Hmod`, `lambda_Hicu`, `lambda_Rmort`) can now be set by the user.
 - Santa Cruz and default parametrization now use a tighter fit to mortality data and a looser fit to case count.
 - Mean of the prior for rate of mortality has been changed from 1% to 0.5% based on observed mortality in Santa Cruz County.  

## Instruction for running the model

The core of the model is written in [Stan](https://mc-stan.org/) and a [Jupyter notebook](https://jupyter.org/) is used for reading in data, running the model and visualizing its output. Installation instructions are provided [here](installation_instructions.md). 

In order to provide a new dataset or change the model parameters, edit the begining of the notebook where each input parameter and the required data files are described in detail.

## Parameters of COVID-19 model for Santa Cruz County

All parameters listed below can be adjusted by the user.

Parameter| Value and Distribution | Literature
---------|------------------------|-----------|
Initial Contact Rate | normal(0.2, 0.02)<sup>[\[1\]](#betafootnote)</sup> | N/A|
Latent Period | normal(5, 2) | [Cascella et al](https://www.ncbi.nlm.nih.gov/books/NBK554776/); [Li et al](https://www.ncbi.nlm.nih.gov/books/NBK554776/) |
Asymptomatic Infections Days to Recover | normal(7, 5) | [He et al](https://www.nature.com/articles/s41591-020-0869-5)|
Mild Infections Days to Recover | normal(7,4)| [He et al](https://www.nature.com/articles/s41591-020-0869-5)|
Days to Hospital| normal(5,1)| [Ferguson et al](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf)|
Non-ICU Cases Days in Hospital |normal(7,1)|[Ferguson et al](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf)|
ICU Cases Days in Hospital | normal(16, 1)|[Ferguson et al](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf)|
Fraction Tested | time-dependent<sup>[\[2\]](#fractionfootnote)</sup> |N/A|
Fraction of moderate cases (Hospitalized but non-ICU) | Dirichlet with mean 0.07<sup>[\[3\]](#dirichletfootnote)</sup> | [Stanford model](https://surf.stanford.edu/covid-19-tools/covid-19/); [Ferguson et al](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf); [Verity et al](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext#seccestitle200)
Fraction of severe cases (ICU and alive) | Dirichlet with mean 0.02<sup>[\[3\]](#dirichletfootnote)</sup> | [Stanford model](https://surf.stanford.edu/covid-19-tools/covid-19/); [Ferguson et al](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf); [Verity et al](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext#seccestitle200) |
Fraction of severe cases (ICU and dead) | Dirichlet with mean 0.005<sup>[\[3\]](#dirichletfootnote)</sup> |N/A|
Fraction Asymptomatic cases | Dirichlet with mean 0.178<sup>[\[3\]](#dirichletfootnote)</sup> |[Nishiura et al](https://www.ncbi.nlm.nih.gov/pubmed/32145466); [Mizumoto et al](https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.10.2000180#r13)|
Fraction of mild cases (non-hospitalized) | Dirichlet with mean 1-sum(0.07, 0.02, 0.005, 0.178)<sup>[\[3\]](#dirichletfootnote)</sup> | N/A|
Population of Santa Cruz County | 273, 213| https://www.census.gov/quickfacts/santacruzcountycalifornia |

<a name="betafootnote">[1]</a> A time-dependent contact rate is then estimated from the data using an AR(1) process and spline interpolation. 

<a name="fractionfootnote">[2]</a> Currently increased over time using key dates, see [notebook](seir_santa_cruz.ipynb) for details. 

<a name="dirichletfootnote">[3]</a> All 5 fractions are drawn from the same Dirichlet distribution. 

## Contributors

 * Jann Paul Mattern [jpmattern](https://github.com/jpmattern)
 * Mikala Caton [mtcaton](https://github.com/mtcaton)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


