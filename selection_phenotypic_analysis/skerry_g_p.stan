//
// Stan function for estimation of selection in the skerry data
// author Roger Butlin

// The input data 
data {
  int g;  // number ofgeneration per year
  int<lower=0> Ng;   // number of years of the experiment
  int<lower=0> Nd;   // number of snails - donor population
  int<lower=0> Ns96;  // number of snails - 1996 sample
  int<lower=0> Ns02;  // number of snails - 2002 sample
  int<lower=0> Ns05;  // number of snails - 2005 sample
  int<lower=0> Ns18;  // number of snails - 2018 sample
  int<lower=0> Ns21;  // number of snails - 2021 sample
  int<lower=0> Nw;  // number of snails - Wave reference
  vector[Nd] don;  // donor phenotypes
  vector[Ns96] S96;  // skerry phenotypes 1996
  vector[Ns02] S02;  // skerry phenotypes 2002
  vector[Ns05] S05;  // skerry phenotypes 2005
  vector[Ns18] S18;  // skerry phenotypes 2018
  vector[Ns21] S21;  // skerry phenotypes 2021
  vector[Nw] wav;   // wave phenotypes
  real p_est; // estimated mean plastic effect, p
  real sdp;   // among transect sd of p
  real Vg_est; // estimated mean Vg
  real sd_Vg;  // among transect sd of Vg
}

// The parameters accepted by the model. 
parameters {
  real p;              // plasticity
  real<lower=0> Vg;    // genetic variance
  real<lower=0, upper=50> Vs_Vp;    // strength of selection, relative to phenotypic variance [flat prior 0-50]
  real crab;         // mean phenotype in the Crab donor population
  real wave_opt;     // mean phenotype in the Wave reference population
  real<lower=0> Ve; // non-additive and environmental variance, assumed equal for donor/skerry/wave
}

// transformed parameters, including marginal effect of generations per year
transformed parameters {
  real lp;  // log probability of the data, given the model
  real<lower=0> sd_phen; // phenotypic standard deviation
  real<lower=0> Vs;      // strength of selection
  vector [g*Ng+1] phen;  // mean phenotype in each generation
  sd_phen = sqrt(Vg+Ve);
  Vs = Vs_Vp*sd_phen^2;
  
    phen[1] = crab + p;
    lp = 0;
    for (y in 2:(g*Ng+1)){
      phen[y] = phen[y-1] + (wave_opt - phen[y-1]) * Vg/(sd_phen^2 + Vs); // expected phenotype in generation y, Charlesworth and Charlesworth (2010) p186
    }
    for (j in 1:Nd){lp = lp + normal_lpdf(don[j] | crab, sd_phen);}
    for (j in 1:Ns96){lp = lp + normal_lpdf(S96[j] | phen[g*4+1], sd_phen);}
    for (j in 1:Ns02){lp = lp + normal_lpdf(S02[j] | phen[g*10+1], sd_phen);}
    for (j in 1:Ns05){lp = lp + normal_lpdf(S05[j] | phen[g*13+1], sd_phen);}
    for (j in 1:Ns18){lp = lp + normal_lpdf(S18[j] | phen[g*26+1], sd_phen);}
    for (j in 1:Ns21){lp = lp + normal_lpdf(S21[j] | phen[g*29+1], sd_phen);}
    for (j in 1:Nw){lp = lp + normal_lpdf(wav[j] | wave_opt, sd_phen);}

}

// The model to be estimated. 
model {
  p ~ normal(p_est,sdp);      // prior on p
  Vg ~ normal(Vg_est,sd_Vg);  // prior on Vg
  target += lp;  
}


generated quantities {
  real w_crab;  // fitness of the Crab phenotype in the Wave environment
  vector[g*Ng+1] fit_phen; // fitted mean phenotype for each generation
  w_crab = exp(0-(wave_opt - crab)^2/(2*Vs));
  for (k in 1:(g*Ng+1)){fit_phen[k] = phen[k];}
}

