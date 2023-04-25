//-----------------------------------------------------------------------------//
// yukonChinookRunRecon.cpp                                                    //
// Simulation functions for Yukon River Chinook run reconstruction             //
//                                                                             //
// Copyright 2023 by Landmark Fisheries Research, Ltd.                         //
//                                                                             //
// This software is provided to DFO in the hope that it will be                //
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                        //
//                                                                             //
// ALL INTELLECTUAL PROPERTY REMAINS WITH LANDMARK FISHERIES RESEARCH, LTD.    //
// THIS SOFTWARE MAY NOT BE REDISTRIBUTED, SUBLICENCED, COPIED, OR SHARED      //
// OUTSIDE OF ESSA TECHNOLOGIES WITHOUT THE EXPRESS WRITTEN CONSENT OF         //
// LANDMARK FISHERIES RESEARCH, LTD.                                           //
//                                                                             //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" //
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   //
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  //
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE    //
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         //
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        //
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    //
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     //
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     //
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  //
// POSSIBILITY OF SUCH DAMAGE.                                                 //
//-----------------------------------------------------------------------------//


#include <TMB.hpp>
#include <iostream>

// isNA() - replicates is.na() from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// square()
template <class Type> 
Type square(Type x){
  return pow(x,2);
}

// Bounded inverse logit
template <class Type> 
Type binvlogit(Type x)
{
  Type lb = -1;
  Type ub = 1;
  return invlogit(x)*(ub-lb) + lb;
}

// posfun() - replicates Dave Fournier's ADMB function
// Compares a quantity x to a threshold epsilon and if below the threshold, 
// increases x above the threshold and increments a penalty
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0.));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2.0)-x/eps));
}

// Objective function ------------------------------------------------------ //
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Access multivariate distributions
  using namespace density;

  // Data
  DATA_ARRAY(obsC_rht)        // Observed catch
  DATA_ARRAY(catchSD_rh);     // Observed catch SD
  DATA_ARRAY(xB_jgt);         // Border age/sex-composition
  DATA_ARRAY(xC_arht);        // Catch age/sex-composition
  DATA_ARRAY(n_pdtg);         // Daily border counts by stock
  DATA_ARRAY(E_dtg);          // Total daily counts
  DATA_VECTOR(I_t);           // Border passage indices from mark-recapture
  DATA_VECTOR(mrCV_t);         // SEs for border passage indices from mark-recapture

  // Constants
  DATA_VECTOR(day_d);         // Model days (ordinal)
  DATA_INTEGER(aMin);         // Minimum age class
  DATA_INTEGER(aMax);         // Maximum age class
  DATA_ARRAY(z_ast);          // Reproductive quality
  DATA_ARRAY(RLM_astm);       // Ratio of fish length to mesh perimeter
  DATA_IARRAY(mesh_rht);      // Mesh size index
  DATA_ARRAY(v_asg);          // Fish wheel selectivity
  DATA_SCALAR(weightI);       // Weight on border passage index NLL
  DATA_SCALAR(etaSD);         // SEs for border passage indices from mark-recapture
  DATA_ARRAY(pearPrior);      // Pearson selectivity prior
  DATA_VECTOR(phiPrior);      // Autocorrelation prior
  DATA_VECTOR(UMSYPrior);     // UMSY prior

  // Model dimensions
  int nP = n_pdtg.dim(0);     // Number of stocks
  int nD = n_pdtg.dim(1);     // Number of days
  int nT = n_pdtg.dim(2);     // Number of calendar years
  int nG = n_pdtg.dim(3);     // Number of survey types (sonar, fish wheel)
  int nR = obsC_rht.dim(0);   // Number of regions (US, CA)
  int nH = obsC_rht.dim(1);   // Number of fishery types (subsistence, commercial)
  int nA = aMax-aMin+1;       // Number of ages
  int nY = nT + nA - 1;       // Number of brood years
  int nM = RLM_astm.dim(3);   // Number of mesh sizes

  // Spawner-recruit parameters
  PARAMETER_ARRAY(lnR_py);       // Realized recruitment
  PARAMETER_ARRAY(lnf_rht);      // Fully-selected fishing mortality
  PARAMETER_VECTOR(lnUMSY_p);    // UMSY
  PARAMETER_VECTOR(lnSMSY_p);    // SMSY
  PARAMETER_VECTOR(lnRecSD_p);   // Recruitment deviation SD
  PARAMETER_VECTOR(logitPhi_p);  // Lag-1 autocorrelation in recruitment
  PARAMETER_VECTOR(pearPars);    // Pearson gillnet selectivity parameters
  PARAMETER_VECTOR(logitpFem_y); // Proportion female
  PARAMETER(lnpFemSD);           // Proportion female SD
  PARAMETER_ARRAY(initEta_asy);  // Raw proportions-at-age

  // Run reconstruction parameters
  PARAMETER_VECTOR(lnArrivMu_p);  // Mean date of arrival in initial year (log scale)
  PARAMETER_VECTOR(lnArrivSD_p);  // SD around mean Julian date of arrival (log scale)
  PARAMETER_ARRAY(arrivErr_pt);   // Process error in mean date of arrival
  PARAMETER_VECTOR(lnErrSD_p);    // Process error standard deviation (log scale)
  PARAMETER_ARRAY(logitCor_pp);   // Correlation matrix for process errors (log scale)
  PARAMETER_ARRAY(lnqE_pg);       // Count catchability (log scale)
  PARAMETER_VECTOR(lnqI_p);       // Index catchability
  PARAMETER_ARRAY(lnDisp_tg);     // Negative binomial dispersion

  // Transformed parameters
  vector<Type> arrivSD_p = exp(lnArrivSD_p);
  vector<Type> errSD_p = exp(lnErrSD_p);
  vector<Type> UMSY_p = exp(lnUMSY_p);
  vector<Type> SMSY_p = exp(lnSMSY_p);
  vector<Type> recSD_p = exp(lnRecSD_p);
  vector<Type> phi_p(nP);
  Type pearLambda = pearPars(0);
  Type pearTheta = exp(pearPars(1));
  Type pearSigma = exp(pearPars(2));
  Type pearTau = pearPars(3);
  vector<Type> alpha_p = exp(UMSY_p) / (1-UMSY_p); // Ricker productivity
  vector<Type> beta_p = UMSY_p / SMSY_p;           // Ricker density dependence
  vector<Type> R0_p = log(alpha_p) / beta_p;       // Unfished recruitment

  // Latent variables
  array<Type> N_asrpt(nA,2,nR,nP,nT);   // Abundance
  array<Type> v_astm(nA,2,nT,nM);       // Selectivity
  vector<Type> borderPass_t(nT);        // Total border passage
  array<Type> eta_aspy(nA,2,nP,nY);     // Age proportions
  array<Type> mu_pt(nP,nT);             // Mean Julian date of arrival
  array<Type> omega_py(nP,nY);          // Mean Julian date of arrival
  array<Type> J_dptg(nD,nP,nT,nG);      // Daily observable abundance by stock
  array<Type> rho_dpt(nD,nP,nT);        // Proportion of stock arriving each day
  matrix<Type> cor_pp(nP,nP);           // Arrival timing correlation matrix
  array<Type> F_arhst(nA,nR,nH,2,nT);   // Fishery-specific mortality
  array<Type> F_arst(nA,nR,2,nT);       // Total mortality
  array<Type> C_asrht(nA,2,nR,nH,nT);   // Catch
  array<Type> xB_asgt(nA,2,nG,nT);      // Age/sex-composition at border
  array<Type> xhatC_arht(nA,nR,nH,nT);  // Age/sex-composition in catch
  array<Type> S_aspt(nA,2,nP,nT);       // Spawning escapement
  array<Type> Rexp_py(nP,nY);           // Expected recruitment  
  array<Type> R_py(nP,nY);              // Realized recruitment
  array<Type> R_spy(2,nP,nY);           // Realized recruitment by sex
  array<Type> Z_pt(nP,nT);              // Total reproductive output
  array<Type> Ihat_dtg(nD,nT,nG);       // Predicted abundance
  vector<Type> mrIhat_t(nT);            // Predicted border passage
  array<Type> Phat_pdtg(nP,nD,nT,nG);   // Predicted proportions by stock
  array<Type> runSize_pt(nP,nT);        // Predicted proportions by stock
  vector<Type> runSize_t(nT);           // Predicted proportions by stock
  vector<Type> pFem_y(nY);              // Predicted proportions by stock
  vector<Type> relReproOutput_t(nT);    // Relative reproductive output
  vector<Type> meanMu_p(nP);            // Mean arrival timing

  // Likelihood components
  array<Type> nllI_tg(nT,nG);       // Total daily count NLL
  vector<Type> nllP_g(nG);          // Stock-composition NLL
  array<Type> nllC_rht(nR,nH,nT);   // Catch NLL
  Type nllxB  = 0;                  // Border age/sex-composition NLL
  array<Type> nllxC_rht(nR,nH,nT);  // Catch age/sex-composition NLL
  Type nllMR  = 0;                  // Mark-capture border passage NLL
  array<Type> nllR_py(nP,nY);       // Recruitment NLL
  Type nlp    = 0;                  // Prior penalty
  Type varPen = 0;                  // posfun penalty

  // Initialize variables
  R_py.fill(0);
  nllI_tg.fill(0);
  nllP_g.fill(0);
  F_arst.fill(0);
  C_asrht.fill(0);
  S_aspt.fill(0);
  Z_pt.fill(0);
  mrIhat_t.fill(0);
  eta_aspy.fill(1);
  xB_asgt.fill(0);
  xhatC_arht.fill(1e-6);
  nllC_rht.fill(0);
  nllxC_rht.fill(0);
  nllR_py.fill(0);
  runSize_t.fill(0);
  borderPass_t.fill(0);
  meanMu_p.fill(0);

  for( int p=0; p<nP; p++ )
    phi_p(p) = binvlogit(logitPhi_p(p));

  // Intialize random walk
  mu_pt.col(0) = exp(lnArrivMu_p);

  for( int y=0; y<nY; y++ )
  {
    pFem_y(y) = invlogit(logitpFem_y(y));
    for( int p=0; p<nP; p++ )
    {
      R_spy(0,p,y) = exp(lnR_py(p,y)) * (1 - pFem_y(y));
      R_spy(1,p,y) = exp(lnR_py(p,y)) * pFem_y(y);

      for( int s=0; s<2; s++ )
      {
        eta_aspy.col(y).col(p).col(s).segment(0,nA-1) = exp(initEta_asy.col(y).col(s));
        eta_aspy.col(y).col(p).col(s) /= eta_aspy.col(y).col(p).col(s).sum();

        if( y>0 )
        {
          vector<Type> tmpeta1 = initEta_asy.col(y).col(s);
          vector<Type> tmpeta0 = initEta_asy.col(y-1).col(s);
          nlp -= dnorm( tmpeta1, tmpeta0, etaSD, TRUE ).sum();
        }

      }
    }
  }

  // In-season abundance dynamics
  for( int t=0; t<nT; t++ )
  {
    // Fishery selectivity by mesh size
    for( int m=0; m<nM; m++ )
    {
      vector<Type> tmp6(2);
      vector<Type> tmp7(nA);
      for( int s=0; s<2; s++ )
      {
        // Pearson selectivity
        Type tmp1 = pow(1 + square(pearLambda)/(4 * square(pearTheta)),pearTheta);
        vector<Type> tmp2 = RLM_astm.col(m).col(t).col(s) - (pearSigma * pearLambda)/(2 * pearTheta) - pearTau;
        vector<Type> tmp3 = 1 + (tmp2*tmp2)/pow(pearSigma,2);
        vector<Type> tmp4(nA);
        for( int a=0; a<nA; a++ )
          tmp4(a)= pow(tmp3(a),-pearTheta);
        vector<Type> tmp5 = exp(-pearLambda * (atan(tmp2/pearSigma) + atan(pearLambda/(2 * pearTheta))));
        tmp7 = tmp1*tmp4*tmp5;
        v_astm.col(m).col(t).col(s) = tmp7;
        tmp6(s) = max(tmp7);
      }  // next s
      matrix<Type> vtmp = v_astm.col(m).col(t);
      v_astm.col(m).col(t) = vtmp / max(tmp6);
    }  // next m

    // Fishing mortality
    for( int r=0; r<nR; r++ )
    {    
      for( int s=0; s<2; s++ )
      {
        for( int h=0; h<nH; h++ )
        {
          F_arhst.col(t).col(s).col(h).col(r) = v_astm.col(mesh_rht(r,h,t)).col(t).col(s)*exp(lnf_rht(r,h,t));

          for( int a=0; a<nA; a++ )
            F_arst(a,r,s,t) += F_arhst(a,r,h,s,t);
        
        }  // next h
      }  // next s
    }  // next r

    for( int p=0; p<nP; p++ )
    {
      
      for( int s=0; s<2; s++ )
      {      
        // Arrival to AK
        for( int a=0; a<nA; a++ )
          N_asrpt(a,s,0,p,t) = R_spy(s,p,t+nA-a-1) * eta_aspy(a,s,p,t+nA-a-1);

        // Arrival to pilot and border
        for( int r=1; r<nR; r++ )
        {
          vector<Type> tmpF0_a = F_arst.col(t).col(s).col(r-1);
          N_asrpt.col(t).col(p).col(r).col(s) = N_asrpt.col(t).col(p).col(r-1).col(s) * exp(-tmpF0_a);
        }

        // Mark-recapture predicted index
        mrIhat_t(t) += exp(lnqI_p(p))*N_asrpt.col(t).col(p).col(nR-1).col(s).sum();

        // Escapement
        vector<Type> tmpF_a = F_arst.col(t).col(s).col(nR-1);
        S_aspt.col(t).col(p).col(s) = N_asrpt.col(t).col(p).col(nR-1).col(s) * exp(-tmpF_a);

        // Total reproductive output
        Z_pt(p,t) += (S_aspt.col(t).col(p).col(s) *
                      z_ast.col(t).col(s)).sum();

        // Catch by region and gear
        for( int r=0; r<nR; r++ )
        {
          tmpF_a = F_arst.col(t).col(s).col(r);
          for( int h=0; h<nH; h++ )
          {
            C_asrht.col(t).col(h).col(r).col(s) += N_asrpt.col(t).col(p).col(r).col(s) *
                                                   F_arhst.col(t).col(s).col(h).col(r) *
                                                   (1-exp(-tmpF_a)) /
                                                   tmpF_a;
          }  // next h
        }  // next r
      }  // next s

      // Border age/sex-comps
      for( int g=0; g<nG; g++ )
        xB_asgt.col(t).col(g) += N_asrpt.col(t).col(p).col(2)*v_asg.col(g);

      // Save quantities of interest for ADREPORTing
      runSize_pt(p,t) = N_asrpt.col(t).col(p).col(0).sum();
      runSize_t(t) += runSize_pt(p,t);
      borderPass_t(t) += N_asrpt.col(t).col(p).col(nR-1).sum();

    }  // next p

    relReproOutput_t(t) = Z_pt.col(t).sum() / S_aspt.col(t).sum();
  
    // Border age/sex-comps NLL
    for( int g=0; g<nG; g++ )
    {
      xB_asgt.col(t).col(g) /= xB_asgt.col(t).col(g).sum();
      if( !isNA(xB_jgt(0,g,t)) )
      {
        vector<Type> obs_j = xB_jgt.col(t).col(g);
        vector<Type> est_j(2*nA);
        est_j.segment(0,nA) = xB_asgt.col(t).col(g).col(0);
        est_j.segment(nA,nA) = xB_asgt.col(t).col(g).col(1);
        est_j /= est_j.sum();
        nllxB -= dmultinom( obs_j, est_j, TRUE );
      }
    }

    // Catch
    for( int h=0; h<nH; h++ )
    {
      for( int r=0; r<nR; r++ )
      {
        // Total catch across ages and sexes
        Type Csum = C_asrht.col(t).col(h).col(r).sum();

        // Predictd catch age/sex composition
        for( int s=0; s<2; s++ )
          xhatC_arht.col(t).col(h).col(r) += C_asrht.col(t).col(h).col(r).col(s);

        xhatC_arht.col(t).col(h).col(r) /= Csum;

        // Catch likelihood
        if( obsC_rht(r,h,t)>0 )
        {
          nllC_rht(r,h,t) -= dnorm( log(obsC_rht(r,h,t)+1e-10), log(Csum+1e-10), catchSD_rh(r,h), TRUE );

          // Catch age/sex composition likelihood
          if( !isNA(xC_arht(0,r,h,t)) )
          {
            vector<Type> obs_j = xC_arht.col(t).col(h).col(r);
            vector<Type> est_j = xhatC_arht.col(t).col(h).col(r);
            nllxC_rht(r,h,t) -= dmultinom( obs_j, est_j, TRUE );
          }
        }
      }  // next r
    }  // next h

    // RUN RECONSTRUCTION

    for( int p=0; p<nP; p++ )
    {

      for( int p2=0; p2<nP; p2++ )
        cor_pp(p,p2) = binvlogit(logitCor_pp(p,p2));

      // Random walk in arrival timing
      if( t>0 )
        mu_pt(p,t) = mu_pt(p,t-1)*exp(arrivErr_pt(p,t-1));

      // Relative daily arrivals
      for( int d=0; d<nD; d++ )
        rho_dpt(d,p,t) = exp( -square(day_d(d)-mu_pt(p,t))/(2*square(arrivSD_p(p))) );
      
      // Convert relative daily arrivals to proportions
      rho_dpt.col(t).col(p) /= rho_dpt.col(t).col(p).sum();

      // Store predicted abundance
      for( int g=0; g<nG; g++ )
      {
        J_dptg.col(g).col(t).col(p) = (N_asrpt.col(t).col(p).col(2)*v_asg.col(g)).sum()*rho_dpt.col(t).col(p);
        Ihat_dtg.col(g).col(t) += exp(lnqE_pg(p,g))*J_dptg.col(g).col(t).col(p);
      }

    } // next p

    meanMu_p += mu_pt.col(t);

  } // next t

  meanMu_p /= nT;

  // Observation model ----------------------------------------------------- //

  for( int t=0; t<nT; t++ )
  {
    for( int g=0; g<nG; g++ )
    {
      for( int d=0; d<nD; d++ )
      {
        // Likelihood for abundance observations
        if( !isNA(E_dtg(d,t,g)) )
        {
          // Neg-binomial will fail if prediction is too close to 0
          // Keep it above 1e-8
          Type mu = posfun( Ihat_dtg(d,t,g), Type(1e-8), varPen );
          // Neg-binomial variance
          Type var = mu+mu*mu*exp(lnDisp_tg(t,g));
          
          // Neg-binomial likelihood
          nllI_tg(t,g) -= dnbinom_robust(E_dtg(d,t,g),log(mu),log(var-mu),TRUE);
        }

        // Predicted stock composition
        for( int p=0; p<nP; p++ )
          Phat_pdtg(p,d,t,g) = exp(lnqE_pg(p,g))*J_dptg(d,p,t,g)/Ihat_dtg(d,t,g);

        // Multinomial likelihood for composition data
        if( !isNA(n_pdtg(0,d,t,g)) )
        {
          vector<Type> N_p = n_pdtg.col(g).col(t).col(d);
          vector<Type> P_p = Phat_pdtg.col(g).col(t).col(d);
          nllP_g(g) -= dmultinom( N_p, P_p, TRUE );
        }
      }
    }

    // Poisson likelihood for mark-recapture indices
    if( !isNA(I_t(t)) )
    {
      Type Isd = sqrt( log(square(mrCV_t(t))+1) );
      nllMR -= dnorm( log(I_t(t)),
                      log(mrIhat_t(t)),
                      Isd,
                      TRUE );
    }

  }

  for( int y=0; y<nY; y++ )
  {
    // Expected recruitment
    if( y < aMax )
      R_py.col(y) = R0_p;
    else
      R_py.col(y) = Z_pt.col(y-aMax)*alpha_p*exp(-beta_p*Z_pt.col(y-aMax));

    if( y > 1 )
      for( int p=0; p<nP; p++ )
        omega_py(p,y) = phi_p(p)*(lnR_py(p,y-1)-log(R_py(p,y-1)));
    
    for( int p=0; p<nP; p++ )
      nllR_py(p,y) -= dnorm( lnR_py(p,y), log(R_py(p,y)) + omega_py(p,y)+1e-10, recSD_p(p), TRUE );  
  }
  
  // Priors ---------------------------------------------------------------- //

  // Diagonal matrix with process error SDs on diagonal
  matrix<Type> D(nP,nP);
  D.fill(0);
  D.diagonal() = exp(lnErrSD_p);

  // Covariance matrix
  matrix<Type> cov_ss = D*cor_pp*D;

  // Multivariate density for process errors
  MVNORM_t<Type> errDens(cov_ss);

  // Negative-log penalty for process errors
  for( int t=0; t<(nT-1); t++ )
    nlp += errDens(arrivErr_pt.col(t));

  nlp -= UMSYPrior(2)   * dnorm( UMSY_p, UMSYPrior(0), UMSYPrior(1), TRUE ).sum();
  nlp -= phiPrior(1)    * dnorm( logitPhi_p, Type(0), phiPrior(0), TRUE ).sum();
  nlp -= pearPrior(0,2) * dnorm( pearLambda, pearPrior(0,0), pearPrior(0,1), TRUE );
  nlp -= pearPrior(1,2) * dnorm( pearTheta,  pearPrior(1,0), pearPrior(1,1), TRUE );
  nlp -= pearPrior(2,2) * dnorm( pearSigma,  pearPrior(2,0), pearPrior(2,1), TRUE );
  nlp -= pearPrior(3,2) * dnorm( pearTau,    pearPrior(3,0), pearPrior(3,1), TRUE );  

  // Total objective function
  Type objFun = nllC_rht.sum() +
                nllxB +
                nllxC_rht.sum() +
                nllR_py.sum() +
                nllMR + 
                nlp +
                weightI*nllI_tg.sum() +
                nllP_g.sum() +
                1e3*varPen;
                

  // REPORT SECTION -------------------------------------------------------- //
  ADREPORT(runSize_pt);
  ADREPORT(runSize_t);
  ADREPORT(borderPass_t);
  ADREPORT(mu_pt);
  ADREPORT(SMSY_p);
  ADREPORT(UMSY_p);
  ADREPORT(pFem_y);
  ADREPORT(phi_p);
  ADREPORT(relReproOutput_t);
  ADREPORT(recSD_p);
  ADREPORT(errSD_p);
  ADREPORT(meanMu_p);

  // Data and controls
  REPORT(n_pdtg);
  REPORT(E_dtg);
  REPORT(I_t);
  REPORT(mrCV_t);
  REPORT(day_d);
  REPORT(catchSD_rh);
  REPORT(v_asg);
  REPORT(RLM_astm);
  REPORT(mesh_rht);
  REPORT(obsC_rht);
  REPORT(xB_jgt);
  REPORT(xC_arht);
  REPORT(z_ast);
  REPORT(etaSD);
  REPORT(phiPrior);
  REPORT(pearPrior);

  // Dimensions
  REPORT(nA);
  REPORT(nP);
  REPORT(nG);
  REPORT(nT);
  REPORT(nY);
  REPORT(nD);
  REPORT(nR);
  REPORT(nH);
  REPORT(nM);

  // Parameters
  REPORT(lnArrivMu_p);
  REPORT(lnArrivSD_p);
  REPORT(arrivErr_pt);
  REPORT(logitCor_pp);
  REPORT(cor_pp);
  REPORT(lnErrSD_p);
  REPORT(lnqI_p);
  REPORT(weightI);
  REPORT(lnDisp_tg);
  REPORT(initEta_asy);
  REPORT(lnR_py);
  REPORT(lnf_rht);
  REPORT(lnUMSY_p);
  REPORT(lnSMSY_p);
  REPORT(pearPars);
  REPORT(logitPhi_p);
  REPORT(logitpFem_y);
  REPORT(pFem_y);
  REPORT(lnpFemSD);

  // Latent parameters
  //REPORT(J_dt);
  REPORT(J_dptg);
  REPORT(mu_pt);
  REPORT(rho_dpt);
  REPORT(lnqE_pg);
  REPORT(eta_aspy);
  REPORT(v_astm);
  REPORT(R_spy);
  REPORT(N_asrpt);
  REPORT(C_asrht);
  REPORT(F_arhst);
  REPORT(F_arst);
  REPORT(R_spy);
  REPORT(xB_asgt);
  REPORT(xhatC_arht);
  REPORT(R_py);
  REPORT(Z_pt);
  REPORT(alpha_p);
  REPORT(beta_p);
  REPORT(S_aspt);
  REPORT(recSD_p);
  REPORT(omega_py);
  REPORT(pearTau);
  REPORT(pearSigma);
  REPORT(pearTheta);
  REPORT(pearLambda);
  REPORT(phi_p);
  REPORT(runSize_pt);
  REPORT(relReproOutput_t);

  // Model predictions
  REPORT(Ihat_dtg);
  REPORT(mrIhat_t);
  REPORT(Phat_pdtg);

  // Objective function
  REPORT(objFun);
  REPORT(nllI_tg);
  REPORT(nllP_g);
  REPORT(nllMR);
  REPORT(nllC_rht);
  REPORT(nllxC_rht);
  REPORT(nllxB);
  REPORT(nlp);
  REPORT(nllR_py);
  REPORT(varPen);

  return objFun;
}




