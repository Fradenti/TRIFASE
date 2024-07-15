#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;



// [[Rcpp::export]]
Rcpp::List TRIFASE_cpp(arma::mat X,
                       arma::mat Fm,
                       arma::mat Gm,
                       arma::mat L,
                       arma::mat iKernel, 
                       double current_tau,
                       arma::mat Linv,
                       arma::mat Dinv,
                       int K, int R,
                       int niter,
                       int swipes = 10,
                       double threshold = 0.01,
                       double nugget = 1e-9,
                       bool estimate_tau = true,
                       std::string upd_row  = "C", // Stoc or Classif 
                       std::string upd_col = "S"  // Stoc, Appr, Classif
){
  
  arma::vec possible_label_col = arma::linspace(0, R-1, R);
  arma::vec possible_label_row = arma::linspace(0, K-1, K);
  
  arma::colvec loss(niter);
  arma::colvec taus(niter);
  loss.zeros();
  arma::mat mu(K,R);
  arma::mat tFm = Fm.t();
  arma::mat tGm = Gm.t();
  arma::mat tX = X.t();
  arma::mat XDinv = X * Dinv;
  arma::mat Linvt = (Linv).t();
  arma::mat X_hat =  X * (Linvt);
  arma::mat tX_hat =  (X_hat.t());
  arma::vec losses(R); 
  double min_loss = 1e+9;
  int ind_min_loss = 0;
  arma::mat final_F=Fm;
  arma::mat final_G=Gm;
  
  
  int ind_row = 0;
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat Eye_R1 = arma::eye(R,R);
  arma::mat Eye_R = arma::eye(R,R) * nugget;
  arma::mat Eye_K = arma::eye(K,K) * nugget;
  int temp = niter-1;
  
  for(int i =0; i<niter; i++){
    Rcout << "Starting iteration "<< i <<";\n";
    Rcpp::checkUserInterrupt();
    
    
    // Update mu (centroids)
    mu = arma::inv( tFm * Fm + Eye_K) * tFm * 
      XDinv * 
      Gm * arma::inv(tGm * Dinv * Gm + Eye_R);
    
    // Update F
    arma::mat G_tilde =  mu * (Gm).t() * Linvt ;
    Fm.zeros();
    for(int r = 0; r < n; r++){
      arma::mat BOH = arma::repelem(  X_hat.row(r), K, 1);
      //Rcout << BOH;
      arma::colvec LossF = arma::sum(pow(BOH - G_tilde,2),1);
      if( upd_row == "S" ){
        arma::vec pp1 = exp( - LossF - max( - LossF ));
        ind_row = RcppArmadillo::sample(possible_label_row, 1, TRUE, pp1)[0]; 
      }else if(upd_row == "C"){
        ind_row = LossF.index_min();
      }else{
        Rcout << "row_upd can be either S or C";
        break;
      }
      //
      
      Fm(r,ind_row) = 1;
    }
    tFm = Fm.t();
    
    arma::mat F_tilde = (Fm * mu);
    
    // Update G
    if(upd_col == "S"){
      for(int repl = 0; repl<swipes; repl++){
        for(int q = p-1; q > -1; q--){
          arma::mat Gtempt = tGm;
          for(int qq = 0; qq < R; qq++){
            Gtempt.col(q) = Eye_R1.col(qq);
            losses(qq) = arma::accu( pow( X_hat - F_tilde * Gtempt * Linvt, 2 ) );
          }
          arma::vec pp2 = exp( -losses - max( - losses ));
          int ind = RcppArmadillo::sample(
            possible_label_col, 1, TRUE, pp2)[0];
          //arma::uword ind2 = losses.index_min();
          (tGm.col(q)).zeros();
          tGm(ind,q) = 1;
        }
      }
      Gm = tGm.t();
    }else if(upd_col == "C"){
      
      for(int repl = 0; repl<swipes; repl++){
        for(int q = p-1; q > -1; q--){
          arma::mat Gtempt = tGm;
          for(int qq = 0; qq < R; qq++){
            Gtempt.col(q) = Eye_R1.col(qq);
            losses(qq) = arma::accu( pow( X_hat - F_tilde * Gtempt * Linvt, 2 ) );
          }
          arma::uword ind2 = losses.index_min();
          (tGm.col(q)).zeros();
          tGm(ind2,q) = 1;
        }
      }
      Gm = tGm.t();
      
    }else if(upd_col == "A"){
      Gm.zeros();
      for(int q = 0; q < p; q++){
        arma::mat BOH2 = arma::repelem(  tX.row(q), R, 1);
        arma::colvec LossG = arma::sum(pow(BOH2 - F_tilde.t(),2),1); 
        arma::uword ind2 = LossG.index_min();
        Gm(q,ind2) = 1;
      }
      tGm = Gm.t();
    }else{
      Rcout << "row_upd can be only S, C, A";
      break;
    }
    arma::mat RES = (X - Fm * mu * tGm);
    if(estimate_tau){
      double tau_new = arma::trace(  RES * 
                                  iKernel * 
                                  RES.t())/(p*n);
    
    // recompute quantities with new tau
      //////
      Dinv =  Dinv * current_tau/tau_new;
      
      XDinv = X * (Dinv);//Dinv * current_tau/tau_new;
      ////////////////////////////////////////////////////
      ////////////////////////////////////////////////////
      Linvt = (Linv*sqrt( current_tau/tau_new)).t();
      X_hat =  X * (Linvt);
      tX_hat =  (X_hat.t());
      
      current_tau = tau_new;
      taus(i) = current_tau;
    }else{
      taus(i) = current_tau;
    }
    loss(i) = arma::accu( 
      pow((X_hat) - F_tilde * tGm * Linvt,2));
    Rcout << "Loss value: "<< loss(i) <<";\n";
    
    if(i>0){
      Rcout << "Delta Loss: "<< loss(i)-loss(i-1) <<";\n";
      Rcout << "------------------------------\n";
      Rcout << "Tau: "<< taus(i-1) <<";\n";
      Rcout << "------------------------------\n";
      
      if(loss(i)<min_loss){
        min_loss = loss(i);
        ind_min_loss = i;
        final_F = Fm;
        final_G = Gm;
      }
      
      if( fabs( loss(i) - loss(i-1) )/fabs( loss(i-1) ) < threshold ){ 
        Rcout << "Convergence reached in " << i << " iterations\n";
        temp = i;
        break;
      }
    }
    
  }
  
  arma::colvec loss2 = loss.rows(1,temp);
  arma::colvec taus2 = taus.rows(1,temp);
  
  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["mu"]  = mu,
    Rcpp::_["ind_minloss"]  = ind_min_loss+1,
    Rcpp::_["F"] = final_F,
    Rcpp::_["G"] = final_G,
    Rcpp::_["loss"] = loss2,
    Rcpp::_["tau_est"] = current_tau,
    Rcpp::_["tau_evol"] = taus2);
  
  return(results);
}



