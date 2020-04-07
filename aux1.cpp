// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
int cat1(double value, NumericVector prob) {
  int res=prob.length();
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function calculates distance between two sets of coordinates
// [[Rcpp::export]]
NumericMatrix GetDistance(NumericMatrix AcCoord,NumericMatrix GridCoord, int Ngrid, int Nac) {
  
  NumericMatrix res(Ngrid,Nac);
  double x2;
  double y2;
  for(int i=0; i<Ngrid;i++){
    for(int j=0; j<Nac; j++){
      x2=pow(GridCoord(i,0)-AcCoord(j,0),2.0);
      y2=pow(GridCoord(i,1)-AcCoord(j,1),2.0);
      res(i,j)=sqrt(x2+y2);
    }
  }
  return res;
}

//' This function helps calculate theta starting from stick-breaking v's
// [[Rcpp::export]]
NumericMatrix GetTheta(NumericVector v, int nac, int ntsegm){
  NumericMatrix theta(ntsegm,nac);
  double tmp;
  for(int i=0; i<ntsegm;i++){
    //for j==0
    theta(i,0)=v(i,0);
    tmp=1-v(i,0);
    for(int j=1; j<nac; j++){
      theta(i,j)=v(i,j)*tmp;
      tmp=tmp*(1-v(i,j));
    }
  }
  return theta;
}

//' This function samples z
// [[Rcpp::export]]

List SampleZ(int ntsegm, int ngrid, int nac, NumericVector z,
             IntegerMatrix dat, NumericMatrix theta, NumericMatrix ProbDist){
  //create array
  arma::cube z1(z.begin(), ntsegm, ngrid, nac, false);
  NumericVector prob(nac);
  int ind;
  //sample z
  for(int i=0; i<ntsegm; i++){
    for(int j=0; j<ngrid; j++){
      if(dat(i,j)>0){
        for(int k=0; k<nac; k++){
          prob[k]=theta(i,k)*ProbDist(k,j);          
        }
        prob=prob/sum(prob);
        
        //sample cluster memberships from multinomial
        NumericVector runif1=runif(dat(i,j));
        for(int m=0; m<dat(i,j); m++){
          ind=cat1(runif1[m], prob);  
          z1(i,j,ind)=z1(i,j,ind)+1;
        }
      }
    }
  }

  List L = List::create(Named("z") =z1);
  return(L);
}