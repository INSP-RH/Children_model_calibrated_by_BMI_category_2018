//
//  child_weight.cpp
//
//  This is a function that calculates
//  weight change for children using the dynamic
//  weight model by Kevin D. Hall et al. and
//  Runge Kutta method to solve the ODE system.
//
//  Input:
//  age             .-  Years since individual first arrived to Earth
//  sex             .-  Either 1 = "female" or 0 = "male"
//  FFM             .-  Fat Free Mass (kg) of the individual
//  FM              .-  Fat Mass (kg) of the individual
//  input_EIntake   .-  Energy intake (kcal) of individual per day
//  days            .-  Days to model (integer)
//  dt              .-  Time step used to solve the ODE system numerically
//  K               .-  Richardson parameter
//  Q               .-  Richardson parameter
//  A               .-  Richardson parameter
//  B               .-  Richardson parameter
//  nu              .-  Richardson parameter
//  C               .-  Richardson parameter
//  Note:
//  Weight = FFM + FM. No extracellular fluid or glycogen is considered
//  Please see child_weight.hpp for additional information
//
//  Authors:
//  Dalia Camacho-García-Formentí
//  Rodrigo Zepeda-Tello
//
// References:
//
//  Deurenberg, Paul, Jan A Weststrate, and Jaap C Seidell. 1991. “Body Mass Index as a Measure of Body Fatness:
//      Age-and Sex-Specific Prediction Formulas.” British Journal of Nutrition 65 (2). Cambridge University Press: 105–14.
//
//  Ellis, Kenneth J, Roman J Shypailo, Steven A Abrams, and William W Wong. 2000. “The Reference Child and Adolescent Models of
//      Body Composition: A Contemporary Comparison.” Annals of the New York Academy of Sciences 904 (1). Wiley Online Library: 374–82.
//
//  Fomon, Samuel J, Ferdinand Haschke, Ekhard E Ziegler, and Steven E Nelson. 1982.
//      “Body Composition of Reference Children from Birth to Age 10 Years.” The American Journal of
//      Clinical Nutrition 35 (5). Am Soc Nutrition: 1169–75.
//
//  Hall, Kevin D, Nancy F Butte, Boyd A Swinburn, and Carson C Chow. 2013. “Dynamics of Childhood Growth
//      and Obesity: Development and Validation of a Quantitative Mathematical Model.” The Lancet Diabetes & Endocrinology 1 (2). Elsevier: 97–105.
//
//  Haschke, F. 1989. “Body Composition During Adolescence.” Body Composition Measurements in Infants and Children.
//      Ross Laboratories Columbus, OH, 76–83.
//
//  Katan, Martijn B, Janne C De Ruyter, Lothar DJ Kuijper, Carson C Chow, Kevin D Hall, and Margreet R Olthof. 2016.
//      “Impact of Masked Replacement of Sugar-Sweetened with Sugar-Free Beverages on Body Weight Increases with Initial Bmi:
//      Secondary Analysis of Data from an 18 Month Double–Blind Trial in Children.” PloS One 11 (7). Public Library of Science: e0159771.
//
//----------------------------------------------------------------------------------------
// License: MIT
// Copyright 2018 Instituto Nacional de Salud Pública de México
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute,
// sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
// is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies
// or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
// BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//----------------------------------------------------------------------------------------


#include "child_weight.h"

//Default (classic) constructor for energy matrix
Child::Child(NumericVector input_age, NumericVector input_sex, NumericVector input_bmiCat, NumericVector input_FFM, NumericVector input_FM, NumericMatrix input_EIntake,
             double input_dt, bool checkValues, double input_referenceValues){
    age   = input_age;
    sex   = input_sex;
    bmiCat = input_bmiCat;
    FM    = input_FM;
    FFM   = input_FFM;
    dt    = input_dt;
    EIntake = input_EIntake;
    check = checkValues;
    generalized_logistic = false;
    referenceValues = input_referenceValues;
    build();
}

//Constructor which uses Richard's curve with the parameters of https://en.wikipedia.org/wiki/Generalised_logistic_function
Child::Child(NumericVector input_age, NumericVector input_sex,  NumericVector input_bmiCat, NumericVector input_FFM, NumericVector input_FM, double input_K,
             double input_Q, double input_A, double input_B, double input_nu, double input_C, 
             double input_dt, bool checkValues, double input_referenceValues){
    age   = input_age;
    sex   = input_sex;
    bmiCat = input_bmiCat;
    FM    = input_FM;
    FFM   = input_FFM;
    dt    = input_dt;
    K_logistic = input_K;
    A_logistic = input_A;
    Q_logistic = input_Q;
    B_logistic = input_B;
    nu_logistic = input_nu;
    C_logistic = input_C;
    check = checkValues;
    referenceValues = input_referenceValues;
    generalized_logistic = true;
    ;
    build();
}

Child::~Child(void){
    
}

void Child::build(){
    getParameters();
}

//General function for expressing growth and eb terms
NumericVector Child::general_ode(NumericVector t, NumericVector input_A, NumericVector input_B,
                                 NumericVector input_D, NumericVector input_tA,
                                 NumericVector input_tB, NumericVector input_tD,
                                 NumericVector input_tauA, NumericVector input_tauB,
                                 NumericVector input_tauD){
    
    return input_A*exp(-(t-input_tA)/input_tauA ) +
            input_B*exp(-0.5*pow((t-input_tB)/input_tauB,2)) +
            input_D*exp(-0.5*pow((t-input_tD)/input_tauD,2));
}

NumericVector Child::Growth_dynamic(NumericVector t){
    return general_ode(t, A, B, D, tA, tB, tD, tauA, tauB, tauD);
}

NumericVector Child::Growth_impact(NumericVector t){
    return general_ode(t, A1, B1, D1, tA1, tB1, tD1, tauA1, tauB1, tauD1);
}

NumericVector Child::EB_impact(NumericVector t){
    return general_ode(t, A_EB, B_EB, D_EB, tA_EB, tB_EB, tD_EB, tauA_EB, tauB_EB, tauD_EB);
}

NumericVector Child::cRhoFFM(NumericVector input_FFM){
    return 4.3*input_FFM + 837.0;
}

NumericVector Child::cP(NumericVector FFM, NumericVector FM){
    NumericVector rhoFFM = cRhoFFM(FFM);
    NumericVector C      = 10.4 * rhoFFM / rhoFM;
    return C/(C + FM);
}

NumericVector Child::Delta(NumericVector t){
    return deltamin + (deltamax - deltamin)*(1.0 / (1.0 + pow((t / P),h)));
}

NumericVector Child::FFMReference(NumericVector t){ 
  /*  return ffm_beta0 + ffm_beta1*t; */
NumericVector under = ifelse(bmiCat == 1, 1.0, 0.0);
NumericVector normales = ifelse(bmiCat == 2, 1.0, 0.0);
NumericVector over = ifelse(bmiCat == 3, 1.0, 0.0);
NumericVector obese = ifelse(bmiCat == 4, 1.0, 0.0);


  if(referenceValues == 0){
  // -------------------------- Mean values
NumericMatrix ffm_ref(17,nind);
ffm_ref(0,_)   = 10.134*(1-sex)+9.477*sex;       // 2 years old
ffm_ref(1,_)   = 12.099*(1 - sex) + 11.494*sex;    // 3 years old
ffm_ref(2,_)   = 14.0*(1 - sex) + 13.2*sex;        // 4 years old
ffm_ref(3,_)   = 15.72*(1 - sex) + 14.86*sex;      // 5 years old
ffm_ref(4,_)   = under*(13.78*(1 - sex) + 15.97*sex)   + normales*(17.59*(1 - sex) + 15.81*sex)   + over*(19.43*(1 - sex) + 18.59*sex)  + obese*(21.87*(1 - sex) + 21.12*sex);   // 6 years old   
ffm_ref(5,_)   = under*(17.59*(1 - sex) + 16.89*sex)   + normales*(18.97*(1 - sex) + 17.96*sex)   + over*(21.84*(1 - sex) + 21.07*sex)  + obese*(24.88*(1 - sex) + 25.64*sex);  // 7 years old
ffm_ref(6,_)   = under*(17.84*(1 - sex) + 18.11*sex)   + normales*(20.72*(1 - sex) + 19.99*sex)   + over*(25.18*(1 - sex) + 22.99*sex)  + obese*(28.81*(1 - sex) + 28.21*sex); // 8 years old
ffm_ref(7,_)   = under*(19.88*(1 - sex) + 16.14*sex)   + normales*(23.46*(1 - sex) + 22.13*sex)   + over*(27.45*(1 - sex) + 27.50*sex)  + obese*(32.39*(1 - sex) + 31.09*sex); // 9 years old
ffm_ref(8,_)   = under*(23.36*(1 - sex) + 23.89*sex)   + normales*(25.35*(1 - sex) + 25.22*sex)   + over*(30.94*(1 - sex) + 31.30*sex)  + obese*(35.98*(1 - sex) + 35.88*sex); // 10 years old
ffm_ref(9,_)   = under*(23.89*(1 - sex) + 21.65*sex)   + normales*(28.65*(1 - sex) + 29.40*sex)   + over*(33.65*(1 - sex) + 35.30*sex)  + obese*(39.31*(1 - sex) + 39.46*sex); // 11 years old
ffm_ref(10,_)  = under*(27.80*(1 - sex) + 26.46*sex)   + normales*(33.08*(1 - sex) + 32.61*sex)   + over*(39.48*(1 - sex) + 37.21*sex)  + obese*(44.78*(1 - sex) + 42.21*sex); // 12 years old
ffm_ref(11,_)  = under*(31.85*(1 - sex) + 28.45*sex)   + normales*(38.71*(1 - sex) + 35.03*sex)   + over*(42.83*(1 - sex) + 39.29*sex)  + obese*(47.03*(1 - sex) + 45.01*sex); // 13 years old
ffm_ref(12,_)  = under*(34.02*(1 - sex) + 34.24*sex)   + normales*(42.24*(1 - sex) + 36.52*sex)   + over*(48.24*(1 - sex) + 41.28*sex)  + obese*(54.66*(1 - sex) + 46.63*sex); // 14 years old
ffm_ref(13,_)  = under*(34.97*(1 - sex) + 33.17*sex)   + normales*(45.14*(1 - sex) + 38.67*sex)   + over*(50.03*(1 - sex) + 43.47*sex)  + obese*(55.64*(1 - sex) + 47.78*sex); // 15 years old
ffm_ref(14,_)  = under*(39.77*(1 - sex) + 31.70*sex)   + normales*(47.04*(1 - sex) + 39.64*sex)   + over*(53.71*(1 - sex) + 45.74*sex)  + obese*(58.05*(1 - sex) + 50.88*sex); // 16 years old
ffm_ref(15,_)  = under*(42.10*(1 - sex) + 33.63*sex)   + normales*(48.25*(1 - sex) + 39.85*sex)   + over*(55.36*(1 - sex) + 45.26*sex)  + obese*(60.13*(1 - sex) + 50.52*sex); // 17 years old
ffm_ref(16,_)  = under*(44.56*(1 - sex) + 35.98*sex)   + normales*(49.11*(1 - sex) + 40.92*sex)   + over*(56.32*(1 - sex) + 46.59*sex)  + obese*(61.05*(1 - sex) + 50.02*sex);    // 18 years old

  }


  if(referenceValues == 1){
    // -------------------------- Median values
NumericMatrix ffm_ref(17,nind);
ffm_ref(0,_)   = 10.134*(1-sex)+9.477*sex;       // 2 years old
ffm_ref(1,_)   = 12.099*(1 - sex) + 11.494*sex;    // 3 years old
ffm_ref(2,_)   = 14.0*(1 - sex) + 13.2*sex;        // 4 years old
ffm_ref(3,_)   = 15.72*(1 - sex) + 14.86*sex;      // 5 years old
ffm_ref(4,_)   = under*(14.58*(1-sex) + 14.61*sex) + normales*(17.28*(1-sex) + 15.67*sex) + over*(19.14*(1-sex) + 19.07*sex) + obese*(21.68*(1-sex) + 20.68*sex);   // 6 years old
ffm_ref(5,_)   = under*(18.82*(1-sex) + 16.14*sex) + normales*(18.78*(1-sex) + 17.94*sex) + over*(22.30*(1-sex) + 20.92*sex) + obese*(24.91*(1-sex) + 25.33*sex);   // 7 years old
ffm_ref(6,_)   = under*(17.26*(1-sex) + 18.20*sex) + normales*(20.44*(1-sex) + 20.16*sex) + over*(24.75*(1-sex) + 22.76*sex) + obese*(28.54*(1-sex) + 27.93*sex);   // 8 years old
ffm_ref(7,_)   = under*(19.30*(1-sex) + 16.31*sex) + normales*(23.42*(1-sex) + 21.85*sex) + over*(26.94*(1-sex) + 27.04*sex) + obese*(31.99*(1-sex) + 30.77*sex);   // 9 years old
ffm_ref(8,_)   = under*(23.89*(1-sex) + 23.89*sex) + normales*(24.99*(1-sex) + 25.32*sex) + over*(31.37*(1-sex) + 31.09*sex) + obese*(35.81*(1-sex) + 35.76*sex);   // 10 years old
ffm_ref(9,_)   = under*(23.74*(1-sex) + 21.20*sex) + normales*(28.19*(1-sex) + 29.95*sex) + over*(33.20*(1-sex) + 35.68*sex) + obese*(38.81*(1-sex) + 39.30*sex);   // 11 years old
ffm_ref(10,_)   = under*(28.13*(1-sex) + 25.50*sex) + normales*(32.71*(1-sex) + 33.00*sex) + over*(38.84*(1-sex) + 36.92*sex) + obese*(46.35*(1-sex) + 42.30*sex);   // 12 years old
ffm_ref(11,_)   = under*(32.61*(1-sex) + 28.45*sex) + normales*(38.70*(1-sex) + 35.08*sex) + over*(43.40*(1-sex) + 38.67*sex) + obese*(47.86*(1-sex) + 44.98*sex);   // 13 years old
ffm_ref(12,_)   = under*(35.03*(1-sex) + 37.22*sex) + normales*(42.27*(1-sex) + 36.28*sex) + over*(47.71*(1-sex) + 41.50*sex) + obese*(54.53*(1-sex) + 46.94*sex);   // 14 years old
ffm_ref(13,_)   = under*(30.64*(1-sex) + 32.87*sex) + normales*(44.69*(1-sex) + 38.99*sex) + over*(50.18*(1-sex) + 43.76*sex) + obese*(54.58*(1-sex) + 47.37*sex);   // 15 years old
ffm_ref(14,_)   = under*(41.86*(1-sex) + 31.44*sex) + normales*(46.71*(1-sex) + 39.61*sex) + over*(53.18*(1-sex) + 46.38*sex) + obese*(57.93*(1-sex) + 50.98*sex);   // 16 years old
ffm_ref(15,_)   = under*(42.27*(1-sex) + 34.11*sex) + normales*(48.75*(1-sex) + 39.49*sex) + over*(55.31*(1-sex) + 45.69*sex) + obese*(60.26*(1-sex) + 50.13*sex);   // 17 years old
ffm_ref(16,_)   = under*(43.32*(1-sex) + 35.98*sex) + normales*(48.78*(1-sex) + 41.66*sex) + over*(57.29*(1-sex) + 46.94*sex) + obese*(59.68*(1-sex) + 49.72*sex);   // 18 years old
  }

NumericVector ffm_ref_t(nind);
int jmin;
int jmax;
double diff;
for(int i=0;i<nind;i++){
  if(t(i)>=18.0){
    ffm_ref_t(i)=ffm_ref(16,i);
  }else{
    jmin=floor(t(i));
    jmin=std::max(jmin,2);
    jmin=jmin-2;
    jmax= std::min(jmin+1,17);
    diff= t(i)-floor(t(i));
    ffm_ref_t(i)=ffm_ref(jmin,i)+diff*(ffm_ref(jmax,i)-ffm_ref(jmin,i));
  } 
}
return ffm_ref_t;
}

NumericVector Child::FMReference(NumericVector t){
   /* return fm_beta0 + fm_beta1*t;*/
NumericVector under = ifelse(bmiCat == 1, 1.0, 0.0);
NumericVector normales = ifelse(bmiCat == 2, 1.0, 0.0);
NumericVector over = ifelse(bmiCat == 3, 1.0, 0.0);
NumericVector obese = ifelse(bmiCat == 4, 1.0, 0.0);


 if(referenceValues == 0){
  // ---------------------------------------- Mean values
NumericMatrix fm_ref(17,nind);
fm_ref(0,_)   = 2.456*(1-sex)+ 2.433*sex;       // 2 years old
fm_ref(1,_)   = 2.576*(1 - sex) + 2.606*sex;    // 3 years old
fm_ref(2,_)   = 2.7*(1 - sex) + 2.8*sex;        // 4 years old
fm_ref(3,_)   = 3.66*(1 - sex) + 4.47*sex;      // 5 years old
fm_ref(4,_)   = under*(2.02*(1 - sex) + 2.77*sex)   + normales*(3.52*(1 - sex) + 3.99*sex)   + over*(4.87*(1 - sex) + 6.01*sex)   + obese*(7.34*(1 - sex) + 9.10*sex);   // 6 years old   
fm_ref(5,_)   = under*(2.43*(1 - sex) + 2.92*sex)   + normales*(3.70*(1 - sex) + 4.51*sex)   + over*(5.43*(1 - sex) + 6.79*sex)   + obese*(8.73*(1 - sex) + 11.60*sex);  // 7 years old
fm_ref(6,_)   = under*(2.19*(1 - sex) + 3.02*sex)   + normales*(3.99*(1 - sex) + 4.89*sex)   + over*(6.30*(1 - sex) + 7.40*sex)   + obese*(10.49*(1 - sex) + 12.71*sex); // 8 years old
fm_ref(7,_)   = under*(2.54*(1 - sex) + 2.37*sex)   + normales*(4.41*(1 - sex) + 5.23*sex)   + over*(6.92*(1 - sex) + 9.09*sex)   + obese*(12.56*(1 - sex) + 14.88*sex); // 9 years old
fm_ref(8,_)   = under*(2.96*(1 - sex) + 4.00*sex)   + normales*(4.63*(1 - sex) + 6.07*sex)   + over*(8.25*(1 - sex) + 10.95*sex)  + obese*(13.87*(1 - sex) + 17.80*sex); // 10 years old
fm_ref(9,_)   = under*(2.83*(1 - sex) + 3.62*sex)   + normales*(5.32*(1 - sex) + 7.33*sex)   + over*(9.06*(1 - sex) + 12.84*sex)  + obese*(16.32*(1 - sex) + 22.61*sex); // 11 years old
fm_ref(10,_)  = under*(3.20*(1 - sex) + 4.35*sex)   + normales*(6.33*(1 - sex) + 8.59*sex)   + over*(11.39*(1 - sex) + 14.45*sex) + obese*(19.77*(1 - sex) + 24.10*sex); // 12 years old
fm_ref(11,_)  = under*(3.44*(1 - sex) + 4.38*sex)   + normales*(7.79*(1 - sex) + 9.73*sex)   + over*(12.66*(1 - sex) + 15.46*sex) + obese*(21.56*(1 - sex) + 29.22*sex); // 13 years old
fm_ref(12,_)  = under*(3.81*(1 - sex) + 5.44*sex)   + normales*(8.71*(1 - sex) + 9.89*sex)   + over*(14.96*(1 - sex) + 16.18*sex) + obese*(26.44*(1 - sex) + 27.85*sex); // 14 years old
fm_ref(13,_)  = under*(3.98*(1 - sex) + 5.17*sex)   + normales*(9.44*(1 - sex) + 10.85*sex)  + over*(16.07*(1 - sex) + 17.81*sex) + obese*(28.15*(1 - sex) + 29.34*sex); // 15 years old
fm_ref(14,_)  = under*(4.45*(1 - sex) + 4.95*sex)   + normales*(10.04*(1 - sex) + 11.16*sex)  + over*(18.43*(1 - sex) + 19.81*sex) + obese*(30.06*(1 - sex) + 32.58*sex); // 16 years old
fm_ref(15,_)  = under*(4.66*(1 - sex) + 5.19*sex)   + normales*(10.25*(1 - sex) + 10.94*sex) + over*(18.50*(1 - sex) + 19.14*sex) + obese*(30.60*(1 - sex) + 30.37*sex); // 17 years old
fm_ref(16,_)  = under*(5.07*(1 - sex) + 5.04*sex)   + normales*(10.78*(1 - sex) + 11.02*sex) + over*(19.24*(1 - sex) + 19.53*sex) + obese*(37.55*(1 - sex) + 31.50*sex);    // 18 years old
 }

 if(referenceValues == 1){
  // ---------------------------------------- Median values
NumericMatrix fm_ref(17,nind);
fm_ref(0,_)   = 2.456*(1-sex)+ 2.433*sex;       // 2 years old
fm_ref(1,_)   = 2.576*(1 - sex) + 2.606*sex;    // 3 years old
fm_ref(2,_)   = 2.7*(1 - sex) + 2.8*sex;        // 4 years old
fm_ref(3,_)   = 3.66*(1 - sex) + 4.47*sex;      // 5 years old
fm_ref(4,_)   = under*(2.22*(1-sex) + 2.57*sex) + normales*(3.44*(1-sex) + 3.97*sex) + over*(4.74*(1-sex) + 6.00*sex) + obese*(6.56*(1-sex) + 8.57*sex);   // 6 years old
fm_ref(5,_)   = under*(2.73*(1-sex) + 2.93*sex) + normales*(3.70*(1-sex) + 4.52*sex) + over*(5.52*(1-sex) + 6.73*sex) + obese*( 8.11*(1-sex) + 10.80*sex);   // 7 years old
fm_ref(6,_)   = under*(2.02*(1-sex) + 3.06*sex) + normales*(4.00*(1-sex) + 4.99*sex) + over*(6.23*(1-sex) + 7.22*sex) + obese*( 9.35*(1-sex) + 11.84*sex);   // 8 years old
fm_ref(7,_)   = under*(2.57*(1-sex) + 2.79*sex) + normales*(4.46*(1-sex) + 4.98*sex) + over*(6.74*(1-sex) + 8.82*sex) + obese*(11.91*(1-sex) + 13.08*sex);   // 9 years old
fm_ref(8,_)   = under*(3.01*(1-sex) + 4.01*sex) + normales*(4.67*(1-sex) + 5.86*sex) + over*( 8.33*(1-sex) + 10.68*sex) + obese*(14.07*(1-sex) + 16.46*sex);   // 10 years old
fm_ref(9,_)   = under*(2.76*(1-sex) + 3.57*sex) + normales*(4.92*(1-sex) + 7.33*sex) + over*( 8.96*(1-sex) + 12.53*sex) + obese*(14.80*(1-sex) + 20.96*sex);   // 11 years old
fm_ref(10,_)   = under*(3.17*(1-sex) + 4.18*sex) + normales*(6.25*(1-sex) + 8.59*sex) + over*(11.46*(1-sex) + 14.09*sex) + obese*(19.12*(1-sex) + 22.63*sex);   // 12 years old
fm_ref(11,_)   = under*(3.64*(1-sex) + 4.38*sex) + normales*(7.67*(1-sex) + 9.94*sex) + over*(12.15*(1-sex) + 14.67*sex) + obese*(22.42*(1-sex) + 28.20*sex);   // 13 years old
fm_ref(12,_)   = under*(3.77*(1-sex) + 5.88*sex) + normales*(8.51*(1-sex) + 9.54*sex) + over*(14.70*(1-sex) + 15.84*sex) + obese*(24.80*(1-sex) + 25.32*sex);   // 14 years old
fm_ref(13,_)   = under*(3.66*(1-sex) + 5.30*sex) + normales*( 9.04*(1-sex) + 10.93*sex) + over*(15.74*(1-sex) + 17.67*sex) + obese*(25.71*(1-sex) + 28.37*sex);   // 15 years old
fm_ref(14,_)   = under*(4.43*(1-sex) + 4.99*sex) + normales*( 9.87*(1-sex) + 11.09*sex) + over*(18.88*(1-sex) + 19.74*sex) + obese*(27.81*(1-sex) + 31.11*sex);   // 16 years old
fm_ref(15,_)   = under*(4.40*(1-sex) + 5.36*sex) + normales*(10.35*(1-sex) + 10.43*sex) + over*(17.69*(1-sex) + 18.50*sex) + obese*(27.69*(1-sex) + 29.70*sex);   // 17 years old
fm_ref(16,_)   = under*(5.18*(1-sex) + 5.05*sex) + normales*(10.44*(1-sex) + 11.10*sex) + over*(19.39*(1-sex) + 18.70*sex) + obese*(31.82*(1-sex) + 28.55*sex);   // 18 years old
 }



  
NumericVector fm_ref_t(nind);
int jmin;
int jmax;
double diff;
for(int i=0;i<nind;i++){
  if(t(i)>=18.0){
    fm_ref_t(i)=fm_ref(16,i);
  }else{
    jmin=floor(t(i));
    jmin=std::max(jmin,2);
    jmin=jmin-2;
    jmax= std::min(jmin+1,17);
    diff= t(i)-floor(t(i));
    fm_ref_t(i)=fm_ref(jmin,i)+diff*(fm_ref(jmax,i)-fm_ref(jmin,i));
  } 
}
return fm_ref_t;
}

NumericVector Child::IntakeReference(NumericVector t){
    NumericVector EB      = EB_impact(t);
    NumericVector FFMref  = FFMReference(t);
    NumericVector FMref   = FMReference(t);
    NumericVector delta   = Delta(t);
    NumericVector growth  = Growth_dynamic(t);
    NumericVector p       = cP(FFMref, FMref);
    NumericVector rhoFFM  = cRhoFFM(FFMref);
    return EB + K + (22.4 + delta)*FFMref + (4.5 + delta)*FMref +
                230.0/rhoFFM*(p*EB + growth) + 180.0/rhoFM*((1-p)*EB-growth);
}

NumericVector Child::Expenditure(NumericVector t, NumericVector FFM, NumericVector FM){
    NumericVector delta     = Delta(t);
    NumericVector Iref      = IntakeReference(t);
    NumericVector Intakeval = Intake(t);
    NumericVector DeltaI    = Intakeval - Iref;
    NumericVector p         = cP(FFM, FM);
    NumericVector rhoFFM    = cRhoFFM(FFM);
    NumericVector growth    = Growth_dynamic(t);
    NumericVector Expend    = K + (22.4 + delta)*FFM + (4.5 + delta)*FM +
                                0.24*DeltaI + (230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p))*Intakeval +
                                growth*(230.0/rhoFFM -180.0/rhoFM);
    return Expend/(1.0+230.0/rhoFFM *p + 180.0/rhoFM*(1.0-p));
}

//Rungue Kutta 4 method for Adult
List Child::rk4 (double days){
    
    //Initial time
    NumericMatrix k1, k2, k3, k4;
    
    //Estimate number of elements to loop into
    int nsims = floor(days/dt);
    
    //Create array of states
    NumericMatrix ModelFFM(nind, nsims + 1); //in rcpp
    NumericMatrix ModelFM(nind, nsims + 1); //in rcpp
    NumericMatrix ModelBW(nind, nsims + 1); //in rcpp
    NumericMatrix AGE(nind, nsims + 1); //in rcpp
    NumericVector TIME(nsims + 1); //in rcpp
    
    //Create initial states
    ModelFFM(_,0) = FFM;
    ModelFM(_,0)  = FM;
    ModelBW(_,0)  = FFM + FM;
    TIME(0)  = 0.0;
    AGE(_,0)  = age;
    
    //Loop through all other states
    bool correctVals = true;
    for (int i = 1; i <= nsims; i++){

        
        //Rungue kutta 4 (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
        k1 = dMass(AGE(_,i-1), ModelFFM(_,i-1), ModelFM(_,i-1));
        k2 = dMass(AGE(_,i-1) + 0.5 * dt/365.0, ModelFFM(_,i-1) + 0.5 * k1(0,_), ModelFM(_,i-1) + 0.5 * k1(1,_));
        k3 = dMass(AGE(_,i-1) + 0.5 * dt/365.0, ModelFFM(_,i-1) + 0.5 * k2(0,_), ModelFM(_,i-1) + 0.5 * k2(1,_));
        k4 = dMass(AGE(_,i-1) + dt/365.0, ModelFFM(_,i-1) + k3(0,_), ModelFM(_,i-1) +  k3(1,_));
        
        //Update of function values
        //Note: The dt is factored from the k1, k2, k3, k4 defined on the Wikipedia page and that is why
        //      it appears here.
        ModelFFM(_,i) = ModelFFM(_,i-1) + dt*(k1(0,_) + 2.0*k2(0,_) + 2.0*k3(0,_) + k4(0,_))/6.0;        //ffm
        ModelFM(_,i)  = ModelFM(_,i-1)  + dt*(k1(1,_) + 2.0*k2(1,_) + 2.0*k3(1,_) + k4(1,_))/6.0;        //fm
        
        //Update weight
        ModelBW(_,i) = ModelFFM(_,i) + ModelFM(_,i);
        
        //Update TIME(i-1)
        TIME(i) = TIME(i-1) + dt; // Currently time counts the time (days) passed since start of model
        
        //Update AGE variable
        AGE(_,i) = AGE(_,i-1) + dt/365.0; //Age is variable in years
    }
    
    return List::create(Named("Time") = TIME,
                        Named("Age") = AGE,
                        Named("Fat_Free_Mass") = ModelFFM,
                        Named("Fat_Mass") = ModelFM,
                        Named("Body_Weight") = ModelBW,
                        Named("Correct_Values")=correctVals,
                        Named("Model_Type")="Children");


}

NumericMatrix  Child::dMass (NumericVector t, NumericVector FFM, NumericVector FM){
    
    NumericMatrix Mass(2, nind); //in rcpp;
    NumericVector rhoFFM    = cRhoFFM(FFM);
    NumericVector p         = cP(FFM, FM);
    NumericVector growth    = Growth_dynamic(t);
    NumericVector expend    = Expenditure(t, FFM, FM);
    Mass(0,_)               = (1.0*p*(Intake(t) - expend) + growth)/rhoFFM;    // dFFM
    Mass(1,_)               = ((1.0 - p)*(Intake(t) - expend) - growth)/rhoFM; //dFM
    return Mass;
    
}

void Child::getParameters(void){
    
    //General constants
    rhoFM    = 9.4*1000.0;
    deltamin = 10.0;
    P        = 12.0;
    h        = 10.0;
    
    //Number of individuals
    nind     = age.size();
    
    //Sex specific constants
    ffm_beta0 = 2.9*(1 - sex)  + 3.8*sex;
    ffm_beta1 = 2.9*(1 - sex)  + 2.3*sex;
    fm_beta0  = 1.2*(1 - sex)  + 0.56*sex;
    fm_beta1  = 0.41*(1 - sex) + 0.74*sex;
    K         = 800*(1 - sex)  + 700*sex;
    deltamax  = 19*(1 - sex)   + 17*sex;
    A         = 3.2*(1 - sex)  + 2.3*sex;
    B         = 9.6*(1 - sex)  + 8.4*sex;
    D         = 10.1*(1 - sex) + 1.1*sex;
    tA        = 4.7*(1 - sex)  + 4.5*sex;       //years
    tB        = 12.5*(1 - sex) + 11.7*sex;      //years
    tD        = 15.0*(1-sex)   + 16.2*sex;      //years
    tauA      = 2.5*(1 - sex)  + 1.0*sex;       //years
    tauB      = 1.0*(1 - sex)  + 0.9*sex;       //years
    tauD      = 1.5*(1 - sex)  + 0.7*sex;       //years
    A_EB      = 7.2*(1 - sex)  + 16.5*sex;
    B_EB      = 30*(1 - sex)   + 47.0*sex;
    D_EB      = 21*(1 - sex)   + 41.0*sex;
    tA_EB     = 5.6*(1 - sex)  + 4.8*sex;
    tB_EB     = 9.8*(1 - sex)  + 9.1*sex;
    tD_EB     = 15.0*(1 - sex) + 13.5*sex;
    tauA_EB   = 15*(1 - sex)   + 7.0*sex;
    tauB_EB   = 1.5*(1 -sex)   + 1.0*sex;
    tauD_EB   = 2.0*(1 - sex)  + 1.5*sex;
    A1        = 3.2*(1 - sex)  + 2.3*sex;
    B1        = 9.6*(1 - sex)  + 8.4*sex;
    D1        = 10.0*(1 - sex) + 1.1*sex;
    tA1       = 4.7*(1 - sex)  + 4.5*sex;
    tB1       = 12.5*(1 - sex) + 11.7*sex;
    tD1       = 15.0*(1 - sex) + 16.0*sex;
    tauA1     = 1.0*(1 - sex)  + 1.0*sex;
    tauB1     = 0.94*(1 - sex) + 0.94*sex;
    tauD1     = 0.69*(1 - sex) + 0.69*sex;
}


//Intake in calories
NumericVector Child::Intake(NumericVector t){
    if (generalized_logistic) {
        return A_logistic + (K_logistic - A_logistic)/pow(C_logistic + Q_logistic*exp(-B_logistic*t), 1/nu_logistic); //t in years
    } else {
        int timeval = floor(365.0*(t(0) - age(0))/dt); //Example: Age: 6 and t: 7.1 => timeval = 401 which corresponds to the 401 entry of matrix
        return EIntake(timeval,_);
    }
    
}
