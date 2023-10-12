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

NumericMatrix ffm_ref(17,nind);
  if(referenceValues == 0){
  // -------------------------- Mean values
ffm_ref(0,_)   = 10.134*(1-sex)+9.477*sex;       // 2 years old
ffm_ref(1,_)   = 12.099*(1 - sex) + 11.494*sex;    // 3 years old
ffm_ref(2,_)   = 14.0*(1 - sex) + 13.2*sex;        // 4 years old
ffm_ref(3,_)   = 15.72*(1 - sex) + 14.86*sex;      // 5 years old
ffm_ref(4,_)   = under*(12.7942*(1-sex) + 13.7957*sex) + normales*(17.0238*(1-sex) + 15.2337*sex) + over*(19.3070*(1-sex) + 17.7866*sex) + obese*(22.2248*(1-sex) + 21.2170*sex);   // 6 years old
ffm_ref(5,_)   = under*(17.8106*(1-sex) + 18.4835*sex) + normales*(19.0775*(1-sex) + 17.5198*sex) + over*(20.3344*(1-sex) + 18.9406*sex) + obese*(23.1765*(1-sex) + 22.2733*sex);   // 7 years old
ffm_ref(6,_)   = under*(20.3597*(1-sex) + 18.5363*sex) + normales*(20.4774*(1-sex) + 19.6317*sex) + over*(22.1128*(1-sex) + 21.6080*sex) + obese*(25.8151*(1-sex) + 25.1641*sex);   // 8 years old
ffm_ref(7,_)   = under*(19.3668*(1-sex) + 17.0314*sex) + normales*(22.3768*(1-sex) + 21.3680*sex) + over*(26.7714*(1-sex) + 26.1791*sex) + obese*(31.3143*(1-sex) + 30.1484*sex);   // 9 years old
ffm_ref(8,_)   = under*(23.9096*(1-sex) + 19.1085*sex) + normales*(24.8998*(1-sex) + 24.0922*sex) + over*(30.4866*(1-sex) + 30.3541*sex) + obese*(34.1717*(1-sex) + 35.2838*sex);   // 10 years old
ffm_ref(9,_)   = under*(23.5033*(1-sex) + 23.3318*sex) + normales*(27.5943*(1-sex) + 28.2737*sex) + over*(32.6556*(1-sex) + 34.1915*sex) + obese*(38.2638*(1-sex) + 37.0428*sex);   // 11 years old
ffm_ref(10,_)   = under*(24.7662*(1-sex) + 25.9357*sex) + normales*(31.5163*(1-sex) + 31.9490*sex) + over*(37.5262*(1-sex) + 37.0654*sex) + obese*(42.3513*(1-sex) + 42.5446*sex);   // 12 years old
ffm_ref(11,_)   = under*(28.9497*(1-sex) + 30.2351*sex) + normales*(36.3432*(1-sex) + 34.3348*sex) + over*(41.6549*(1-sex) + 39.1559*sex) + obese*(48.1398*(1-sex) + 44.0205*sex);   // 13 years old
ffm_ref(12,_)   = under*(33.9297*(1-sex) + 33.6380*sex) + normales*(40.9730*(1-sex) + 36.1797*sex) + over*(48.0671*(1-sex) + 40.9960*sex) + obese*(50.1084*(1-sex) + 46.0726*sex);   // 14 years old
ffm_ref(13,_)   = under*(35.2601*(1-sex) + 33.0539*sex) + normales*(43.7795*(1-sex) + 38.1065*sex) + over*(49.3493*(1-sex) + 42.8965*sex) + obese*(55.6289*(1-sex) + 48.6841*sex);   // 15 years old
ffm_ref(14,_)   = under*(40.5041*(1-sex) + 32.9676*sex) + normales*(46.9540*(1-sex) + 40.1114*sex) + over*(52.9435*(1-sex) + 45.6216*sex) + obese*(58.9917*(1-sex) + 49.7917*sex);   // 16 years old
ffm_ref(15,_)   = under*(42.0445*(1-sex) + 32.3827*sex) + normales*(47.8972*(1-sex) + 39.6064*sex) + over*(55.8888*(1-sex) + 46.1784*sex) + obese*(58.7117*(1-sex) + 51.0534*sex);   // 17 years old
ffm_ref(16,_)   = under*(44.0779*(1-sex) + 35.5248*sex) + normales*(49.6930*(1-sex) + 41.2798*sex) + over*(56.5725*(1-sex) + 45.9979*sex) + obese*(61.7620*(1-sex) + 49.8746*sex);   // 18 years old

  }


  if(referenceValues == 1){
    // -------------------------- Median values
ffm_ref(0,_)   = 10.134*(1-sex)+9.477*sex;       // 2 years old
ffm_ref(1,_)   = 12.099*(1 - sex) + 11.494*sex;    // 3 years old
ffm_ref(2,_)   = 14.0*(1 - sex) + 13.2*sex;        // 4 years old
ffm_ref(3,_)   = 15.72*(1 - sex) + 14.86*sex;      // 5 years old
ffm_ref(4,_)   = under*(14.4641*(1-sex) + 13.8627*sex) + normales*(17.1430*(1-sex) + 15.1282*sex) + over*(19.2280*(1-sex) + 17.6859*sex) + obese*(21.9501*(1-sex) + 20.4992*sex);   // 6 years old
ffm_ref(5,_)   = under*(16.3729*(1-sex) + 16.6347*sex) + normales*(18.2285*(1-sex) + 17.2507*sex) + over*(21.7099*(1-sex) + 20.0341*sex) + obese*(24.9713*(1-sex) + 23.4162*sex);   // 7 years old
ffm_ref(6,_)   = under*(18.0019*(1-sex) + 17.2583*sex) + normales*(19.9148*(1-sex) + 19.4286*sex) + over*(24.6404*(1-sex) + 22.1758*sex) + obese*(27.4774*(1-sex) + 26.8346*sex);   // 8 years old
ffm_ref(7,_)   = under*(19.2548*(1-sex) + 17.5150*sex) + normales*(21.9058*(1-sex) + 21.2721*sex) + over*(26.5243*(1-sex) + 25.6952*sex) + obese*(30.8636*(1-sex) + 29.2900*sex);   // 9 years old
ffm_ref(8,_)   = under*(23.9096*(1-sex) + 20.1493*sex) + normales*(24.8603*(1-sex) + 23.6199*sex) + over*(29.9298*(1-sex) + 29.5716*sex) + obese*(34.1859*(1-sex) + 34.1346*sex);   // 10 years old
ffm_ref(9,_)   = under*(23.7557*(1-sex) + 24.0089*sex) + normales*(27.4756*(1-sex) + 28.2708*sex) + over*(32.4980*(1-sex) + 32.8672*sex) + obese*(38.1778*(1-sex) + 37.5833*sex);   // 11 years old
ffm_ref(10,_)   = under*(24.1310*(1-sex) + 25.5209*sex) + normales*(31.2494*(1-sex) + 32.2679*sex) + over*(37.7967*(1-sex) + 36.7435*sex) + obese*(42.8213*(1-sex) + 42.2971*sex);   // 12 years old
ffm_ref(11,_)   = under*(28.2941*(1-sex) + 32.6849*sex) + normales*(36.0685*(1-sex) + 33.7855*sex) + over*(41.4671*(1-sex) + 38.6218*sex) + obese*(48.1462*(1-sex) + 43.5195*sex);   // 13 years old
ffm_ref(12,_)   = under*(33.7396*(1-sex) + 37.2420*sex) + normales*(40.9866*(1-sex) + 35.9762*sex) + over*(47.9945*(1-sex) + 40.9744*sex) + obese*(50.9872*(1-sex) + 45.6421*sex);   // 14 years old
ffm_ref(13,_)   = under*(35.7472*(1-sex) + 32.2773*sex) + normales*(44.0430*(1-sex) + 38.2639*sex) + over*(49.7454*(1-sex) + 43.1117*sex) + obese*(54.9071*(1-sex) + 48.1360*sex);   // 15 years old
ffm_ref(14,_)   = under*(41.8846*(1-sex) + 33.0258*sex) + normales*(46.8444*(1-sex) + 39.6752*sex) + over*(53.3482*(1-sex) + 45.7056*sex) + obese*(58.5851*(1-sex) + 48.9594*sex);   // 16 years old
ffm_ref(15,_)   = under*(42.6661*(1-sex) + 31.6275*sex) + normales*(48.2625*(1-sex) + 39.5399*sex) + over*(55.9614*(1-sex) + 47.2530*sex) + obese*(58.4194*(1-sex) + 50.7464*sex);   // 17 years old
ffm_ref(16,_)   = under*(42.8578*(1-sex) + 37.5435*sex) + normales*(49.4174*(1-sex) + 41.5349*sex) + over*(56.7387*(1-sex) + 45.9623*sex) + obese*(63.6968*(1-sex) + 50.0229*sex);   // 18 years old
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

NumericMatrix fm_ref(17,nind);
 if(referenceValues == 0){
  // ---------------------------------------- Mean values

fm_ref(0,_)   = 2.456*(1-sex)+ 2.433*sex;       // 2 years old
fm_ref(1,_)   = 2.576*(1 - sex) + 2.606*sex;    // 3 years old
fm_ref(2,_)   = 2.7*(1 - sex) + 2.8*sex;        // 4 years old
fm_ref(3,_)   = 3.66*(1 - sex) + 4.47*sex;      // 5 years old
fm_ref(4,_)   = under*(1.7764*(1-sex) + 2.5951*sex) + normales*(3.4540*(1-sex) + 3.8303*sex) + over*(4.8055*(1-sex) + 5.7014*sex) + obese*(7.9672*(1-sex) + 9.3883*sex);   // 6 years old
fm_ref(5,_)   = under*(2.3398*(1-sex) + 2.8164*sex) + normales*(3.5859*(1-sex) + 4.2782*sex) + over*(5.4625*(1-sex) + 6.5960*sex) + obese*( 8.4350*(1-sex) + 10.4148*sex);   // 7 years old
fm_ref(6,_)   = under*(3.2767*(1-sex) + 3.0828*sex) + normales*(4.1138*(1-sex) + 5.2226*sex) + over*(5.5455*(1-sex) + 7.3667*sex) + obese*( 9.3266*(1-sex) + 12.0550*sex);   // 8 years old
fm_ref(7,_)   = under*(2.3902*(1-sex) + 2.6538*sex) + normales*(4.1705*(1-sex) + 5.0218*sex) + over*(6.6958*(1-sex) + 8.6945*sex) + obese*(11.5896*(1-sex) + 14.1436*sex);   // 9 years old
fm_ref(8,_)   = under*(2.9954*(1-sex) + 3.1389*sex) + normales*(4.5465*(1-sex) + 5.7742*sex) + over*( 8.1191*(1-sex) + 10.6667*sex) + obese*(13.4114*(1-sex) + 17.3329*sex);   // 10 years old
fm_ref(9,_)   = under*(2.6803*(1-sex) + 3.8049*sex) + normales*(5.0225*(1-sex) + 6.9162*sex) + over*( 8.7335*(1-sex) + 12.3291*sex) + obese*(15.2821*(1-sex) + 19.0058*sex);   // 11 years old
fm_ref(10,_)   = under*(2.8835*(1-sex) + 4.2002*sex) + normales*(5.9324*(1-sex) + 8.2706*sex) + over*(10.5608*(1-sex) + 14.4379*sex) + obese*(18.3024*(1-sex) + 24.9390*sex);   // 12 years old
fm_ref(11,_)   = under*(3.1579*(1-sex) + 4.7942*sex) + normales*(7.0763*(1-sex) + 9.1606*sex) + over*(12.3945*(1-sex) + 15.0401*sex) + obese*(21.7342*(1-sex) + 28.2547*sex);   // 13 years old
fm_ref(12,_)   = under*(3.6857*(1-sex) + 5.3309*sex) + normales*( 8.3966*(1-sex) + 10.0249*sex) + over*(15.0498*(1-sex) + 17.1050*sex) + obese*(24.2628*(1-sex) + 29.7700*sex);   // 14 years old
fm_ref(13,_)   = under*(3.9803*(1-sex) + 5.2442*sex) + normales*( 9.0181*(1-sex) + 10.5653*sex) + over*(15.5611*(1-sex) + 17.5730*sex) + obese*(27.0142*(1-sex) + 29.9077*sex);   // 15 years old
fm_ref(14,_)   = under*(4.6019*(1-sex) + 4.8228*sex) + normales*(10.0921*(1-sex) + 11.4444*sex) + over*(18.1619*(1-sex) + 19.9088*sex) + obese*(30.8170*(1-sex) + 31.2351*sex);   // 16 years old
fm_ref(15,_)   = under*(4.8405*(1-sex) + 4.8583*sex) + normales*(10.0547*(1-sex) + 10.6654*sex) + over*(19.2423*(1-sex) + 19.4731*sex) + obese*(30.7942*(1-sex) + 31.1807*sex);   // 17 years old
fm_ref(16,_)   = under*(4.6858*(1-sex) + 5.3332*sex) + normales*(10.7726*(1-sex) + 11.3437*sex) + over*(19.1356*(1-sex) + 19.0598*sex) + obese*(35.6945*(1-sex) + 30.3288*sex);   // 18 years old

 }

 if(referenceValues == 1){
  // ---------------------------------------- Median values

fm_ref(0,_)   = 2.456*(1-sex)+ 2.433*sex;       // 2 years old
fm_ref(1,_)   = 2.576*(1 - sex) + 2.606*sex;    // 3 years old
fm_ref(2,_)   = 2.7*(1 - sex) + 2.8*sex;        // 4 years old
fm_ref(3,_)   = 3.66*(1 - sex) + 4.47*sex;      // 5 years old
fm_ref(4,_)   = under*(2.0359*(1-sex) + 2.5660*sex) + normales*(3.4642*(1-sex) + 3.7042*sex) + over*(4.6220*(1-sex) + 5.6735*sex) + obese*(7.1058*(1-sex) + 8.7339*sex);   // 6 years old
fm_ref(5,_)   = under*(2.3771*(1-sex) + 2.9560*sex) + normales*(3.6030*(1-sex) + 4.1865*sex) + over*(5.5651*(1-sex) + 6.4374*sex) + obese*(8.0501*(1-sex) + 9.3100*sex);   // 7 years old
fm_ref(6,_)   = under*(2.1231*(1-sex) + 3.0917*sex) + normales*(3.6729*(1-sex) + 4.8531*sex) + over*(5.8971*(1-sex) + 7.0172*sex) + obese*( 8.9372*(1-sex) + 11.5469*sex);   // 8 years old
fm_ref(7,_)   = under*(2.4068*(1-sex) + 2.9027*sex) + normales*(4.0597*(1-sex) + 4.8707*sex) + over*(6.5720*(1-sex) + 8.7112*sex) + obese*(10.8084*(1-sex) + 12.7559*sex);   // 9 years old
fm_ref(8,_)   = under*(2.9954*(1-sex) + 3.1757*sex) + normales*(4.5932*(1-sex) + 5.4455*sex) + over*( 8.0701*(1-sex) + 10.6143*sex) + obese*(12.3133*(1-sex) + 15.7121*sex);   // 10 years old
fm_ref(9,_)   = under*(2.7443*(1-sex) + 3.8911*sex) + normales*(4.7619*(1-sex) + 6.9604*sex) + over*( 8.6445*(1-sex) + 11.7518*sex) + obese*(14.4743*(1-sex) + 17.4123*sex);   // 11 years old
fm_ref(10,_)   = under*(2.8190*(1-sex) + 4.1099*sex) + normales*(5.5671*(1-sex) + 8.3722*sex) + over*(10.2431*(1-sex) + 14.7437*sex) + obese*(17.3155*(1-sex) + 22.9359*sex);   // 12 years old
fm_ref(11,_)   = under*(3.0059*(1-sex) + 5.3651*sex) + normales*(6.7689*(1-sex) + 9.2549*sex) + over*(12.0232*(1-sex) + 14.6163*sex) + obese*(21.0382*(1-sex) + 26.6716*sex);   // 13 years old
fm_ref(12,_)   = under*(3.7104*(1-sex) + 5.8580*sex) + normales*(8.4317*(1-sex) + 9.8827*sex) + over*(15.2507*(1-sex) + 16.2256*sex) + obese*(22.9540*(1-sex) + 27.6643*sex);   // 14 years old
fm_ref(13,_)   = under*(4.4546*(1-sex) + 5.2493*sex) + normales*( 8.7820*(1-sex) + 10.3785*sex) + over*(15.6754*(1-sex) + 17.3977*sex) + obese*(25.5113*(1-sex) + 28.0559*sex);   // 15 years old
fm_ref(14,_)   = under*(4.6585*(1-sex) + 4.8742*sex) + normales*( 9.5728*(1-sex) + 11.4776*sex) + over*(18.3549*(1-sex) + 19.7533*sex) + obese*(29.9916*(1-sex) + 30.6943*sex);   // 16 years old
fm_ref(15,_)   = under*(4.8189*(1-sex) + 4.7975*sex) + normales*(10.3426*(1-sex) + 10.3454*sex) + over*(18.9543*(1-sex) + 19.3869*sex) + obese*(27.2116*(1-sex) + 29.9799*sex);   // 17 years old
fm_ref(16,_)   = under*(4.5259*(1-sex) + 5.7815*sex) + normales*(10.7497*(1-sex) + 10.9042*sex) + over*(18.9053*(1-sex) + 19.1592*sex) + obese*(31.9253*(1-sex) + 28.3702*sex);   // 18 years old
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
