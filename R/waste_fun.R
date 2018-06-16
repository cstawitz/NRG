#' Waste function; Fish Bioenergetics model
#' This function calculates the temperature and weight specific estimates of lost energey (waste) for an individual fish
#'
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#'
#' #' Based on the Fish Bioenergetics Model 3.0, see update 4.0 by Deslauriers, D.,
#' Chipps, S.R., Breck, J.E., Rice, J.A. & Madenjian, C.P. 2017. Fish Bioenergetics 4.0:
#' An R- based modeling application. Fisheries, 42(11): 586–596.
#'
#' This function calculates the weight and tempature specific physiological paramaters given a set of parameters and data inputs.
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#'
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$W  Weight of the fish in grams
#' @param data$TempC  Temperature
#' @param data$Weq  allometric waste equation (1,2 or 3)
#' @keywords Temperature, scaling, Respiration
#' @export Waste_fun
#' @examples
#' Waste_fun()

Waste_fun<-function(par,data){
  TempC<-data$TempC
  C_in<-data$C_in
  Weq<-par$Weq

  if (Weq==1){
    # F<-FA*C # egestion
    F<-par$FA*C_in  #g/g/d
    U<-par$UA*(C_in-F)  # g/g/d
  }else if (Weq==2){
    FB<-par$FB
    FG<-par$FG
    UB<-par$UB
    UG<-par$UG
    F<-FA*(TempC^FB)*exp(FG*RFR)*C_in  #g/g/d
    F<-UA*(TempC^UB)*exp(UG*RFR)*(C_in-F)  #g prey/g fish/d
  }else if (Weq==3){
    FB<-par$FB
    FG<-par$FG
    UB<-par$UB
    UG<-par$UG
    PFF<-sum(data$indgst*data$diet) # where indgst and diet are vectors of indgst= indigestable proportion, and diet = proportion of prey in the diet
    PE<-FA*(TempC^FB)*exp(FG*RFR)
    PF<-((PE-.1)/.9)*(1-PFF)+PFF
    F<-PF*C_in
    U<-UA*(TempC^RB)*exp(UG*RFR)*(C_in-F)
  }
  return(list(F=F,U=U))
}
