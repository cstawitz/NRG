#' Respiration function; Fish Bioenergetics model
#' This function calculates the temperature and weight specific metabolic demand for an individual fish
#'
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#'
#' #' Based on the Fish Bioenergetics Model 3.0, see update 4.0 by Deslauriers, D.,
#' Chipps, S.R., Breck, J.E., Rice, J.A. & Madenjian, C.P. 2017. Fish Bioenergetics 4.0:
#' An R- based modeling application. Fisheries, 42(11): 586–596.
#'
#' This function calculates the weight and tempature specific physiological rates given a set of parameters and data inputs.
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#'
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$W  weight of the fish in grams
#' @param data$Am  weight of the fish in grams
#' @param data$TempC  Temperature
#' @keywords Temperature, scaling, Respiration
#' @export Resp_fun
#' @examples
#' Resp_fun()

Resp_fun<-function(par,data){
  TempC<-data$TempC
  Vel<-NA
  Req<-par$Req
  W<-data$W
  Am<-par$Am
  if (Req==0){
    fTr<-data$fTrmodel(TempC)
  }else if (Req==1){
    RK1<-par$RK1
    RK4<-par$RK4
    Bact<-par$Bact
    Trl<-par$Trl
    Tro<-par$Tro
    fTr<-exp(par$QR*TempC)
    # Vel<-rep(0,length(TempC))
    # Vel[TempC>Trl]<-(RK1*(W^RK4))*exp(Bact*Trl)
    # Vel[TempC<=Trl]<-(Am*(W^RK4))*exp(Bact*TempC)
    Vel<-rep(0,length(TempC))
    if(any(TempC>Trl)) Vel[TempC>Trl]<-(RK1*(W^RK4))*exp(Bact*Trl)
    if(any(TempC<=Trl)) Vel[TempC<=Trl]<-(Am*(W^RK4))*exp(Bact*TempC[TempC<=Trl])
    Act<-exp(Tro*Vel)
  }else if (Req==2){
    Tro<-par$Tro
    Trm<-par$Trm
    Vr<-(Trm-TempC)/(Trm-Tro)
    Zr<-log(par$QR)*(Trm-Tro)
    Yr<-log(par$QR)*(Trm-Tro+2)
    Xr<-(Zr^2*(1+(1+(40/Yr))^0.5)^2)/400
    fTr<-Vr^Xr*exp(Xr*(1-Vr))
    Act<-par$Am
  }else{
    message("ERROR: Req does not match criteria [0,2]")
  }
  if(data$fitLL=="vel"){
    LL<-dnorm((Vel)-(data$velobs),par$sigma,log=data$log)
    LLuse<--sum(LL)
    if(is.na(LLuse)){LLuse<-1e6}
    # if(Tcm>35){LLuse<-1e6}
    return(LLuse)
  }else if(data$fitLL=="fT"){
    LL<-dnorm((fTr)-(data$fT_obs),par$sigma,log=data$log)
    LLuse<--sum(LL)
    if(is.na(LLuse)){LLuse<-1e6}
    # if(Tcm>35){LLuse<-1e6}
    return(LLuse)
  }else {
    return(list(Act=Act,fTr=fTr,Vel=Vel))
  }
}
