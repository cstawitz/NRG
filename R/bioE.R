
#' Fish Bioenergetics model
#'
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#' Based on the Fish Bioenergetics Model 3.0, see update 4.0 by Deslauriers, D.,
#' Chipps, S.R., Breck, J.E., Rice, J.A. & Madenjian, C.P. 2017. Fish Bioenergetics 4.0:
#' An R- based modeling application. Fisheries, 42(11): 586–596.
#'
#' This function calculates the weight and tempature specific physiological rates given a set of parameters and data inputs. Predicts consumption or growth from a fit or set pvalue using model.type
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#' @param par is a list of parameters used for the model
#' @param par$Qox oxycaloric value of 13.56 J O2mg-1 from Brett & Groves 1979 Qox for Fish
#' @param Ceq Consumption equation
#' @param RFR Relative foraging rate
#' @param RA Respiration intercept
#' @param RB Respiration slope
#' @param SA Specific dynamic action
#' @param CA Intercept of the consumption function
#' @param CB Slope of the consumption function
#' @param par$Vel velocity parameter
#' @param data is a list with inputs for the model including a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param W  weight of the fish in grams
#' @param TempC  is the temperature the fish experiences (in deg C)
#' @param Eprey  is the energey density of the prey in Joules per gram of prey
#' @param Epred  is the energey density of the fish predator in Joules per gram of pred
#' @param indgst  is the indigestable proportion of each prey item in the diet
#' @param diet  is the proportion by weight of each prey item in the diet
#' @param fTCmodel  is the optional user specified function for the temperature scaling component of C; set to NA by default
#' @param fTrmodel  the optional user specified function for the temperature scaling component of R; set to NA by default
#' @param velmodel  the optional user specified function for the temperature scaling component of swim velocity; set to NA by default
#' @keywords Temperature, scaling, consumption
#' @export bioE
#' @examples
#' plk_par<-data.frame(RFR=1, Qox=13560,Ceq=2,Req=2,Weq=1,Tco=10,Tcm=15,QC=2.6,CA=0.119,CB=-0.46,
#'                     RA=0.0075,RB=-0.251, QR=2.6,Tro=13,Trm=18,SA=0.125, Am=2,FA=0.15,UA=0.11)
#' ebs_data<-list(W=100,TempC=0:10,Eprey=5539.6,Epred=4184,indgst=0,diet=0,fTCmodel=NA,fTrmodel=NA,velmodel=NA)
#' bioE(data=ebs_data,par=plk_par)
#' fTfun()

bioE<-function(par,data=list(W,TempC,Eprey,Epred,indgst,diet,fTCmodel=NA,fTrmodel=NA,velmodel=NA)){
  # data
  W<-data$W
  # TempC<-data$TempC
  Eprey<-data$Eprey
  Epred<-data$Epred
  indgst<-data$indgst
  diet<-data$diet

  ### PARAMETERS
  RFR<-par$RFR
  RA<-par$RA
  RB<-par$RB
  SA<-par$SA
  CA<-par$CA
  CB<-par$CB
  Vel<-Act<-NA
  Qox<-par$Qox

  # Consumption
  fTc<-fTC_fun(par=par,data=data)
  #Cmax_ggd<-CA*(W^CB)*fTc 	      # max cosumption g prey per g fish/d
  Cmax_jgd<-CA*(W^CB)*fTc*Eprey 	# max cosumption joules per g fish/d
  C_jgd<-CA*(W^CB)*fTc*Eprey*RFR  # consumption joules per g fish/d
  # Waste
  wdat<-data;wdat$C_in<-C_jgd
  Waste<-Waste_fun(par=par,data=wdat)
  F_jgd<-Waste$F 		#eggestion g prey/g fish/d
  U_jgd<-Waste$U		#excretion g prey/g fish/d
  # specific dynamic action j/g fish/ d
  SDA_jgd<-SA*(C_jgd-F_jgd)
  # Metabolism
  r_data<-data;r_data$fitLL<-FALSE
  Resp<-Resp_fun(par=par,data=r_data)
  fTr<-Resp$fTr		# Temperature function of resp
  Act<-Resp$Act  # Activity
  Vel<-Resp$Vel
  Rmax_gO2_g_d<-RA*W^RB 		# max Resp in g O2/g fish /d
  R_gO2_gd<-Rmax_gO2_g_d*fTr # max Resp in g O2/g fish /d
  R_jgd<-R_gO2_gd*(Qox)*Act  		# Resp in j per g fish / d

  # Growth
  G_jgd<-(C_jgd-(R_jgd+SDA_jgd+F_jgd+U_jgd))  # Growth J fish /g fish /day
  G_ggd<-G_jgd/Epred  							# Growth in g fish /g fish /day
  # add a term to allocate lipids to Epred<-fn()
  # add a term for the remainder to be G
  # k as a function of temperature ? rather than constant

  # return(G)
  return(list(
    J_per_gd=data.frame(
      RFR=RFR,
      Eprey=Eprey,
      Epred=Epred,
      TempC=data$TempC,
      W=data$W,
      fTc=fTc,
      fTr=fTr,
      Act=Act,
      Vel=Vel,
      G_jgd=G_jgd,
      C_jgd=C_jgd,
      Cmax_jgd=Cmax_jgd,
      R_jgd=R_jgd,
      R_Act_jgd=R_jgd*Act,
      SDA_jgd=SDA_jgd,
      F_jgd=F_jgd,
      U_jgd=U_jgd,
      def="Joules of prey per gram of pred per day"),
    gFish_per_gFish_d=data.frame(
      RFR=RFR,
      Eprey=Eprey,
      Epred=Epred,
      TempC=data$TempC,
      W=data$W,
      fTc=fTc,
      fTr=fTr,
      Act=Act,
      Vel=Vel,
      G_ggd=G_jgd/Epred,
      C_ggd=C_jgd/Epred,
      Cmax_ggd=Cmax_jgd/Epred,
      R_ggd=R_jgd/Epred,
      SDA_ggd=SDA_jgd/Epred,
      F_ggd=F_jgd/Epred,
      U_ggd=U_jgd/Epred,
      def="gram of pred per gram of pred per day"),
    gPrey_per_gFish_d=data.frame(
      RFR=RFR,
      Eprey=Eprey,
      Epred=Epred,
      TempC=data$TempC,
      W=data$W,
      fTc=fTc,
      fTr=fTr,
      Act=Act,
      Vel=Vel,
      G_ggd=G_jgd/Eprey,
      C_ggd=C_jgd/Eprey,
      Cmax_ggd=Cmax_jgd/Eprey,
      R_ggd=R_jgd/Eprey,
      SDA_ggd=SDA_jgd/Eprey,
      F_ggd=F_jgd/Eprey,
      U_ggd=U_jgd/Eprey,
      def="gram of prey per gram of pred per day"))
  )
}









