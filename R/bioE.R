
#' Fish Bioenergetics model
#'
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#' This function calculates the weight and tempature specific physiological rates given a set of parameters and data inputs. Predicts consumption or growth from a fit or set pvalue using model.type
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#' model.type =
#' 1) fit to observed growth by adjusting P-value
#' 2) predict growth from set p-value (in par)
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$W  weight of the fish in grams
#' @keywords Temperature, scaling, consumption
#' @export fTc
#' @examples
#' plk_par<-data.frame(RFR=1, Qox=13560,Ceq=2,Req=2,Weq=1,Tco=10,Tcm=15,QC=2.6,CA=0.119,CB=-0.46, 
#'                     RA=0.0075,RB=-0.251, QR=2.6,Tro=13,Trm=18,SA=0.125, Am=2,FA=0.15,UA=0.11)
#' ebs_data<-list(W=100,TempC=0:10,Eprey=5539.6,Epred=4184,indgst=0,diet=0)
#' bioE(data=ebs_data,par=plk_par)

#' fTfun()
#13.56 J mg-1 Brett & Groves 1979 Qox for Fish 
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


#' Respiration function; Fish Bioenergetics model
#' This function calculates the temperature and weight specific metabolic demand for an individual fish
#' 
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#' This function calculates the weight and tempature specific physiological rates given a set of parameters and data inputs. 
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#' 
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$W  weight of the fish in grams
#' @keywords Temperature, scaling, Respiration
#' @export 
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

#' Waste function; Fish Bioenergetics model
#' This function calculates the temperature and weight specific estimates of lost energey (waste) for an individual fish
#' 
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#' This function calculates the weight and tempature specific physiological paramaters given a set of parameters and data inputs. 
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#' 
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$W  weight of the fish in grams
#' @keywords Temperature, scaling, Respiration
#' @export 
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

#' Consumption function; Fish Bioenergetics model
#' This function calculates the temperature and weight specific consumption rate for an individual fish
#' 
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#' This function calculates the weight and tempature specific physiological paramaters given a set of parameters and data inputs.
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#' 
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param Ceq is a par list object (par$Ceq) that specifies the type of consumption equation to use: 0 = specify function using data$fTCmodel, 1= exponential (type 1 from Fish Bioenergetics), 2 = type 2, 3 = type 3
#' @keywords Temperature, scaling, Respiration
#' @export 
#' @examples
#' c_data<-ebs_data
#' c_data$fTCmodel<-function(TempC){-.5*TempC}
#' fTC_fun(par=plk_par,data=c_data)


fTC_fun<-function(par,data){
  TempC<-data$TempC
  Ceq<-par$Ceq
  if(Ceq==0){
    fTc<-data$fTCmodel(TempC)
  }else if  (Ceq==1){
    fTc<-exp(par$QC*TempC)
  }else if (Ceq==2){
    Tco<-par$Tco
    Tcm<-par$Tcm
    Yc<-log(par$QC)*(Tcm-Tco+2)
    Zc<-log(par$QC)*(Tcm-Tco)
    Vc<-(Tcm-TempC)/(Tcm-Tco)
    Xc<-((Zc^2)*(1+((1+40/Yc)^0.5))^2)/400
    fTc<-(Vc^Xc)*exp(Xc*(1-Vc))
  }else if (Ceq==3){
    CK1<-par$CK1
    CK4<-par$CK4
    Tcl<-par$Tcl
    Tco<-par$Tco
    Tcm<-par$Tcm
    G2<-(1/(Tcl-Tcm))*log((0.98*(1-CK4))/(CK4*0.02))
    L2<-exp(G2*(Tcl-TempC))
    Kb<-(CK4*L2)/(1+CK4*(L2-1))
    G1<-(1/(Tco-par$QC))*log((0.98*(1-CK1))/(CK1*0.02))
    L1<-exp(G1*(TempC-par$QC))
    Ka<-(CK1*L1)/(1+CK1*(L1-1))
    fTc<-Ka*Kb
  }else{
    message("ERROR: Ceq does not match criteria [0,3]")
  }
  return(fTc=fTc)
}



