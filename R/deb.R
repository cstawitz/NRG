#' Fish Dynamic Energy Budget Model
#'
#' Rcode developed by: Kirstin Holsman  
#' kirstin.holsman@noaa.gov
#' This function calculates the dynamic energy budget as outlined by B. Einarsson et al. / Journal of Theoretical Biology 281 (2011) 1–8
#' @param par is a list of parameters used for the model
#' @param input is a list with inputs for the model including a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param Em  maximum energy density (J cm-3)
#' @param Lm  is the max volumetric length (cm)
#' @param v  is the conductance (cm d-1)
#' @param f  denotes the functional food response (dimensionless)
#' @param g (dimensionless) is the energy investment ratio
#' @param L is the volumetric length
#' s is the shape coeff that related volumetric length to actual length often obtained by fitting a weight–length relationship of the type M=(s*length)^3 to
#' @keywords Temperature, scaling, consumption
#' @export bioE
#' @examples
#' plk_par<-data.frame(RFR=1, Qox=13560,Ceq=2,Req=2,Weq=1,Tco=10,Tcm=15,QC=2.6,CA=0.119,CB=-0.46, 
#'                     RA=0.0075,RB=-0.251, QR=2.6,Tro=13,Trm=18,SA=0.125, Am=2,FA=0.15,UA=0.11)
#' ebs_data<-list(W=100,TempC=0:10,Eprey=5539.6,Epred=4184,indgst=0,diet=0,fTCmodel=NA,fTrmodel=NA,velmodel=NA)
#' bioE(data=ebs_data,par=plk_par)
#' fTfun()
#'13.56 J mg-1 Brett & Groves 1979 Qox for Fish 
par<-data.frame(
  k=.4,
  v=0.2,
  JEAm=.23,
  yVE=0.8,
  kj=0.001,
  EpH=3.93,
  gamma=0.20,
  TA=9100,
  Tr=6.5+273,
  Mv=4.4,
  mu=500,
  s=0.161,
  Lm=2.82,
  dv=1,
  pE=39.3,
  pR=10.00-8.33
)
deb<-function(par,input=list(l,TempC,e,L){
  # data
  W<-data$W
  # TempC<-data$TempC

  ### PARAMETERS
    k<-par$k
    v<-par$v
    JEAm<-par$JEAm
    yVE<-par$yVE
    kj<-par$kj
    EpH<-par$EpH
    gamma<-par$gamma
    TA<-par$TA
    Tr<-par$Tr
    Mv<-par$Mv
    mu<-par$mu
    s<-par$s
    Lm<-par$Lm
    dv<-par$dv
    pE<-par$pE
    pR<-par$pR
  
  
    Em<-(UE*JEAm)/v
    g=(v*Mv)/(k*JEAm*yVE)
    # JEAm maximum assimilation rate per surface area, JEAm (mmol d-1 cm-2)
    # yVE yield of structure from reserve in growth, yVE (dimensionless)
    
    
  # Equation 1: 
  E<-Em*Lm^3*e*l^3
  # Equation 2:
  V=(Lm*l)^3
  # Equation 3:
  EH=Em*Lm^3*UH
  # Equation 4:
  ER=Em*Lm^3*UR
  # Equation 5-7:
  dedt<-(v/Lm*l)*(f-e)
  if(l<e) dldt<-(v/(3*Lm))*((e-1)/(e+g))
  if(l>=e) dldt<-0
  # Equation 8:
  if(UH<UpH){
    dURdt<-0
  }else{
    dURdt<-(v/Lm)*(1-k)*e*l^2*((l+g)/(e+g))*KJ*UpH
  }
  
  #Roe maturity 
  derdt<-gamma*(UR-er)*er
  #Arrhenius temperature
  pT<-PTr*exp((TA/Tr)-(TA/TempC))
  #Roe and fat
  
  Wfat<-((E+ER-Er)/pE)
  Wroe<-(Er/pr)
  M<-dv*V+Wfat+Wroe
  #percentage of the roe weight of the total body weight
  R<-100*(Wroe/M)
  #percentage of the fat weight of the total body weight
  F<-100*(Wfat/M)
  
  # shape coeff that links volumetric length to actual length
  length<-(1/s)*L
  
  
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

