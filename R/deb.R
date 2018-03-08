#' Fish Dynamic Energy Budget Model
#'
#' Rcode developed by: Kirstin Holsman  
#' kirstin.holsman@noaa.gov
#' This function calculates the dynamic energy budget as outlined by B. Einarsson et al. / Journal of Theoretical Biology 281 (2011) 1–8
#' According to DEB theory, each individual allocates a fixed
#' fraction k of utilized energy from reserves to growth and somatic
#' maintenance. The rest, (1-k), is then allocated to maturity maintenance and reproduction. 
#' @param par is a list of parameters used for the model
#' @param input is a list with inputs for the model including a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param k  Fraction of energy to somatic growth and maintenance
#' @param v  Energy conductance (cm d-1), temperature dependent
#' @param JEAm  Maximum assimilation rate per surface area (mmol d-1 cm-2)
#' @paran yVE   Yield of structure from reserve in growth, yVE (dimensionless)
#' @param kj  Fraction of maturity maintenancec (d-1), temperature dependent
#' @param EpH Maturity energy at puberty (kJ)
#' @param gamma Growth rate of roe (d-1)
#' @param TA Arrhenius temperature (K)
#' @param Tr Reference temperature (K)
#' @param Mv Volume specific structural mass (mmol cm-3)
#' @param muE  Chemical potential (J mmol-3)
#' @param s Shape coefficient that relates volumetric length to actual length; often obtained by fitting a weight–length relationship of the type M=(s*length)^3 
#' @param Lm Maximum structural (volumetric) length (cm)
#' @param dv Density of structural volume (g cm-3)
#' @param pE Energy reserve density (kJ g-1)
#' @param pR Energy density of roe (kJ g-1)
#'
#' @keywords Temperature, scaling, consumption
#' @export bioE
#' @examples
#' plk_par<-data.frame(RFR=1, Qox=13560,Ceq=2,Req=2,Weq=1,Tco=10,Tcm=15,QC=2.6,CA=0.119,CB=-0.46, 
#'                     RA=0.0075,RB=-0.251, QR=2.6,Tro=13,Trm=18,SA=0.125, Am=2,FA=0.15,UA=0.11)
#' ebs_data<-list(W=100,TempC=0:10,Eprey=5539.6,Epred=4184,indgst=0,diet=0,fTCmodel=NA,fTrmodel=NA,velmodel=NA)
#' bioE(data=ebs_data,par=plk_par)
#' fTfun() 
par<-c(
        k=.4,           # Fraction of energy to somatic growth and maintenance
        v=0.02,         # Energy conductance (cm d-1), temperature dependent
        JEAm=.23,       # Maximum assimilation rate per surface area (mmol d-1 cm-2)
        yVE=0.8,        # Yield of structure from reserve growth 
        kj=0.001,       # Fraction of maturity maintenance (d-1), temperature dependent 
        EpH=3.93,       # Maturity energy at puberty (kJ)
        gamma=0.20,     # Growth rate of roe (d-1)
        TA=9100,        # Arrhenius temperature (K)
        Tr=6.5+273.15,  # Reference temperature (K)
        Mv=4.4,         # Volume specific structural mass (mmol cm-3)
        muE=500,        # Chemical potential (J mmol-3)
        s=0.161,        # Shape coefficient that relates volumetric length to actual length
        Lm=2.82,        # Maximum structural (volumetric) length (cm)
        dv=1,           # Density of structural volume (g cm-3)
        pE=39.3,        # Energy reserve density (kJ g-1)
        pR=10.00-8.33   # Energy density of roe (kJ g-1)
      )
# Initial values:
V0    <- 1e-6
state <- c(E = 0.000124/V0, # reserve density, J/cm^3
           V = V0)          # structural volume, cm^3

deb<-function(par,input=list(l,TempC,e,W,L=NA)){
  with(as.list(c(state, par,input)), {
  # W<-input$W
  # if(is.na(input$L)) L<-W^(1/3)  # structural volume
  # 
  # # data
  # W<-data$W
  # TempC<-input$TempC
  # TempK<-TempC+273.15
  # # TempC<-data$TempC
  # 
  # ### PARAMETERS
  #   k<-par$k
  #   v<-par$v
  #   JEAm<-par$JEAm
  #   yVE<-par$yVE
  #   kj<-par$kj
  #   EpH<-par$EpH
  #   gamma<-par$gamma
  #   TA<-par$TA
  #   Tr<-par$Tr
  #   Mv<-par$Mv
  #   muE<-par$mu
  #   s<-par$s
  #   Lm<-par$Lm
  #   dv<-par$dv
  #   pE<-par$pE
  #   pR<-par$pR
  
    #cals
    
    Em<-(muE*JEAm)/v                       # Max ED
    g=(v*Mv)/(k*JEAm*yVE)                  # energy investment ratio.
    Eg<-(muE*Mv)/yVE                       # volume specific cost of structural volume, J/cm^3
    pM<-(k*muE*JEAm)/Lm                    # volume specific somatic maintenance costs, J/(cm^3*d)
    
    
    # structural volume V is the amount of biomass
    # dynamics are such that maintenance is assumed to take precedence over growth. 
    # L=V^(1/3) denote the structural (volu- metric) length of an individual,
    
    # maximum energy density is a bilogical constant characteristic to each spp.
    # for capelin Em = 586670.43 (J g-11) so muE and JEAm and v were selected to match this value
    
    f <- X/(X_K + X)                       # scaled functional response (food limitation)
    
    tfun<-PTr*exp((TA/Tr)-(TA/TempK))    # temperature function
    
    # rates
    dEdt <- (p_Am*f - ec*E) / V^(1/3)*tfun                 # reserve density
    dVdt <- (kap*ec*E*V^(2/3) - p_M*V)/(kap*E + E_G)*tfun  # structural volume
    
       
  #from Kooijman et al. (2008) 
  # Eq 1: Reserve energy is the energy available to the individual. 
    # Its source is food uptake and it is the energy an organism 
    # utilizes for growth and somatic maintenance on one hand, 
    # and maturity, reproduction and maturity maintenance on 
    # the other hand. 
  
    E<-Em*Lm^3*e*l^3                    # Reserve Energy
  
    e = mE/mEm <-Mev/(L^3*JEAm)          # scaled reserve density
    
    mE<-(ME*L^3)/Mv                      # reserve density
    
    # output variable
    length <- V^(1/3)/del_M                                # physical length, cm
  
  # Equation 2:
    V=(Lm*l)^3
  
  # Equation 3:EH denote a maturity energy. It is important 
    # vto note that this variable is abstract and does not 
    # contribute directly to the weight of the fish. Initially, 
    # energy is allocated to this variable, and the maturity 
    # maintenance will be a fraction of this energy, kJEH.
    # When EH exceeds a certain threshold, EHp , the fish is mature and
    # allocation of energy to EH ceases
  
    EH=Em*Lm^3*UH                  # maturity energy
  
  # Equation 4: After puberty has been reached, the energy starts to flow to ER,
    # which is the total energy available for reproduction. 
    # We note that the dynamics of the energy flow to maturity is 
    # the same as that to reproduction. This energy will, in turn, be converted into roe.
  
    ER=Em*Lm^3*UR
  
  # Equation 5-7:f (dimension- less) denotes the functional food response, 
      # g (dimensionless) is the energy investment ratio.
  dedt<-(v/Lm*l)*(f-e)              # total energy available for reproduction
  
  if(l<e) dldt<-(v/(3*Lm))*((e-1)/(e+g))
  if(l>=e) dldt<-0
  
  # Equation 8:
  if(UH<UpH){
    dURdt<-0
  }else{
    dURdt<-(v/Lm)*(1-k)*e*l^2*((l+g)/(e+g))*KJ*UpH
  }
  
  
  #Roe maturity
  er<-Er/(Em*Lm^3) 
  
  derdt<-gamma*(UR-er)*er
  #Arrhenius temperature

  #Roe and fat
  
  # Er as the energy converted from the reproduction energy to eggs. 
  Wfat<-((E+ER-Er)/pE)
  Wroe<-(Er/pr)
  M<-dv*V+Wfat+Wroe
  #percentage of the roe weight of the total body weight
  R<-100*(Wroe/M)
  #percentage of the fat weight of the total body weight
  F<-100*(Wfat/M)

  # shape coeff that links volumetric length to actual length
  length<-(1/s)*L
  
  }
}

