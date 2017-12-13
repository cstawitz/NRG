#' Halibut BIOE model
#'
#' This function calculates the weight and tempature specific physiological paramaters given a set of parameters and data inputs. Predicts consumption or growth from a fit or set pvalue using model.type
#' model.type =
#' 1) fit to observed growth by adjusting P-value
#' 2) predict growth from set p-value (in par)
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$W  weight of the fish in grams
#' @keywords Temperature, scaling, consumption
#' @export f.T
#' @examples
#' par<-c(2)
#'   pars<-data.frame(CA=Cparm2[1],CB=Cparm2[2],Tco=exp(Ft_model$par)[1],Tcm=18,Qc=exp(Ft_model$par)[2],sigma=exp(Ft_model$par)[3])
#'	pars2<-data.frame(CA=Cparm[1],CB=Cparm[2],Tco=exp(Ft_model$par)[1],Tcm=18,Qc=exp(Ft_model$par)[2],sigma=exp(Ft_model$par)[3])
#' 	pars3<-data.frame(CA=Cparm[1],CB=Cparm[2],Tco=exp(Ft_model$par)[1],Tcm=18,CQ=exp(Ft_model$par)[2],O2cal=13560)
#'	Halibut_parms<-list(first_choice=newPAR2,second_choice=Cparm2)
#' 	PARMS<-c(pars,Rpar,UA=UA,UA.sigma=UA.sigma,FA=FA,SA=SA, SA.sigma=SA.sigma)
#'	PARMS_USE<-c(pars2,Rpar,UA=UA,UA.sigma=UA.sigma,FA=FA,SA=SA, SA.sigma=SA.sigma)
#' fish_nrg(par=PARMS_USE,data=list(TempC=10,W=100,Eprey=4000,Pvalue=1))

fish_nrg<-function(par,data){
	TempC<-data$TempC
	Pvalue<-data$Pvalue
	W<-data$W
	Eprey<-data$Eprey

	CA<-par$CA
	CB<-par$CB
	Tco<-par$Tco
	Tcm<-par$Tcm
	Qc<-par$Qc
	RA<-par$RA
	RB<-par$RB
	RQ<-par$RQ
	UA<-par$UA
	FA<-par$FA
	SA<-par$SA
	RK4<-par$RK4
	logAm<-par$logAm
	Bact<-par$Bact
	Trl<-par$Trl
	Tro<-par$Tro

	Yc<-log(Qc)*(Tcm-Tco+2)
	Zc<-log(Qc)*(Tcm-Tco)
	Vc<-(Tcm-TempC)/(Tcm-Tco)
	Xc<-((Zc^2)*(1+((1+40/Yc)^0.5))^2)/400
	f.T<-(Vc^Xc)*exp(Xc*(1-Vc))
	Cmax_g_g_d<-CA*(W^CB)  #g_g_d
	C_g_g_d<-Cmax_g_g_d*f.T*Pvalue
	CmaxFt_g_g_d<-Cmax_g_g_d*f.T
	F_g_g_d<-FA*C_g_g_d
	U_g_g_d<-UA*(C_g_g_d-F_g_g_d)

	Act<-vel_fun(par=c(RK4=RK4,logAm=logAm,Bact=Bact,logsigma=log(.2)),data=(list(Tro=Tro,Trl=Trl,actdat=data.frame(TempC=TempC,W=W,speed=NA),fitLL=FALSE)))$Activity
	
	f.Tr<-exp(RQ*TempC)
	Rmax_gO2_g_d<-RA*W^RB

	R_gO2_g_d<-Rmax_gO2_g_d*f.Tr
	R_g_g_d<-R_gO2_g_d*(13560)*(1/Eprey)
	SDA_g_g_d<-SA*(C_g_g_d-F_g_g_d)

	G_g_g_d<-(C_g_g_d-(R_g_g_d*Act+SDA_g_g_d+F_g_g_d+U_g_g_d))/Efish

	return(list(
		W=W,
		Eprey=Eprey,
		TempC=TempC,
		f.T=f.T,
		f.Tr=f.Tr,
		G_g_g_d=G_g_g_d,
		Cmax_g_g_d=Cmax_g_g_d,
		C_g_g_d=C_g_g_d,
		CmaxFt_g_g_d=CmaxFt_g_g_d,
		Cmax_g_g_d=Cmax_g_g_d,
		R_g_g_d=R_g_g_d,
		R_g_g_dACT=R_g_g_d*Act,
		ACT=Act,
		SDA_g_g_d=SDA_g_g_d,
		F_g_g_d=F_g_g_d,
		U_g_g_d=U_g_g_d))
}

#' simple simulation model
#'
#' This function calculates the weight and tempature specific physiological paramaters given a set of parameters and data inputs. Predicts consumption or growth from a fit or set pvalue using model.type
#' model.type =
#' 1) fit to observed growth by adjusting P-value
#' 2) predict growth from set p-value (in par)
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$W  weight of the fish in grams
#' @keywords Temperature, scaling, consumption
#' @export f.T
#' @examples
#' nrg_sim()

nrg_sim<-function(Tdat,Pvalue, preyEnergy,fishEnergy=20050,W0){
	
	nd<-length(Tdat$date)
	model<-data.frame(Cmax=rep(0,nd))
	model$G<-model$allR<-model$R<-model$W<-model$C<-model$Temp<-0
	model$date<-Tdat$date[1]
	d<-1
	#W0<-1000
	
	for(d in 1:(nd)){
		model$Temp[d]<-Tdat$Temp[d]
		model$date[d]<-Tdat$date[d]
		if(d==1){
			tt.dat<-list(TempC=Tdat$Temp[d],W=W0,Eprey=preyEnergy,Pvalue=Pvalue)
			tt<-Halibut_E(par=PARMS_USE,data=tt.dat)
			model$Cmax[d]<-tt$CmaxFt_g_g_d
			model$G[d]<-tt$G_g_g_d*preyEnergy/fishEnergy
			model$W[d]<-W0+model$G[d]*W0
			model$C[d]<-tt$C_g_g_d
			model$R[d]<-tt$R_g_g_d
			model$allR[d]<-tt$R_g_g_d*tt$ACT+tt$SDA_g_g_d
		}else{
			tt.dat<-list(TempC=Tdat$Temp[d],W=model$W[d-1],Eprey=preyEnergy,Pvalue=Pvalue)
			tt<-Halibut_E(par=PARMS_USE,data=tt.dat)
			model$G[d]<-tt$G_g_g_d*preyEnergy/fishEnergy
			model$Cmax[d]<-tt$CmaxFt_g_g_d
			model$W[d]<-model$W[d-1]+model$G[d]*model$W[d-1]
			model$C[d]<-tt$C_g_g_d
			model$R[d]<-tt$R_g_g_d
			model$allR[d]<-tt$R_g_g_d*tt$ACT+tt$SDA_g_g_d
		}
		
	}
	return(model)
}