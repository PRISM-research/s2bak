library(mgcv)
library(foreach)
library(doParallel)

#----------default single species SDMs (for sightings only and S2 models), based on GAMs. These can be replaced with other user-defined SDMs
so_sdm<-function(dat,env_names,numer_env)
{
	names_xso=paste("s(",env_names[numer_env],",k=5)",sep='')   
	if(length(env_names[!numer_env])>0)
		names_xso=c(names_xso,env_names[!numer_env])

	return(gam(formula(paste("pa ~",paste(names_xso,collapse='+'))),family = binomial, data=dat,select = TRUE, method="GCV.Cp"))
}
sdm_predict<-function(sdm,dat)
{
	return(predict.gam(sdm,type="response",newdata=dat))
}
S2_sdm<-function(dat,env_names,numer_env)
{
	names_xf=paste("s(",env_names[numer_env],",so,k=5)",sep='')  
	if(length(env_names[!numer_env])>0) 
		names_xf=c(names_xf,env_names[!numer_env]) 

	return(gam(formula(paste("pa ~",paste(names_xf,collapse='+'))),family = binomial, data=dat,select = TRUE, method="GCV.Cp"))

}


#**********INPUT OPTIONS (see readme file for descriptions)
#options
opts<-list(scale_var=T,min_rec=10,npseudoab=10000,mk_sightings=F,mk_pseudoab=T,comp_S2=T,mk_S2=F,mk_projections=T, mk_val=T,ncores=NA,prjt_integer=10000)

nm_dat<-list(so_dat="./so_dat.rds", surv_dat="./surv_dat.rds",env_dat="./env_dat.rds", env_names='env_names.rds', trait_dat="traits.rds", trait_names="traitnames.rds", val_dat="./val_dat.rds", prjt_dat="./prjt_dat.rds", na_insert="na_insert.rds")

fit_files<-list(so_fit="./so_fit.rds",so_prjt="so_prjt.rds",so_val="so_val.rds", S2_fit="./S2_fit.rds", S2_prjt="S2_prjt.rds", S2_val="S2_val.rds", pseudoab="pseudoab.rds")

sdms_functions<-list(so=so_sdm,S2=S2_sdm,so_predict=sdm_predict,S2_predict=sdm_predict) #default functions defined above




S2BaK<-function(o,nd, ff, sdms)
{
#	o<-opts
#	nd=nm_dat
#	ff=fit_files
#	sdms=sdms_functions
	l<-list() #will contain the output variables/data
#--------------------GENERAL FUNCTIONS
	sum1<-function(x)
	{
		return(sum(x,na.rm=T))
	}
	#--truncation values
	truncatev<-function(v,tp)
	{
		if(tp=="zso")#applied to vector of zso_only
		{
			v[v>15]=15
			v[v< -15]=-15
		}
		if(tp=="prob")#applied to vector of probability values
		{
			v[v<.0001]=0.0001
			v[v>.9999]=0.9999
		}
		if(tp=="pred")
		{
			#truncate x-values outside range of fitting set
			for(i in env_names[numer_env])
			{
				m<-max(pred_surv[,i])
				v[v[,i]>m,i]=m
				m<-min(pred_surv[,i])
				v[v[,i]<m,i]=m
			}
			for(i in env_names[!numer_env])
			{
				v[which(!v[,i] %in% unique(pred_surv[,i])),i]=0
			}
	
		}
		return(v)
	}
	conv_integer<-function(dat)
	{
		dat<-as.matrix(dat)
		dat[]<-as.integer(dat*o$prjt_integer)
		dat<-as.data.frame(dat)
		return(dat)
	}
	#---------for so models
	get_so<-function(sp_nm)
	{
		#get rid of overlap - pseudoabsences do not overlap with presences
		sp_l <-so_dat[so_dat[,2] == sp_nm,3]
		pos=which(l$pseudoab[,1] %in% pred_so[sp_l,1])
		pos1<-1:(o$npseudoab+length(pos))
		if(length(pos)>0)
			pos1<-pos1[-pos]
		#rbind the pseudoabsences data to that (-1 gets rid of index, which is no longer needed

		return(rbind(cbind(pos1,l$pseudoab[pos1,-1]), cbind(pos1=sp_l,pred_so[sp_l,-1])))
	}

	#---------for S2 models
	get_surv<-function(sp_nm)
	{
		pa=surv_dat[,sp_nm]
		return(cbind(pos1=1:nrow(pred_surv),pred_surv[,-1],pa))
	}


	get_full<-function(sp_nm)
	{
		return(rbind(get_surv(sp_nm),get_so(sp_nm)))
	}
	ch_so<-function(d,so)
	{
		d$so=so
		return(d)
	}


	#---------for Validation
	auc<-function(pa,x)
	{
		pos<-which(is.na(pa))
		if(length(pos)>0)
		{
			pa=pa[-pos]
			x=x[-pos]
		}	
		pos<-which(is.na(x))
		if(length(pos)>0)
		{
			pa=pa[-pos]
			x=x[-pos]
		}	

		if(sum(pa)>0)
		{
			pa <- pa[order(x, decreasing=TRUE)]
			return(mean(cumsum(pa)/sum(pa)))
		}
		else
		{
			return(NA)
		}
	}
	dev_exp<-function(pa,x,null)
	{
		aic<-sum1(log(1-abs(pa-x)))
		aic_null<-sum1(log(1-abs(pa-null)))
		return((aic_null-aic)/aic_null)
	}
	registerDoParallel()
	NC=max(detectCores()-1,1)
	options(cores=min(NC,o$ncores,na.rm=T))
	
	print("loading/processing data")
	#--------INPUT DATA HERE
	#note: NAs in surv_dat and val_dat can be introduced if some species/locations were not sampled (e.g., locations where only trees sampled). If so, so_fit should likewise be modified. No location or species should be fully NA. 
	#note: species names must be consistent between so_dat, surv_dat, val_dat
		so_dat<-readRDS(nd$so_dat) #presence_only - contains location index (in column 1) and species names (column 2)
		surv_dat<-readRDS(nd$surv_dat) #survey data, contains location index (in column 1), and species pres/abs
		if(o$mk_val==T)
			val_dat<-readRDS(nd$val_dat)

	#note: predictors should have no NAs, or will cause mismatches below. Impute values prior to analysis.
		env_dat<-readRDS(nd$env_dat) 
		env_names=readRDS(nd$env_names)
		#traits
		trait_dat<-readRDS(nd$trait_dat)
		trait_names=readRDS(nd$trait_names)
		if(o$mk_projections==T)
		{
			prjt_dat<-readRDS(nd$prjt_dat)# this is the environmental conditions. Note: env_dat & prjt_dat column names must be consistent
			if(is.vector(prjt_dat)) #indices from env_dat, otherwise input matrix
				prjt_dat=env_dat[prjt_dat,]
		}
		if(!is.na(nd$na_insert))
			na_dat=readRDS(nd$na_insert)

	#-----------PROCESSING INFORMATION

		numer_env<-unlist(lapply(env_dat[1,env_names],is.numeric))
		numer_tr<-unlist(lapply(trait_dat[1,trait_names],is.numeric))

		if(o$scale_var==T)
		{
			if(o$mk_projections==T) #need to scale prjt to original env_dat
			{
				for(i in env_names[numer_env])
				{
					m=mean(env_dat[,i])
					s=sd(env_dat[,i])
					env_dat[,i]=(env_dat[,i]-m)/s
					prjt_dat[,i]=(prjt_dat[,i]-m)/s
				}
			}else{
			env_dat[,env_names[numer_env]]=as.numeric(scale(env_dat[,env_names[numer_env]]))
			}
			trait_dat[,trait_names[numer_tr]]=as.numeric(scale(trait_dat[,trait_names[numer_tr]]))
		}

		#making polynomial predictor variable names: note can replace with more complex models
		names_xbe=c(env_names[numer_env],paste("I(",env_names[numer_env],"^2)",sep=''))
		names_xbt=c(trait_names[numer_tr],paste("I(",trait_names[numer_tr],"^2)",sep=''))
		#add in factors or characters
		if(length(env_names[!numer_env])>0)
			names_xbe=c(names_xbe,env_names[!numer_env])
		if(length(trait_names[!numer_tr])>0)
			names_xbt=c(names_xbt,trait_names[!numer_tr])

		#filtering steps - min size, and surveyed species restricted to ones in sightings data
		so_datsp=tapply(so_dat[,1],so_dat[,2],length)
		sp2=names(so_datsp[so_datsp>o$min_rec])
		sp2=sp2[!is.na(sp2)] 
		sp2=sp2[order(sp2)]
		sp1=colnames(surv_dat)
		sp=sp1[which( sp1 %in% sp2)] #only analyzing surveyed species also observed in presence-only
		l$sp2=sp2
		l$sp=sp
		if(o$mk_val==T)
		{
			sp_val=colnames(val_dat)
			sp_val=sp_val[which( sp_val %in% sp2)]
			l$sp_val=sp_val
		}

	#-------------CREATE DATAFRAMES NEEDED FOR FITTING

		#make survey data
		pred_surv=env_dat[surv_dat[,1],]
		pred_surv$so=0 #not presence only
		#make presence only
		pos=unique(so_dat[,1]) #gbif indices from original list
		index2=rep(NA,nrow(env_dat))
		index2[pos]=1:length(pos) #match the right row number to the new matrix pred_so
		pred_so=env_dat[pos,]
		so_dat[,3]=index2[so_dat[,1]]
		pred_so$so=1 #presence only
		pred_so$pa=1 #present/not absent

		if(o$mk_val==T)
		{
			pred_val=env_dat[val_dat[,1],]
			pred_val$so=0 #not sightings only
		}



	if (o$mk_pseudoab == T){
	#Choose 10000 pseudoabsences from terrestrial cells - remove survey sites
		pos=(1:nrow(env_dat))[-pred_surv$index]
		l$pseudoab <- env_dat[sample(pos, o$npseudoab +max(so_datsp,na.rm=T), replace=F),] # add extra so that can have absences, rather than reference sites
		l$pseudoab$so <- 1 #sightings only
		l$pseudoab$pa <- 0 #absent


	}else{
		l$pseudoab=readRDS(ff$pseudoab)
	}


	l$so=list()
	
	tmpA=0
	tmpB=0
	if(o$mk_val==T)
		tmpA=nrow(pred_val)
	if(o$mk_projections==T)
		tmpB=nrow(prjt_dat)
	if (o$mk_sightings==T)
	{
	print("making sightings only")

		tmp<-as.data.frame(matrix(NA,nrow=nrow(pred_surv)+tmpA+tmpB,ncol=length(sp2)))
		names(tmp)=sp2
	#-----------run SDM for each species, fill dataframes
		print(Sys.time())
		 
		tmp[]<-foreach(i = sp2) %dopar%
		{
			print(i)
			tmp_na=NULL
			if(!is.na(nd$na_insert)) #allow for NAs, and removal of species from consideration in surveyed locations
			{
				pos=which(na_dat[,2] == i)
				if(length(pos)>0) #the species exists
					tmp_na=na_dat[pos,1] #provides the indices, corresponding to env_dat,and column1 of data
			}


			tmp_sdm<-sdms$so(get_so(i),env_names,numer_env)
			tmp1<-sdms$so_predict(tmp_sdm,pred_surv)
			tmp2=NULL
			tmp3=NULL
			if(o$mk_val==T)
				tmp2<-sdms$so_predict(tmp_sdm,pred_val)
			if(o$mk_projections==T)
				tmp3<-sdms$so_predict(tmp_sdm,prjt_dat)

			if(!is.null(tmp_na))
			{
				tmp1[which(pred_surv[,1] %in% tmp_na)]=NA
				if(o$mk_val==T)
					tmp2[which(pred_val[,1] %in% tmp_na)]=NA
			}

			return(c(tmp1,tmp2,tmp3))
		}
		ep=nrow(pred_surv)
		l$so$fit=tmp[1:ep,]
		if(o$mk_val==T)
		{
			l$so$val=tmp[(ep+1):(ep+nrow(pred_val)),]
			ep=ep+nrow(pred_val)
		}
		if(o$mk_projections==T)
		{
			l$so$prjt=tmp[(ep+1):nrow(tmp),]
			if(!is.na(o$prjt_integer))
			{
				l$so$prjt=conv_integer(l$so$prjt)
			}
		}

		print(Sys.time())

#	saveRDS(l$so,"lso.rds")
	}else{

#		l$so=readRDS("lso.rds")
		l$so$fit=readRDS(ff$so_fit)
		if(o$mk_projections==T)
			l$so$prjt<-readRDS(ff$so_prjt)
		if(o$mk_val==T)
			l$so$val=readRDS(ff$so_val)
	}
	
	if(o$comp_S2==T)
	{
		l$S2=list()
		if(o$mk_S2==T)
		{
			print("making S2")
			#-------------data frames filled
			tmp<-as.data.frame(matrix(NA,nrow=nrow(pred_surv)+tmpA+tmpB,ncol=length(sp)))
			names(tmp)=sp
			print(Sys.time())
			tmp[]<-foreach(i = sp) %dopar%
			{
				print(i)
				tmp_sdm<-sdms$S2(get_full(i),env_names,numer_env)
				tmp1<-sdms$S2_predict(tmp_sdm,pred_surv)
				tmp2=NULL
				tmp3=NULL
				if(o$mk_val==T)
					tmp2<-sdms$S2_predict(tmp_sdm,pred_val)
				if(o$mk_projections==T)
					tmp3<-sdms$S2_predict(tmp_sdm,cbind(prjt_dat,so=0))
				return(c(tmp1,tmp2,tmp3))
			}
			ep=nrow(pred_surv)
			l$S2$fit=tmp[1:ep,]
			if(o$mk_val==T)
			{
				l$S2$val=tmp[(ep+1):(ep+nrow(pred_val)),]
				ep=ep+nrow(pred_val)
			}
			if(o$mk_projections==T)
			{
				l$S2$prjt=tmp[(ep+1):nrow(tmp),]
				if(!is.na(o$prjt_integer))
				{
					l$S2$prjt=conv_integer(l$S2$prjt)
				}
			}

			print(Sys.time())
#			saveRDS(l$S2,"lS2.rds")

		}else{

#			l$S2=readRDS("lS2.rds")
			l$S2$fit=readRDS(ff$S2_fit)
			if(o$mk_projections==T)
				l$S2$prjt<-readRDS(ff$S2_prjt)
			if(o$mk_val==T)
				l$S2$val=readRDS(ff$S2_val)

		}

	}

	#--------CALCULATE BIAS:make dataframes based on summed expectations vs observations
	print("calculating bias kernels")
	l$bak=list()
	fit_l=pred_surv
	fit_l$pa=apply(surv_dat[,sp],1,sum1)#get the summed occurrences
	fit_l$so_only=apply(l$so$fit,1,sum1)#get the summed expected sightings
	fit_l$lr=log((fit_l$pa+1)/(fit_l$so_only+1))

	pos=tapply(1:nrow(trait_dat),trait_dat[,1],function(x){return(x[1])})
	fit_sp=data.frame(sp=sp)
	fit_sp=cbind(fit_sp,trait_dat[pos[sp],-1]) 
	fit_sp$sp=as.character(fit_sp$sp)
	#needs to be for traits - set missing species to the mean -0
	fit_sp[is.na(fit_sp)]=0
	fit_sp$pa=apply(surv_dat[,fit_sp[,1]],2,sum1)#get the summed occurrences
	fit_sp$so_only=apply(l$so$fit[,fit_sp[,1]],2,sum1)#get the summed expected sightings
	fit_sp$lr=log((fit_sp$pa+1)/(fit_sp$so_only+1))
	#--------get bias kernels
	
	l$bak$bias_l<-gam(formula(paste("lr ~",paste(names_xbe,collapse='+'))), data=fit_l,select = TRUE, method="GCV.Cp")
	l$bak$bias_sp<-gam(formula(paste("lr ~",paste(names_xbt,collapse='+'))), data=fit_sp,select = TRUE, method="GCV.Cp")

	print("making bak")
	#Can calculate once across species and re-use for all instances
	fit_sp2=data.frame(sp=sp2)
	fit_sp2=cbind(fit_sp2,trait_dat[pos[sp2],-1]) 
	fit_sp2$sp=as.character(fit_sp2$sp)
	fit_sp2[is.na(fit_sp2)]=0
	for(i in trait_names[!numer_tr]) #any non-fit factor levels need to set to zero
	{
		fit_sp2[which(!fit_sp2[,i] %in% unique(fit_sp[,i])),i]=0	
	}

	fit_sp2$scale_so<-predict.gam(l$bak$bias_sp,newdata=fit_sp2)


	#---------apply bias kernel and so_only model to fit final model
	mk_bak<-function(dat,so_only,pred,sp_nm) #pass species info, so_only predictions, and env pred
	{
		nrs=nrow(dat)
		bak<-as.data.frame(matrix(nrow=nrs*length(sp2),ncol=6)) #need only location, species name, so,zso,and observations, then add bias adjustments, truncate
		names(bak)=c('loc','sp','pa','zso_only',"scale_so_l","scale_so_sp")

		bak$zso_only=truncatev(as.vector(as.matrix(-log((1-so_only)/so_only))),"zso")
		bak$sp=rep(sp2,each=nrs)
		bak$loc=rep(1:nrs,length(sp2))
		bak$scale_so_l=predict.gam(l$bak$bias_l,newdata=pred) 
		bak$scale_so_sp=rep(fit_sp2$scale_so,each=nrs) 
		bak$pa=0
		#do it multistep, since sp and sp2 may not be in order
		osp<-tapply(1:nrow(bak),bak$sp,function(x){return(x)}) #this puts sp2 in order
		nsp2=names(osp) #this tracks ordered sp2 
		pos=which(nsp2 %in% sp_nm) #This fills in pa - i.e., observed species
		x<-as.vector(as.matrix(dat[,nsp2[pos]]))
		bak$pa[unname(unlist(osp[pos]))]=x
		return(bak)
	}
	fill_S2<-function(S2_pred, dat)#dat is either fit or val, S2 are predictions (which is done for all sp)
	{
		#do it multistep, since sp and sp2 may not be in order
		osp<-tapply(1:nrow(dat),dat$sp,function(x){return(x)}) #this puts sp2 in order
		nsp2=names(osp) #this tracks ordered sp2 
		pos=which(nsp2 %in% sp) #This fills in S2, which is calculated for all sp species, regardless of whether it is fit or val 
		x<-as.vector(as.matrix(S2_pred[,nsp2[pos]]))
		tmpS2=rep(NA,nrow(dat))
		tmpS2[unname(unlist(osp[pos]))]=x
		return(tmpS2)
	}
	 
	#	S2 based on sp (fitted species)

	l$eval=list() 
	l$fit=mk_bak(surv_dat,l$so$fit,pred_surv,sp)
	l$bak$bias_adj<-glm(pa~zso_only+scale_so_sp+scale_so_l, family="binomial", data=l$fit)
	l$bak$coef=l$bak$bias_adj$coefficients
	l$fit$bak=truncatev(predict.glm(l$bak$bias_adj,type="response",newdata=l$fit,na.action=na.pass),"prob")
	
	null=mean(l$fit$pa,na.rm=T) #the expectation base on fitting data, if no model were used
	print("Fitting")
	print(paste("num species:",length(unique(l$fit$sp))))
	l$eval$fit_auc_bak=auc(l$fit$pa,l$fit$bak)
	l$eval$fit_dev_bak=dev_exp(l$fit$pa,l$fit$bak, null)
	print(paste("fit_bak: auc: ",l$eval$fit_auc_bak,"dev_exp: ",l$eval$fit_dev_bak))

	if(o$comp_S2==T)
	{
		l$fit$S2=truncatev(fill_S2(l$S2$fit,l$fit),"prob") #can then compare this against pa, for non NA rows
		nullS2=mean(as.matrix(surv_dat[,sp]),na.rm=T) #only applied to fitted species
		pos=which(!is.na(l$fit$S2))
		#included for comparison: bak applied only for surveyed species 
		l$eval$fit_auc_bak_surv=auc(l$fit$pa[pos],l$fit$bak[pos])
		l$eval$fit_dev_bak_surv=dev_exp(l$fit$pa[pos],l$fit$bak[pos],nullS2)
		print(paste("num surveyed species:",length(unique(l$fit$sp[pos]))))
		print(paste("fit_bak_surv: auc: ",l$eval$fit_auc_bak_surv,"dev_exp: ",l$eval$fit_dev_bak_surv))


		l$eval$fit_auc_S2=auc(l$fit$pa[pos],l$fit$S2[pos])
		l$eval$fit_dev_S2=dev_exp(l$fit$pa[pos],l$fit$S2[pos],nullS2)
		print(paste("fit_S2: auc: ",l$eval$fit_auc_S2,"dev_exp: ",l$eval$fit_dev_S2))

	}
	if(o$mk_val==T)
	{
	print("")
	print("Validation")
		l$val=mk_bak(val_dat,l$so$val,truncatev(pred_val,"pred"),sp_val)
		l$val$bak=truncatev(predict.glm(l$bak$bias_adj,type="response",newdata=l$val,na.action=na.pass),"prob")
		l$eval$val_auc_bak=auc(l$val$pa,l$val$bak)
		l$eval$val_dev_bak=dev_exp(l$val$pa,l$val$bak, null)
		print(paste("num species:",length(unique(l$fit$sp))))
		print(paste("val_bak: auc: ",l$eval$val_auc_bak,"dev_exp: ",l$eval$val_dev_bak))
		if(o$comp_S2==T)
		{
			l$val$S2=truncatev(fill_S2(l$S2$val,l$val),"prob") #can then compare this against pa, for non NA rows
			pos1=which(!is.na(l$val$S2))

		#included for comparison: bak applied only for surveyed species 
			l$eval$val_auc_bak_surv=auc(l$val$pa[pos1],l$val$bak[pos1])
			l$eval$val_dev_bak_surv=dev_exp(l$val$pa[pos1],l$val$bak[pos1],nullS2)
			print(paste("num surveyed species:",length(unique(l$val$sp[pos1]))))
			print(paste("val_bak_surv: auc: ",l$eval$val_auc_bak_surv,"dev_exp: ",l$eval$val_dev_bak_surv))

			l$eval$val_auc_S2=auc(l$val$pa[pos1],l$val$S2[pos1])
			l$eval$val_dev_S2=dev_exp(l$val$pa[pos1],l$val$S2[pos1], nullS2)
			print(paste("val_S2: auc: ",l$eval$val_auc_S2,"dev_exp: ",l$eval$val_dev_S2))

		}
	}


	#for forecasting
	if(o$mk_projections==T)
	{
		print("Forecasting")
		sc=l$bak$bias_adj$coefficients
		scale_so_l=predict.gam(l$bak$bias_l,newdata=truncatev(prjt_dat,"pred")) #apply to all locations
		if(!is.na(o$prjt_integer)) 
		{
		#assumes that sightings only projections are scaled - if zero, gives truncated value
			zso_only=truncatev(as.vector(-log((o$prjt_integer-l$so$prjt)/l$so$prjt)),"zso")
		}else{
			zso_only=truncatev(as.vector(-log((1-l$so$prjt)/l$so$prjt)),"zso")
		}
		z=sc[1]+sc[2]*zso_only+sc[3]*rep(fit_sp2$scale_so,each=nrow(prjt_dat))+sc[4]*scale_so_l
		l$bak$prjt=truncatev(1/(1+exp(-z)),"prob")
		if(!is.na(o$prjt_integer))
		{
			l$bak$prjt=conv_integer(l$bak$prjt)
		}
	}
	#returns S2, so_only, bak, validation, fitting, projection, evaluation - see readme file.
	#saveRDS(l,"all.rds")
	return(l) 

}	
#---------------Running S2BaK

#a<-S2BaK(opts,nm_dat,fit_files,sdms_functions)


