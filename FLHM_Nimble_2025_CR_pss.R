
#Nimble
WGBASTCode<-nimbleCode({ 
  
  #here i index denotes current year
  pmat[1,1]<-0                                                                                                      
  pimm[1,1]<-1
  
  for(j in 2:6){
    pmat[1,j]<-LW[1,(j-1)]*pimm[1,(j-1)]      
    pimm[1,j]<-(1-LW[1,(j-1)])*pimm[1,(j-1)]
    
  }
  
  for(i in 2:(m)){
    pmat[i,1]<-0
    pimm[i,1]<-1
    
    for(j in 2:6){
      pmat[i,j]<-LW[i,(j-1)]*pimm[(i-1),(j-1)]     
      pimm[i,j]<-(1-LW[i,(j-1)])*pimm[(i-1),(j-1)]
      
    }
  }
  
  for(s in 1:stocks){
    #here i index denotes smolt year
    
    simm[1,1,s]<-exp(-(11*MpsW[1,s]/Tstep))*exp(-(MpsW[1,s]*sealMort[1,1,1]/Tstep))     #note AU 1 sealMort for all change later
    smat[1,1,s]<-1    #not used anywhere
    
    for(j in 2:6){
      simm[1,j,s]<-exp(-(12*MW/Tstep))
      smat[1,j,s]<-exp(-(3*MW/Tstep))*exp(-(2*MW*sealMort[1,j,AU[s]]/Tstep))*p.ladder[(1+j-1),s]*surv_migr[(1+j-1),s]
    }
    
    for(i in 2:(m)){
      
      simm[i,1,s]<-exp(-(11*MpsW[i,s]/Tstep))*exp(-(MpsW[i,s]*sealMort[i,1,1]/Tstep))   #note AU 1 sealMort for all change later
      smat[i,1,s]<-1    #not used anywhere
      
      for(j in 2:6){
        simm[i,j,s]<-exp(-(12*MW/Tstep))
        smat[i,j,s]<-exp(-(3*MW/Tstep))*exp(-(2*MW*sealMort[i,j,AU[s]]/Tstep))*p.ladder[(i+j-1),s]*surv_migr[(i+j-1),s]
      }
    }
    
    #EPR calculation
    
    for(i in 1:5){
      
      EPR[i,s]<-0
      EPR_M74[i,s]<-EPR[i,s]*(1-M74[i,s])
      
      z[i,s]<-(slope[s]*EPR_M74[i,s])/(4+slope[s]*EPR_M74[i,s])	  
    }
    
    for(i in 6:(m)){    
      
      EPR[i,s]<-pmat[i,2]*prop_fem[(i-1),1,s]*fec[1]*simm[(i-1),1,s]*smat[(i-1),2,s]+pmat[i,3]*prop_fem[(i-2),2,s]*fec[2]*simm[(i-2),1,s]*simm[(i-1),2,s]*smat[(i-2),3,s]+pmat[i,4]*prop_fem[(i-3),3,s]*fec[3]*simm[(i-3),1,s]*simm[(i-2),2,s]*simm[(i-1),3,s]*smat[(i-3),4,s]+pmat[i,5]*prop_fem[(i-4),4,s]*fec[4]*simm[(i-4),1,s]*simm[(i-3),2,s]*simm[(i-2),3,s]*simm[(i-1),4,s]*smat[(i-4),5,s]+pmat[i,6]*prop_fem[(i-5),5,s]*fec[5]*simm[(i-5),1,s]*simm[(i-4),2,s]*simm[(i-3),3,s]*simm[(i-2),4,s]*simm[(i-1),5,s]*smat[(i-5),6,s]     #5SW
      
      EPR_M74[i,s]<-EPR[i,s]*(1-M74[i,s])
      
      z[i,s]<-(slope[s]*EPR_M74[i,s])/(4+slope[s]*EPR_M74[i,s])  
    }
  }
  
  # Northern stocks
  # prior: mean(slope)=0.04, sd_slope=0.017
  mu_a[1]~dnorm(-2.784,2.388)
  sd_a[1]~dlnorm(-0.2653,3.1529)
  tau_a[1]<-1/(sd_a[1]*sd_a[1])
  
  # Southern stocks
  # prior: mean(slope)=0.04, sd_slope=0.017
  mu_a[2]~dnorm(-2.784,2.388)
  sd_a[2]~dlnorm(-0.2653,3.1529)
  tau_a[2]<-1/(sd_a[2]*sd_a[2])
  
  # From eggs to smolts	
  # Total egg production (Eggstot) is a function of sex ratio, fecundity and spawner abundance		
  for (s in 1:stocks){
    
    a_slope[s]~dnorm(mu_a[SR_unit[s]],tau_a[SR_unit[s]])
    logit(slope[s])<-a_slope[s]
    alphaSR[s] <- 1/slope[s]				
    K[s]~dlnorm(M_R0[s],tau_R0[s])     #M_R0 and tau_R0 come from the data file K_prior
    betaSR[s] <- 1/K[s]
    for(i in 1:(6+e_delay[s]-1)){
      smoltPred[i,s]<-0
    }
    
    for(i in 1:5){ 
      
      NspWtot[i,s]<-0
      NrW_msw[i,s]<-0
      NrW_tot[i,s]<-0
      NladderW_tot[i,s]<-0
      NrRsp_msw[i,s]<-0
      NrRsp_tot[i,s]<-0
      NrAll_msw[i,s]<-0
      NrAll_tot[i,s]<-0  
      
      PropWgrilse[i,s]<-0
      PropWMSW[i,s]<-0
      
      mean_sp[i,s]<-0
      var_sp[i,s]<-0
      cv_sp[i,s]<-0
      tau_sp[i,s]<-0
      mu_sp[i,s]<-0
      sp_countX[i,s]<-0
      
      mean_la[i,s]<-0
      var_la[i,s]<-0
      cv_la[i,s]<-0
      tau_la[i,s]<-0
      mu_la[i,s]<-0
      
      eta_star[i,s]<-0
      alpha_msw[i,s]<-1
      beta_msw[i,s]<-1
      NrWmsw[i,s] <-0
      NrWtot[i,s]<-0 
      probMSW[i,s] <-0
      
      Eggstot[i,s]<-0
      Eggstot_M74[i,s]<-0
      
      M74[i,s]~dbeta(M74_alpha[i,s], M74_beta[i,s])   #M74 alpha and beta from M74 data file
      error_SR[i,s] ~dnorm(mu_SR,tau_SR)
      
      
    }
  }
  
  for(s in 1:2){     #Torne, Simo
    for (i in 6:(m)){
      
      Eggstot[i,s] <- (NspW[(i-1),2,s]+NspRsp[(i-1),2,s])*prop_fem[(i-1),1,s] * fec[1] +(NspW[(i-2),3,s]+NspRsp[(i-2),3,s])* prop_fem[(i-2),2,s] * fec[2] +(NspW[(i-3),4,s]+NspRsp[(i-3),4,s])* prop_fem[(i-3),3,s] * fec[3] + (NspW[(i-4),5,s]+NspRsp[(i-4),5,s]) * prop_fem[(i-4),4,s] * fec[4] + (NspW[(i-5),6,s]+NspRsp[(i-5),6,s])* prop_fem[(i-5),5,s] * fec[5]
      
      Eggstot_M74[i,s]<- Eggstot[i,s] * (1-M74[i,s])	  #Eggs surviving after M74
      
      # Beverton Holt SR model for predicted smolt abundance
      #In AU 1-3 it takes 4 years from egg to smolt (1st BH smolts in yr 10)		
      # In AU 4 it takes 3 years from egg to smolt (1st BH smolts in yr 9)
      
      smoltPred[(i+e_delay[s]),s] <-log(Eggstot_M74[i,s]/(alphaSR[s]+betaSR[s]*Eggstot_M74[i,s]))+error_SR[i,s]
      
      M74[i,s]~dbeta(M74_alpha[i,s], M74_beta[i,s])   #M74 alpha and beta from M74 data file
      error_SR[i,s] ~dnorm(mu_SR,tau_SR)
      
      NspWtot[i,s] <- NspW[(i-1), 2,s]+NspW[(i-2), 3,s] +NspW[(i-3), 4,s] +NspW[(i-4), 5,s] + NspW[(i-5),6,s]             +NspRsp[(i-1), 2,s]+NspRsp[(i-2), 3,s] +NspRsp[(i-3), 4,s] + NspRsp[(i-4), 5,s] +NspRsp[(i-5), 6,s] 
      
      NrW_msw[i,s] <- NrW[(i-2),3,s] +NrW[(i-3),4,s] +NrW[(i-4),5,s] +NrW[(i-5),6,s]
      NrW_tot[i,s] <- NrW[(i-1),2,s]+NrW_msw[i,s] 
      NladderW_tot[i,s]<- NladderW[(i-1),2,s]+NladderW[(i-2),3,s]+NladderW[(i-3),4,s]+NladderW[(i-4),5,s]+NladderW[(i-5),6,s] 
      
      NrRsp_msw[i,s]<-NrRsp[(i-2),3,s]+ NrRsp[(i-3),4 ,s]+NrRsp[(i-4),5,s]+NrRsp[(i-5), 6,s]  
      NrRsp_tot[i,s]<-NrRsp[(i-1),2,s]+NrRsp_msw[i,s]
      
      NrAll_msw[i,s] <- NrW_msw[i,s]+NrRsp_msw[i,s] 
      NrAll_tot[i,s] <- NrW_tot[i,s]+NrRsp_tot[i,s]  #NrWtot1 winbugs (1000s)
      
      PropWgrilse[i,s]<-NrW[(i-1), 2,s] /(NrW[(i-1), 2,s]+NrRsp[(i-1), 2,s])   #propn of grilse that are wild
      PropWMSW[i,s]<-NrW_msw[i,s]/NrAll_msw[i,s]   #propn MSW wild
      
      NrWtot[i,s] <- round((NrW_tot[i,s]+NrRsp_tot[i,s])*1000)+1  
      #doesn't make a diff for Ume here if NrW or NladderW used to calc below ratio as p.ladder same all ages 
      probMSW[i,s] <- NrAll_msw[i,s]/NrAll_tot[i,s] 
      
      
      
    }
    for(i in smolt_year[s]:(m+2)){
      SmoltWobs[i,s]~dlnorm(smoltPred[i,s],tau_SmoltW[i,s])
      R0[i,s]<-(K[s]*(EPR_M74[(i-e_delay[s]),s]-alphaSR[s]))/EPR_M74[(i-e_delay[s]),s] 
    }
  }
  
  for(s in 3:stocks){
    for (i in 6:(m)){
      
      Eggstot[i,s] <- NspW[(i-1),2,s]* prop_fem[(i-1),1,s] * fec[1] + NspW[(i-2),3,s]* prop_fem[(i-2),2,s] * fec[2] + NspW[(i-3),4,s]* prop_fem[(i-3),3,s] * fec[3] + NspW[(i-4),5,s] * prop_fem[(i-4),4,s] * fec[4] + NspW[(i-5),6,s]* prop_fem[(i-5),5,s] * fec[5]
      
      Eggstot_M74[i,s]<- Eggstot[i,s] * (1-M74[i,s])	  #Eggs surviving after M74
      
      # Beverton Holt SR model for predicted smolt abundance
      #In AU 1-3 it takes 4 years from egg to smolt (1st BH smolts in yr 10)		
      # In AU 4 it takes 3 years from egg to smolt (1st BH smolts in yr 9)
      
      smoltPred[(i+e_delay[s]),s] <-log(Eggstot_M74[i,s]/(alphaSR[s]+betaSR[s]*Eggstot_M74[i,s]))+error_SR[i,s]   
      M74[i,s]~dbeta(M74_alpha[i,s], M74_beta[i,s])   #M74 alpha and beta from M74 data file
      error_SR[i,s] ~dnorm(mu_SR,tau_SR)
      
      NspWtot[i,s] <- NspW[(i-1), 2,s] +NspW[(i-2), 3,s] +NspW[(i-3), 4,s] +NspW[(i-4), 5,s] + NspW[(i-5),6,s] 
      
      NrW_msw[i,s] <- NrW[(i-2),3,s] +NrW[(i-3),4,s] +NrW[(i-4),5,s] +NrW[(i-5),6,s]
      NrW_tot[i,s] <- NrW[(i-1),2,s]+NrW_msw[i,s]  #NrWtot1 winbugs (1000s)
      NladderW_tot[i,s]<- NladderW[(i-1),2,s]+NladderW[(i-2),3,s]+NladderW[(i-3),4,s]+NladderW[(i-4),5,s]+NladderW[(i-5),6,s] 
      
      
      NrRsp_msw[i,s]<-0  
      NrRsp_tot[i,s]<-0
      
      NrAll_msw[i,s] <- NrW_msw[i,s]+NrRsp_msw[i,s] 
      NrAll_tot[i,s] <- NrW_tot[i,s]+NrRsp_tot[i,s]  
      
      NrWtot[i,s] <- round((NrW_tot[i,s]+NrRsp_tot[i,s])*1000)+1    
      #doesn't make a diff for Ume here if NrW or NladderW used to calc below ratio as p.ladder same all ages 
      probMSW[i,s] <- NrAll_msw[i,s]/NrAll_tot[i,s] 
      
      PropWgrilse[i,s]<-0
      PropWMSW[i,s]<-0
      
      
    }
    for(i in smolt_year[s]:(m+2)){
      SmoltWobs[i,s]~dlnorm(smoltPred[i,s],tau_SmoltW[i,s])
      R0[i,s]<-(K[s]*(EPR_M74[(i-e_delay[s]),s]-alphaSR[s]))/EPR_M74[(i-e_delay[s]),s] 
    }
  }
  
  cv_SR~dlnorm(-1.5,2)
  tau_SR<-1/log(cv_SR*cv_SR+1)
  mu_SR<- -0.5/tau_SR 
  
  #index is number of sea winters
  fec[1]~dlnorm(8.266503956,	54.78192169)
  fec[2]~dlnorm(9.139093315,	168.451323)
  fec[3]~dlnorm(9.504315997,	109.3565359)
  fec[4]~dlnorm(9.505447712,	61.70147645)
  fec[5]~dlnorm(9.68393554,	38.95507182)
  
  # Stocked as parr/smolts
  #mu_Parr, tau_Parr, Smolt_Rsp from ReleaseSimoTorne data file
  for(s in 1:stocks){ 
    for (i in 1:m){ 
      Parr[i,s] ~dlnorm(mu_Parr[i,s], tau_Parr[i,s])  
    }
    for (i in (m+1):(m+2)){ 
      Parr[i,s]<-0
    }
    
    for (i in 1:(m-1)){ 
      SmoltRsp[i,s] <-  Smolt_Rsp[i,s]      #was mu_SmoltT and mu_SmoltS
    }
    
    # Number of smolt releases in assessment year are assumed the same as year before
    SmoltRsp[m,s] <- SmoltRsp[(m-1),s]
  }
  
  # Smolt abundance	
  
  # Fit the latter part of the time series to the prior for the mean amount of observed
  # smolts. => smoltPred gets updated and contains afterwards information from both
  # model predicted amount of smolts (SR-part) and the mean amount of smolts based 
  # on the smolt mark-recapture model.
  # Data_SmoltsW.odc
  
  for (s in 1:stocks){
    
    # Early part of the time series, before the whole life history has gone through:
    for (i in 1:(smolt_year[s]-1)){
      
      R0[i,s]<-exp(M_R0[s])
      SmoltW[i,s] ~ T(dlnorm(mu_SmoltW[i,s], tau_SmoltW[i,s]),0.000001,3000000)  # mu_SmoltW and tau_SmoltW from Data_SmoltsW.odc
      
    }
    
    # From year 10-> we use predicted smolt abundance:
    for (i in smolt_year[s]:(m+2)){
      
      SmoltW[i,s] <- exp(smoltPred[i,s])
      
    }
    
    #     for (i in 1:(m)){
    # 	    IBSFC[i,s] <- SmoltW[i,s]/ R0[i,s]  
    # 	    prob.IBSFC50[i,s] <- step(IBSFC[i,s] - 0.5)	
    # 	    prob.IBSFC75[i,s] <- step(IBSFC[i,s] - 0.75)	
    # 	    prob.IBSFC100[i,s] <- step(IBSFC[i,s] - 1)										
    #   }
    
    ####################################Population dynamics########################################
    
    # Ncc: number of salmon on May 1st
    # Ndo: number of salmon on Jan 1st (imm)
    # Nl: number of salmon on Feb 1st (imm)
    # Nt: number of salmon on April 1st (imm)										 
    
    # Nc: number of salmon in coastal areas on June 1st (mat)
    # Nr: number of salmon in the river on Aug 1st (mat)
    # Nsp: number of spawners on Oct 1st (mat)
    
    
    #surv[,,1,] #natural mortality: 1 month for j= 1 (post-smolts), 0 months for j=2-6
    #surv[,,2,] #natural mortality + longline fishery: 2 months for j= 1 (post-smolts) and j=2-6 (immature)
    #surv[,,3,] #natural mortality at sea: 0 months for j=1, 8 months for j=2-6 (immature)
    #surv[,,4,] #natural mortality + offshore driftnet fishery: 1 month for j= 1 (post-smolts) and j=2-6 (immature)
    #surv[,,5,] #natural mortality + river fishery: 1 month for j= 1 (post-smolts), 2 months for j=2-6 (mature)
    #surv[,,6,] #natural mortality + offshore trolling fishery: 1 month for j= 1 (post-smolts) and j=2-6 (immature)
    
    #survc[,,]  #natural mortality + coastal trap net and gillnet fisheries: 6 months for j= 1 (post-smolts), 2 months for j=2-6 (mature)
    #survdc[,,] #natural mortality + coastal driftnet  fishery: 0 months for j=1, 1 month for j=2-6 (mature)
    
    
    for (i in 1:(m)){              
      
      #Post-smolts
      
      # Fish in the river in  the beginning of June in year 1 
      NrW[i,1,s] <- SmoltW[i,s]*survMps[i,1,s]   #SmoltW beginning of May (surv[i,1,1,1] 1 month)
      
      # Coastal fish in the beginning of July in year 1      
      NcW[i,1,s] <- NrW[i,1,s]*survMps[i,5,s]
      
      # Offshore salmon in the beginning of January in year 1 (calendar year 2)
      NdoW[i,1,s] <- NcW[i,1,s]*survMps[i,3,s]
      
      # Offshore salmon in the beginning of February in year 1  (calendar year 2)
      NlW[i,1,s] <- NdoW[i,1,s]*survMps[i,4,s]
      
      # Offshore salmon in the beginning of April in year 1  (calendar year 2)
      NtW[i,1,s] <- NlW[i,1,s]*survMps[i,2,s]		
      
      NccW[i,1,s]<-0
      NdcW[i,1,s]<-0  
      NdcWI[i,1,s]<-0
      NladderW[i,1,s]<-0  
      NspW[i,1,s]<-0
      
 # Offshore salmon in the beginning of May in year 2 
    NccW[i,2,s]<-NtW[i,1,s]*survMps[i,6,s] #*kE[i,1,6]	 
	 
    for (j in 3:6){  
    
        # Offshore salmon in the beginning of May in year 3-6  
		    NccW[i,j,s]<-NtW[i,(j-1),s]*surv[i,(j-1),6,1]*kE[i,(j-1),6]	 
	 } 
   for (j in 2:6){ 
        
        #Immature salmon in coastal areas in the beginning of May    
        NdcWI[i,j,s] <- NccW[i,j,s]* (1-LW[(i+j-1),(j-1)])  #index 1 of LW is current year, 2 of LW is sea winters  
        
        # Offshore salmon in the beginning of January 
        NdoW[i,j,s] <- NdcWI[i,j,s]*surv[i,j,3,1]
        
        #Offshore salmon in the beginning of February
        NlW[i,j,s] <- NdoW[i,j,s]*surv[i,j,4,1]
        
        # Offshore salmon in the beginning of April        
        NtW[i,j,s] <- NlW[i,j,s]*surv[i,j,2,1]										   
        
        #Salmon in coastal areas in the beginning of May.     #Mature
        NdcW[i,j,s] <- NccW[i,j,s]*LW[(i+j-1),(j-1)]   #Maturation
        
        # Salmon in the coastal area in the beginning of June 
        NcW[i,j,s] <- NdcW[i,j,s]*survdc[i,j,1,AU[s]]
        
        # Salmon in the river in the beginning of August                  #Mature
        NrW[i,j,s] <- NcW[i,j,s]*survc[i,j,1,AU[s]]
        
        # Number of salmon available to river fishery and spawning        #Mature
        #added to deal with different situation in Ume, p.ladder=1 for all rivers except Ume
        NladderW[i,j,s]<-NrW[i,j,s]*p.ladder[(i+j-1),s]	   
        
        # Number of spawners in the river in the beginning of October    #Mature

        NspW[i,j,s]<-NladderW[i,j,s]*surv.riv[i,j,s]*surv_migr[(i+j-1),s] 	 
        
      }
    }  
    
  } #s
  
  for (i in 1:(m)){             
    for(s in 1:2){    #Torne, Simo
      
      #Post-smolts
      
      # Fish in the river in  the beginning of June in year 1 
      NrRsp[i,1,s] <- (SmoltRsp[i,s] + Parr[i,s])*surv[i,1,1,2]
      
      # Coastal fish in the beginning of July in year 1      
      NcRsp[i,1,s] <- NrRsp[i,1,s]*(1-HrW[i,1])*exp(-MpsR[i]*sealMort[i,1,1]/Tstep)
       
       
      # Offshore salmon in the beginning of January in year 1 (calendar year 2)
      NdoRsp[i,1,s] <- NcRsp[i,1,s]*survc[i,1,2,AU[s]]   
      
      # Offshore salmon in the beginning of February in year 1  (calendar year 2)
      NlRsp[i,1,s] <- NdoRsp[i,1,s]*surv[i,1,4,2]  
      
      
      # Offshore salmon in the beginning of April in year 1  (calendar year 2)
      NtRsp[i,1,s] <- NlRsp[i,1,s]*surv[i,1,2,2]	  
      
      NccRsp[i,1,s]<-0
      NdcRsp[i,1,s]<-0  
      NdcRspI[i,1,s]<-0  
      NladderRsp[i,1,s]<-0
      NspRsp[i,1,s]<-0
      
      for (j in 2:6){ 
        
        # Offshore salmon in the beginning of May in year 2-6  

        NccRsp[i,j,s] <- NtRsp[i,(j-1),s]*surv[i,(j-1),6,2]
        
        
        #Immature salmon in coastal areas in the beginning of May    
        NdcRspI[i,j,s] <-NccRsp[i,j,s]*(1-LR[(i+j-1),(j-1)])
        
        # Offshore salmon in the beginning of January 
        NdoRsp[i,j,s]<- NdcRspI[i,j,s]*surv[i,j,3,2]  
        
        #Offshore salmon in the beginning of February
        NlRsp[i,j,s] <- NdoRsp[i,j,s]*surv[i,j,4,2]
        
        #Offshore salmon in the beginning of April
        NtRsp[i,j,s] <- NlRsp[i,j,s]*surv[i,j,2,2]														 
        
        #Salmon in coastal areas in the beginning of May.     #Mature
        NdcRsp[i,j,s] <- NccRsp[i,j,s]*LR[(i+j-1),(j-1)]	            
        
        # Salmon in the coastal area in the beginning of June 
        NcRsp[i,j,s] <- NdcRsp[i,j,s]*survdc[i,j,2,AU[s]]        
        
        # Salmon in the river in the beginning of August              #Mature
        NrRsp[i,j,s] <- NcRsp[i,j,s]*survc[i,j,2,AU[s]]
        
        # Number of salmon available to river fishery and spawning         #Mature
        #added to deal with different situation in Ume. p.ladder=1 for all rivers except Ume
        NladderRsp[i,j,s] <- NrRsp[i,j,s]*p.ladder[(i+j-1),s]		 	    
        
        # Number of spawners in the river in the beginning of October    #Mature
        NspRsp[i,j,s] <- NladderRsp[i,j,s]*(1-HrW[i,j])* exp(-2*MR/Tstep)*surv_migr[(i+j-1),s]

      }
    }
  }
  
  
  
  for(s in 1:AUS){              #s index denotes AU here 
    for (i in 1:m){ 
      
      #Post-smolts
      
      # Fish in the river in  the beginning of June in year 1 
      
      SmoltR[i,s] <- SmoltRdata[i,s] * Usmolt               #SmoltRdata from SmoltR data file Usmolt from
      NrR[i,1,s] <- SmoltR[i,s]*surv[i,1,1,2]                                                    
      
      # Coastal fish in the beginning of July in year 1
      NcR[i,1,s] <- NrR[i,1,s]*surv[i,1,5,2]

      # Offshore salmon in the beginning of January in year 1
      NdoR[i,1,s] <- NcR[i,1,s]*survc[i,1,2,s]
      
      # Offshore salmon in the beginning of February in year 1
      NlR[i,1,s] <- NdoR[i,1,s]*surv[i,1,4,2]
      
      # Offshore salmon in the beginning of April in year 1
      NtR[i,1,s] <- NlR[i,1,s]*surv[i,1,2,2]											   
      
      NccR[i,1,s]<-0
      NdcR[i,1,s]<-0
      NdcRI[i,1,s]<-0
      NspR[i,1,s]<-0

     for (j in 2:6){ 
        
        # Offshore salmon in the beginning of May in year 2
     
        NccR[i,j,s] <- NtR[i,(j-1),s]*surv[i,(j-1),6,2]

        ##Immature salmon in coastal areas in the beginning of May       
        NdcRI[i,j,s] <- NccR[i,j,s]*(1-LR[(i+j-1),(j-1)])

        # Offshore salmon in the beginning of January 
        NdoR[i,j,s] <-  NdcRI[i,j,s]*surv[i,j,3,2]
        
        # Offshore salmon in the beginning of February
        NlR[i,j,s] <- NdoR[i,j,s]*surv[i,j,4,2]

        # Offshore salmon in the beginning of April
        NtR[i,j,s] <- NlR[i,j,s]*surv[i,j,2,2]															 
        
        #Salmon in coastal areas in the beginning of May.     #Mature
        NdcR[i,j,s] <- NccR[i,j,s]*LR[(i+j-1),(j-1)]
 
        # Salmon in the coastal area in the beginning of June 
        NcR[i,j,s] <- NdcR[i,j,s]*survdc[i,j,2,s]
    
        # Salmon in the river in the beginning of August              #Mature
        NrR[i,j,s]<- NcR[i,j,s]*survc[i,j,2,s]
        
        # Number of spawners in the river in the beginning of October    #Mature
        NspR[i,j,s]  <- NrR[i,j,s]*surv[i,j,5,2]	
 
      }			
    }
  }
  
  for (i in 1:m){                    ### Tagging data likelihoods ### 
      
    for(j in 1:5){     #check - seems ok because all immature fish
      
      # Offshore driftnet fishery
      HdoW[i,j] <- 1-exp(-qdW[i,j]*Edo[i,j]) 
      HdoR[i,j] <- 1-exp(-qdR[i,j]*Edo[i,j]) 

      # Offshore longline fishery 
      HlW[i,j]<-1-exp(-qlW[i,j]*El[i,j]) 
      HlR[i,j]<-1-exp(-qlR[i,j]*El[i,j]) 

      
      ################## Total catch calculations ######################
      
      NdoW_all[i,j]<-sum(NdoW[i,j,1:stocks])
      NlW_all[i,j]<- sum(NlW[i,j,1:stocks])
      NtW_all[i,j]<- sum(NtW[i,j,1:stocks])									 
      
      NdoR_all[i,j]<- sum(NdoR[i,j,1:AUS])+sum(NdoRsp[i,j,1:2])
      NlR_all[i,j]<- sum(NlR[i,j,1:AUS])+sum(NlRsp[i,j,1:2])
      NtR_all[i,j]<- sum(NtR[i,j,1:AUS])+sum(NtRsp[i,j,1:2])
      
      # Estimated catches of non-returning salmon in the offshore fishery
      nco_W[i,j] <- PropCW[i] * (HdoW[i,j]*NdoW_all[i,j]+HlW[i,j]*NlW_all[i,j])   
      nco_R[i,j] <- PropCR[i] * (HdoR[i,j]*NdoR_all[i,j]+HlR[i,j]*NlR_all[i,j])   
      nco[i,j] <- nco_W[i,j] + 	nco_R[i,j] 	
      
      # Estimated catches separately for trolling, longline and driftnet 
      nc_otr[i,j] <- PropCW[i] * (HtW[i,j]*NtW_all[i,j]) + PropCR[i] *(HtR[i,j]*NtR_all[i,j])
      nc_otrW[i,j] <- PropCW[i] * (HtW[i,j]*NtW_all[i,j]) #wild only 2024																					  
      nc_oll[i,j] <- PropCW[i] * (HlW[i,j]*NlW_all[i,j]) + PropCR[i] *(HlR[i,j]*NlR_all[i,j])	    #PropCW and PropCR from PropAU16 data file
      nc_odn[i,j] <- PropCW[i] * (HdoW[i,j]*NdoW_all[i,j]) + PropCR[i] * (HdoR[i,j]*NdoR_all[i,j]) 
      
    }    
    
    HdcW[i,1]<-0
    HdcR[i,1]<-0
    
    for(j in 2:6){ 
      HdcW[i,j] <- 1-exp(-qdW[i,j-1]*Edc[i,j])   #WinBUGS j 2:6  HdcW[i,j-1] <- 1-exp(-qdW[j-1]*Edc[i,j]) 
      HdcR[i,j] <- 1-exp(-qdR[i,j-1]*Edc[i,j]) 
  
    }
    
    for (j in 1:6){
      for(au in 1:AUS){   
  
        HcW[i,j,au] <- (1-exp(-(F_sea[i,j,1,au])))	
        HcR[i,j,au] <- (1-exp(-(F_sea[i,j,2,au])))
        
      }
   
      ################## Total catch calculations ######################

            # Morrum and Eman not included why not included in river?
      #NrW_all[i,j]<- sum(NladderW[i,j,1:13])+sum(NladderW[i,j,16:17])+sum(NladderRsp[i,j,1:2])
      NrW_all[i,j]<- sum(NladderW[i,j,1:17]*(HrW[i,j]*rivHR[(i+j-1),1:17]))+sum(NladderRsp[i,j,1:2]*(HrW[i,j]*rivHR[(i+j-1),1:2]))	  #1:13, 16  #08.06.2017 Ume included

      NdcW_all[i,j]<-sum(NdcW[i,j,1:13])+sum(NdcW[i,j,16:17])	
      
      NrR_all[i,j]<- sum(NrR[i,j,1:AUS])
      NdcR_all[i,j]<- sum(NdcR[i,j,1:(AUS-1)])+sum(NdcRsp[i,j,1:2]) #not AU4     
      
      NcR_all[i,j,1]<- NcR[i,j,1]+sum(NcRsp[i,j,1:2])		
      NcR_all[i,j,2]<- NcR[i,j,2]  					
      NcR_all[i,j,3]<- NcR[i,j,3]  		
      
      # Estimated catches of wild and hatchery-reared salmon in the river
      ncr[i,j]<- NrW_all[i,j]+HrR[i,j]*NrR_all[i,j]	
      
      # Estimated catches of hatchery-reared salmon in the river
      #ncrR[i,j]<- HrR[i,j]*NrR_all[i,j]  
      
      # Estimated catches of salmon in the coastal areas
      for(s in 1:stocks){
        nccs[i,j,s]<-HcW[i,j,AU[s]]*NcW[i,j,s] 
      }
#      for(s in 1:2){
#         ncrW[i,j,s]<-HrW[i,j]*NladderW[i,j,s]+NladderRsp[i,j,s]*HrW[i,j]    #rivHR needed?
#      }
#      for(s in 3:stocks){
#         ncrW[i,j,s]<-HrW[i,j]*NladderW[i,j,s]
#      }
      
      ncc[i,j] <- sum(nccs[i,j,1:stocks]) + inprod(HcR[i,j,1:(AUS-1)],NcR_all[i,j,1:(AUS-1)])+ HdcW[i,j]*NdcW_all[i,j]+HdcR[i,j]*NdcR_all[i,j]  #not AU4	   
      
    }	#j       
  }  #i
  
  
  for(i in 1:5){             #move above later
    
    ncrR_Tot[i]<-0
    ncc_Tot[i]<-0
    nco_Tot[i]<-0
    nct_Tot[i]<-0
    nctW_Tot[i]<-0
    nctW_rel[i]<-0
    
    nct_ObsTotX[i]<-0				 
    nco_ObsTotX[i]<-0
    ncc_ObsTotX[i]<-0
    ncr_ObsTotX[i]<-0

     for(s in 1:stocks){
       ncrW_tot[i,s]<-0
       ncrW_ObsTotX[i,s]<-0
  }  
  }
  for (i in 6:(m-1)){
  
  ##NEW stock-specific likelihood river catches
#  for(s in 1:stocks){
#       ncrW_Tot[i,s]<-ncrW[(i-1),2,s] + ncrW[(i-2),3,s] + ncrW[(i-3),4,s] + ncrW[(i-4),5,s] + ncrW[(i-5),6,s]
#       muCRW[i,s] <- log(ncrW_Tot[i,s] /ureport_r[i]) - 0.5/tauCRW
#    ncrW_ObsTot[i,s] ~ dlnorm(muCRW[i,s], tauCRW)	 #change data name and make it reared rivers only
#    ncrW_ObsTotX[i,s] ~ dlnorm(muCRW[i,s], tauCRW)	  
#  }  
 
    
    nct_Tot[i] <- nc_otr[(i-1),2]+nc_otr[(i-2),3]+nc_otr[(i-3),4]+nc_otr[(i-4),5]   #Estimated total trolling catches
    nctW_Tot[i] <- nc_otrW[(i-1),2]+nc_otrW[(i-2),3]+nc_otrW[(i-3),4]+nc_otrW[(i-4),5] 
    nctW_rel[i] <-  nctW_Tot[i]*p.rel[i]  
    muCT[i] <-log(nct_Tot[i] /ureport_o[i]) - 0.5/tauCT	
    nct_ObsTot[i]~dlnorm(muCT[i], tauCT)
    nct_ObsTotX[i]~dlnorm(muCT[i], tauCT)																											 						 
    
    nco_Tot[i] <- nco[(i-1),2]+nco[(i-2),3]+nco[(i-3),4]+nco[(i-4),5]                 #Estimated total offshore catches Nl and Ndo
    muCO[i] <-log(nco_Tot[i] /ureport_o[i]) - 0.5/tauCO	
    nco_ObsTot[i] ~dlnorm(muCO[i], tauCO)
    nco_ObsTotX[i] ~dlnorm(muCO[i], tauCO)
    
    nc_oll_Tot[i] <- nc_oll[(i-1),2]+nc_oll[(i-2),3]+nc_oll[(i-3),4]+nc_oll[(i-4),5]         #Estimated total offshore longline catches Nl only 
    nc_odn_Tot[i] <- nc_odn[(i-1),2]+nc_odn[(i-2),3]+nc_odn[(i-3),4]+nc_odn[(i-4),5]         #Estimated total offshore driftnet catches Ndo only
    
    ncr_Tot[i] <- ncr[(i-1),2]+ncr[(i-2),3]+ncr[(i-3),4]+ncr[(i-4),5]+ ncr[(i-5),6]     # Estimated total river catches   Nr
    muCR[i] <- log(ncr_Tot[i] /ureport_r[i]) - 0.5/tauCR
    ncr_ObsTot[i] ~dlnorm(muCR[i], tauCR)	
    ncr_ObsTotX[i] ~dlnorm(muCR[i], tauCR)
    
#    ncrR_Tot[i] <- ncrR[(i-1),2]+ncrR[(i-2),3]+ncrR[(i-3),4]+ncrR[(i-4),5]+ ncrR[(i-5),6]     # Estimated total river catches   Nr
#    #ncr_Tot[i,s]<-  ncrR_Tot[i] + sum(ncrW_Tot[i,1:stocks])
#    #muCR[i] <- log(ncr_Tot[i] /ureport_r[i]) - 0.5/tauCR
#    
#    muCR[i] <- log(ncrR_Tot[i] /ureport_r[i]) - 0.5/tauCR
#    ncrR_ObsTot[i] ~dlnorm(muCR[i], tauCR)	 #change data name and make it reared rivers only
#    ncrR_ObsTotX[i] ~dlnorm(muCR[i], tauCR)
    
    ncc_Tot[i] <- ncc[(i-1),2]+ncc[(i-2),3]+ncc[(i-3),4]+ncc[(i-4),5]+ ncc[(i-5),6]     # Estimated total coastal catches Nc Ndc
    muCC[i] <- log(ncc_Tot[i] /ureport_c[i]) - 0.5/tauCC	
    ncc_ObsTot[i] ~dlnorm(muCC[i], tauCC)
    ncc_ObsTotX[i] ~dlnorm(muCC[i], tauCC)
    
    
    ######### Spawner counting observation models (Didson and fishladders) #########
    
    mu.sp[i,1]<-get_mu_tau(NrWtot[i,1],p.detect[i,1])[1]
    tau.sp[i,1]<-get_mu_tau(NrWtot[i,1],p.detect[i,1])[2]
    
    mu.l[i,1]<-get_mu_tau(NrWtot[i,1],p.ladder[i,1])[1]
    tau.l[i,1]<-get_mu_tau(NrWtot[i,1],p.ladder[i,1])[2]
    
    sp_count[i,1]~dnorm(mu.sp[i,1],tau.sp[i,1]) # Torne Didson count 
    sp_countX[i,1]~dnorm(mu.sp[i,1],tau.sp[i,1])
    
    ladder_count[i,1]~dnorm(mu.l[i,1],tau.l[i,1])   #Ume only (propn that finds the ladder)
    
    
    #sp_count[i,1]~dbin(p.detect[i,1],NrWtot[i,1]) # Torne Didson count 
    #sp_countX[i,1]~dbin(p.detect[i,1],NrWtot[i,1]) 
    
    #ladder_count[i,1]~dbin(p.ladder[i,1],NrWtot[i,1])   #Ume only (propn that finds the ladder)
    
    eta_star[i,1]<-N_sp_count[i,1]*(eta_msw[1]+1)/(eta_msw[1]+N_sp_count[i,1])-1    #N_sp_count=sp_count (observation)
    alpha_msw[i,1]<-probMSW[i,1]*eta_star[i,1]+1
    beta_msw[i,1]<-(1-probMSW[i,1])*eta_star[i,1]+1
    MSWprop[i,1]~dbeta(alpha_msw[i,1],beta_msw[i,1])
    
    sp_count[i,2]~dlnorm(muDS[i], tauDS) # Simojoki Didson count
    sp_countX[i,2]~dlnorm(muDS[i], tauDS) 
    
    #muDS[i]<-log(NrWtot[i,2]/coefDS)-0.5*(1/tauDS)
   	NrW_msw_Simo[i]<-round((NrW_msw[i,2]+NrRsp_msw[i,2])*1000)+1   
	  muDS[i]<-log(NrW_msw_Simo[i]/coefDS)-0.5*(1/tauDS) # Without hierarchy
    
    for(s in 3:stocks){
      
      mu.sp[i,s]<-get_mu_tau(NrWtot[i,s],p.detect[i,s])[1]
      tau.sp[i,s]<-get_mu_tau(NrWtot[i,s],p.detect[i,s])[2]
      
      mu.l[i,s]<-get_mu_tau(NrWtot[i,s],p.ladder[i,s])[1]
      tau.l[i,s]<-get_mu_tau(NrWtot[i,s],p.ladder[i,s])[2] 
      
      sp_count[i,s]~dnorm(mu.sp[i,s],tau.sp[i,s]) # Torne Didson count 
      sp_countX[i,s]~dnorm(mu.sp[i,s],tau.sp[i,s]) 
      
      ladder_count[i,s]~dnorm(mu.l[i,s],tau.l[i,s])   #Ume only (propn that finds the ladder) 
      
      
      
      eta_star[i,s]<-N_sp_count[i,s]*(eta_msw[s]+1)/(eta_msw[s]+N_sp_count[i,s]) #-1    #N_sp_count=sp_count (observation)
      alpha_msw[i,s]<-probMSW[i,s]*eta_star[i,s]+1
      beta_msw[i,s]<-(1-probMSW[i,s])*eta_star[i,s]+1
      MSWprop[i,s]~dbeta(alpha_msw[i,s],beta_msw[i,s])
      
    }
    
    WGrilse[i,1]~dbin(PropWgrilse[i,1], Grilse_all[i,1]) #proportion of wild grilse
    WMSW[i,1]~dbin(PropWMSW[i,1], MSW_all[i,1])
    
    Wprop[i,1]<-nco_W[i,2]/(nco_W[i,2]+nco_R[i,2])  
    Wprop[i,2]<-nco_W[i,3]/(nco_W[i,3]+nco_R[i,3])

    for (j in 1:2){     #1SW (j=2 elsewhere) and 2SW (j=3 elsewhere) salmon
       eta_wr[i,j]<-(Wprop[i,j]*(1-Wprop[i,j]))/(sd_wr[i,j]*sd_wr[i,j])#-1  
  
       log_WpropObs[i,j]~dnorm(log(Wprop[i,j])-0.5*log(1/(Wprop[i,j]*eta_wr[i,j])+1),
              1/log(1/(Wprop[i,j]*eta_wr[i,j])+1))
       log_RpropObs[i,j]~dnorm(log(1-Wprop[i,j])-0.5*log(1/((1-Wprop[i,j])*eta_wr[i,j])+1),
              1/log(1/((1-Wprop[i,j])*eta_wr[i,j])+1))
    }  
    
    for(rs in 1:rstocks){  #1=Lule?lven, 2=Dal?lven
      
      NrRtot.vul[i,rs]<-round((NrR[(i-1),2,AUR[rs]]*RProp[(i-1),rs]+NrR[(i-2),3,AUR[rs]]*RProp[(i-2),rs]+
                                 NrR[(i-3),4,AUR[rs]]*RProp[(i-3),rs]+NrR[(i-4),5,AUR[rs]]*RProp[(i-4),rs]+
                                 NrR[(i-5),6,AUR[rs]]*RProp[(i-5),rs])*1000)
      
      HRNrR[i,rs]<-min(CatchR[i,rs]/NrRtot.vul[i,rs],0.9999)
      NrRtot[i,rs]<-NrRtot.vul[i,rs]*(1-HRNrR[i,rs])
      
   #   cv.trap[i,rs]<-(sqrt(NrRtot[i,rs]*pTrap[i,rs]*(1-pTrap[i,rs])))/(NrRtot[i,rs]*pTrap[i,rs])
#      s.trap[i,rs]<-sqrt(log(cv.trap[i,rs]*cv.trap[i,rs]+1))
#      mu.trap[i,rs]<-log((NrRtot[i,rs]*pTrap[i,rs])/exp(0.5*s.trap[i,rs]^2))
#      tau.trap[i,rs]<-1/(s.trap[i,rs]*s.trap[i,rs])
      
      pTrap[i,rs]~dbeta(aTrap[rs],bTrap[rs])
      
      mu.trap[i,rs]<-get_mu_tau(NrRtot[i,rs],pTrap[i,rs])[1]
      tau.trap[i,rs]<-get_mu_tau(NrRtot[i,rs],pTrap[i,rs])[2]

      
      #TrapTot[i,rs]~dbin(pTrap[i,rs],NrRtot[i,rs])
      TrapTot[i,rs]~dnorm(mu.trap[i,rs],tau.trap[i,rs])
      
    }
  }
  
  for(i in 1:5){
    NLuleRec[i]~dbin(pTrap[yLule[i],1],NLuleRel[i])
  }	
  
  
  ###################MORTALITY AND SURVIVAL RATES#####################
  
  ## Instantaneous adult natural mortality rate
  
  MW ~T(dlnorm(-2.3, 4.3),0.025,0.35)
  MR ~T(dlnorm(-2.3, 4.3),0.025,0.35)
#  
#  early_MpsW~T(dlnorm(0.23,19),0.5,5)

#  MW ~dlnorm(-2.3, 4.3)
#  MR ~dlnorm(-2.3, 4.3)
  
 
  
   for(st in 1:stocks){
  
        logit_survMps[1,st] ~ T(dnorm(0,0.1),-3,-0.4)  #max ~40%
        MpsW[1,st]<- -log((exp(logit_survMps[1,st])/(1+exp(logit_survMps[1,st]))))
    } #st
    for(t in 1:(m-1)){
        mu_surv[t,1:stocks]<-logit_survMps[t,1:stocks]
        logit_survMps[(t+1),1:stocks]~dmnorm(mu_surv[t,1:stocks],tau_Mps[1:stocks,1:stocks]) 
        for(st in 1:stocks){
            MpsW[(t+1),st]<- -log((exp(logit_survMps[(t+1),st])/(1+exp(logit_survMps[(t+1),st]))))
       
        }
    }  #t
    
    ##############################################
    ## Model correlation matrix based on cholesky onion method ##
    ##############################################
  
        alpha[1]   <- eta + (stocks - 2)/2  
        corY[1]~dbeta(alpha[1], alpha[1])
        r12   <- 2 * corY[1] - 1
        ##
        R[1,1]     <- 1                 
        R[1,2]     <- r12               
        R[2,2]     <- sqrt(1 - r12^2)   
        R[2:stocks,1]   <- 0    #with > 2 stocks (plus lines below)
          
          for (i in 2:(stocks-1)) {         #NB stocks needs to be >2, otherwise just use lines above!!
            ## Draw beta random variable
            alpha[i] <- alpha[(i-1)] - 0.5 
            corY[i] ~ dbeta(i / 2, alpha[i]) 
            ## Draw uniformly on a hypersphere
            for (jj in 1:i) {
              corZ[i, jj] ~ dnorm(0, 1)
            }
            scZ[i, 1:i]  <- corZ[i, 1:i] / sqrt(inprod(corZ[i, 1:i], corZ[i, 1:i]))
            R[1:i,(i+1)]   <- sqrt(corY[i]) * scZ[i,1:i]
            R[(i+1),(i+1)]   <- sqrt(1 - corY[i])
            R[(i+1):stocks,i] <- 0
          }  #i
        
        Rnew[1:stocks,1:stocks]<-t(R[1:stocks,1:stocks]) %*% R[1:stocks,1:stocks] 
        ##
        for(st in 1:stocks){
            sigma_stock[st]  ~ T(dlnorm(-1.61,4),0.001,0.25)    #change upper bound later?
            for(j in 1:stocks){
                sigma_st[st,j] <- Rnew[st,j] * sigma_stock[st] * sigma_stock[j]
            }
        }
 
    tau_Mps[1:stocks, 1:stocks]<- inverse(sigma_st[1:stocks, 1:stocks])
  
  Reff_mu~T(dbeta(1,2),0.01,0.99)
  Reff_eta~dunif(0.01,0.5) 
  Ra<-Reff_mu/Reff_eta
  Rb<-(1-Reff_mu)/Reff_eta   

  for(i in 1:(m)){
    for(s in 1:stocks){
            survMps[i,1,s]<-exp(-MpsW[i,s]/Tstep) 		          #survmortW
            survMps[i,2,s] <- exp(-qlW[i,1]*El[i,1]) * exp(-2*MpsW[i,s]/Tstep) 	       #survlW  
            survMps[i,3,s] <- exp(-F_sea[i,1,1,AU[s]])* exp(-(6*MpsW[i,s]/Tstep))   	      #was survc
            survMps[i,4,s] <- exp(-qdW[i,1]*Edo[i,1]) * exp(-MpsW[i,s]/Tstep)         #survdoW
            survMps[i,5,s] <- (1-HrW[i,1])*exp(-MpsW[i,s]*sealMort[i,1,1]/Tstep)  #was surv.riv
            survMps[i,6,s]<-(1-HtW[i,1]*(p.rel[(i+1-1)]*p.mort+(1-p.rel[(i+1-1)])))* exp(-MpsW[i,s]/Tstep) 
        }
    
    RMps[i]~dbeta(Ra,Rb) 
    ReffectMps[i] <- (RMps[i] * 1.5)+1                 # hatchery-reared effect between 1 and 2.5
    MpsR[i] <- MpsW[i,1] * ReffectMps[i]
    
    surv[i,1,1,2]<-exp(-MpsR[i]/Tstep) 	            #survmortR
    surv[i,1,2,2] <- exp(-qlR[i,1]*El[i,1]) * exp(-2*MpsR[i]/Tstep)          #survlR
    surv[i,1,3,2] <- 0.99
    surv[i,1,4,2] <- exp(-qdR[i,1]*Edo[i,1]) * exp(-MpsR[i]/Tstep) 	 
    surv[i,1,5,2]<-(1-HrR[i,1])*exp(-MpsR[i]*sealMort[i,1,1]/Tstep)    #survrR
    surv[i,1,6,2]<-(1-HtR[i,1]) * exp(-MpsR[i]/Tstep)	
    
    #pmort is the sum of p.release * M.inc + (1-p.release)
    #p.release in 2022 is 100%, M.inc is 25%																											   
    
    for (au in 1:3){	
      
      F_sea[i,1,1,au]<-qctnW[1,au]*Ectn[i,1,au]+qcgnW[1,au]*Ecgn[i,1,au]														
      
      F_sea[i,1,2,au]<-qctnR[1,au]*Ectn[i,1,au]+qcgnR[1,au]*Ecgn[i,1,au]		
      survc[i,1,2,au] <-exp(-F_sea[i,1,2,au]) * exp(-(6*MpsR[i]/Tstep))   
      
      survdc[i,1,1,au] <- 0.99    #not used anywhere
      survdc[i,1,2,au] <- 0.99
    }
    
    F_sea[i,1,1,4]<-0
    F_sea[i,1,2,4]<-0	
      
    survc[i,1,2,4] <-exp(-(6*MpsR[i]/Tstep))      #reared 
    
    survdc[i,1,1,4] <- 0.99   
    survdc[i,1,2,4] <- 0.99
    
    for(j in 2:6){
      
      surv[i,j,1,1]<-0.99 		          #survmortW  not used anywhere
      surv[i,j,1,2]<-0.99 	            #survmortR
      
      #surv[i,j,2,1] <- exp(-qlW[j]*El[i,j]) * exp(-3*MW/Tstep)                        #survlW  
      #surv[i,j,2,2] <- exp(-qlR[j]*El[i,j]) * exp(-3*MR/Tstep)		 
      surv[i,j,2,1] <- exp(-qlW[i,j]*El[i,j]) * exp(-2*MW/Tstep)                        #survlW  
      surv[i,j,2,2] <- exp(-qlR[i,j]*El[i,j]) * exp(-2*MR/Tstep)
      
      surv[i,j,3,1] <- exp(-8*MW/Tstep)	#QUERY 2 time steps in process error calc!	      #survnohW
      surv[i,j,3,2] <- exp(-8*MR/Tstep)
      
      surv[i,j,4,1] <- exp(-qdW[i,j]*Edo[i,j]) * exp(-MW/Tstep)       #survdoW
      surv[i,j,4,2] <- exp(-qdR[i,j]*Edo[i,j]) * exp(-MR/Tstep) 
      
      surv[i,j,5,1] <-1 #(1-HrW[i,j])*exp(-2*MW/Tstep)     #survrW                
      surv[i,j,5,2] <- (1-HrR[i,j])* exp(-2*MR/Tstep)  	
      
      for(s in 1:stocks){
           surv.riv[i,j,s] <- (1-HrW[i,j])*exp(-2*MW/Tstep)     #survrW
        }
         

      
      surv[i,j,6,1]<-(1-HtW[i,j]*(p.rel[(i+j-1)]*p.mort+(1-p.rel[(i+j-1)]))) * exp(-MW/Tstep) 	         
      surv[i,j,6,2]<-(1-HtR[i,j]) * exp(-MR/Tstep)
                        
      for(au in 1:3){
        
        # Survival from natural and fishing mortality during June and July
        
        F_sea[i,j,1,au]<-qctnW[j,au]*Ectn[i,j,au]+qcgnW[j,au]*Ecgn[i,j,au]		
        survc[i,j,1,au] <-exp(-F_sea[i,j,1,au]) * exp(-(2*MW*sealMort[i,j,au]/Tstep))
        
        # Survival from natural and fishing mortality during June and July	
        F_sea[i,j,2,au]<-qctnR[j,au]*Ectn[i,j,au]+qcgnR[j,au]*Ecgn[i,j,au]				
        survc[i,j,2,au] <-exp(-F_sea[i,j,2,au])  * exp(-(2*MR*sealMort[i,j,au]/Tstep)) 
        
        survdc[i,j,1,au] <- exp(-qdW[i,j-1]*Edc[i,j]) * exp(-MW/Tstep)      
        survdc[i,j,2,au] <- exp(-qdR[i,j-1]*Edc[i,j]) * exp(-MR/Tstep)
      }
      
      F_sea[i,j,1,4]<-0
      F_sea[i,j,2,4]<-0
      
      survc[i,j,1,4] <- exp(-(2*MW*sealMort[i,j,4]/Tstep))   #wild   #seal M AU4 = 1   
      survc[i,j,2,4] <- exp(-(2*MR*sealMort[i,j,4]/Tstep))   #reared 
      
      survdc[i,j,1,4] <- exp(-MW/Tstep)    
      survdc[i,j,2,4] <- exp(-MR/Tstep)  	      					
    }
  }
  

  ###################################### Catchabilities & harvest rates ####################################
  
  # Mean reverting AR(1)-model for 1SW (==MSW) harvest rate in offshore trolling fishery
  # =======================================================================
  # Note that autocorrelation takes place between years but in HtW[i,j] index i is smolt cohort 
  
  # Recreational offshore trolling 
  for(i in 1:(m+4)){ # i: smolt cohort
    # Post-smolts are assumed released without mortality
    HtW[i,1]<-0 
    HtR[i,1]<-0
    
    # i+1=i+2-1: cohort i age 2 transformed to year i+1
    logit(HtW[i,2])<-logitHtW2[i+1] 
    logit(HtR[i,2])<-logitHtW2[i+1]#<-logitHtR2[i+1] 
  }
  
  for(i in 1:m){ # i: smolt cohort 
    for(j in 3:6){ 
      # Harvest rate of MSW salmon of different age must be the same within years:
      #
      # The cohort i salmon are of age 3 in year y=i+3-1=i+2 and 
      # in the same year, the age 2 salmon originate from cohort y-1=i+2-1= i+1 (=i+3-2)
      #
      # The cohort i salmon are of age 4 in year y=i+4-1=i+3 and
      # in the same year the age 2 salmon originate from cohort y-1=i+3-1=i+2 (=i+4-2)
      # and so on...
      HtW[i,j]<-HtW[i+j-2,2] 
      HtR[i,j]<-HtW[i,j]#HtR[i+j-2,2]
    }
  }
  
  for(i in 1:(m+4)){ # i: calendar year
    logitHtW2[i+1]~dnorm(mu_trW[i+1],tau_tr)
    mu_trW[i+1]<-phi_tr*logitHtW2[i]+(1-phi_tr)*mean_trW # mean reverting AR(1)
    
    #logitHtR2[i+1]~dnorm(mu_trR[i+1],tau_tr)
    # mu_trR[i+1]<-phi_tr*logitHtR2[i]+(1-phi_tr)*mean_trR
    
  }
  
  tau_tr<-1/((1-pow(phi_tr,2))*(sd_tr*sd_tr))  
  # Marginal variance chosen to give uniform when mean=0, otherwise unimodal
  sd_tr~dunif(0.01,1.6)
  #phi_tr~dunif(0.01,0.99) #positive autocorrelation
  phi_tr~dunif(0.50,0.99) #positive autocorrelation
  mean_trW~dnorm(0,0.39) # implies uniform[0,1] prior for mean harvest rate
  #mean_trR~dnorm(0,0.39) 
  logitHtW2[1]~dnorm(0,0.39) 

  
   
  for(i in 1:m){		  
    # Harvest rate of smolts in the river/river mouth		
    
    HrW[i,1] <- 1-exp(-(qrW[1]*Er[i,1]))
    HrR[i,1] <- 1-exp(-(qrR[1]*Er[i,1])) 	
    
    for (j in 2:3){ 
      HrR[i,j]<-1-exp(-(qrR[j]*Er[i,j]))	    # Harvest rate river fishery reared salmon
    }
    for (j in 4:6){ 
      HrW[i,j] <- HrW[i+j-3,3] #!!!! equal by year not by smolt cohort
      HrR[i,j]<-1-exp(-(qrR[j]*Er[i,j]))	 
    }
  }
  for(i in 1:(m+3)){ # HrW age 3 must go +3 further to match dimensions above		  
    for (j in 2:3){ 
      HrW[i,j]~dbeta(1.6,6.4)	# Harvest rate river fishery in rivers with natural reproduction
    }}



  #post-smolts, j=1
  # Harvest rates in 1987
  # ==================
  HRR[1]  ~dbeta (1,20)	# Harvest rate smolts in river fishery
  qrW[1] ~ dlnorm(mqr[1], tauqr)	# Catchability of wild smolts in the river			
  qrR[1] ~ dlnorm(mqr[1], tauqr)	# Catchability of reared smolts in the river			
  
  mqr[1] <- log(-log(1-HRR[1])/Er[1,1])	# Mean catchability river fishery		
  
  HRL[1] ~ dbeta(1,20)   # Harvest rate of smolts longline fishery	
  HRD[1] ~ dbeta(1,20)   # Harvest rate of smolts in driftnet fishery
  HRCTN[1,1]  ~dbeta (1,20)	# Harvest rate smolts in coastal trapnet fishery
  HRCGN[1,1] ~dbeta (1,20)	# Harvest rate smolts in coastal gillnet fishery
  
  qlW1~ dlnorm(mql[1], tauqd)	# Catchability of reproductive salmon in the longline fishery
  qlR1~ dlnorm(mql[1], tauqd)	# Catchability of non-reproductive salmon in the longline 
  for(i in 1:m){
    qlW[i,1]<-qlW1
    qlR[i,1]<-qlR1
  }
  mql[1] <- log(-log(1-HRL[1])/El[1,1])	# Mean catchability coefficient longline fishery	
  
  qdW1 ~ dlnorm(mqd[1], tauqd)	# Catchability of reproductive salmon in the driftnet 	
  qdR1 ~ dlnorm(mqd[1], tauqd)	# Catchability of non-reproductive salmon in the driftnet 
  for(i in 1:m){
    qdW[i,1]<-qdW1
    qdR[i,1]<-qdR1
  }
  # qdW[1] ~ dlnorm(mqd[1], tauqd)	# Catchability of reproductive salmon in the driftnet 	
  # qdR[1] ~ dlnorm(mqd[1], tauqd)	# Catchability of non-reproductive salmon in the driftnet 
  mqd[1] <- log(-log(1-HRD[1])/Edo[1,1])	# Mean catchability coefficient driftnet fishery	
  
  # Catchability coefficient of salmon by coastal fisheries		
  qctnW[1,1] ~ dlnorm(mqctn[1,1], tauqctn)	# Catchability wild salmon coastal trapnet
  qctnR[1,1] ~ dlnorm(mqctn[1,1], tauqctn)	# Catchability reared salmon coastal trapnet
  mqctn[1,1] <- log(-log(1-HRCTN[1,1])/Ectn[1,1,1])# Mean catchability coastal trapnet fishery
  
  qcgnW[1,1] ~ dlnorm(mqcgn[1,1], tauqcgn)	# Catchability wild salmon coastal gillnet	
  qcgnR[1,1] ~ dlnorm(mqcgn[1,1], tauqcgn)	# Catchability reared salmon coastal gillnet 
  mqcgn[1,1] <- log(-log(1-HRCGN[1,1])/Ecgn[1,1,1])# Mean catchability coastal gillnet fishery
  
  for (au in 2:3){
    
    HRCTN[1,au]  ~dbeta (1,20)	# Harvest rate smolts in coastal trapnet fishery
    HRCGN[1,au] ~dbeta (1,20)	# Harvest rate smolts in coastal gillnet fishery
    
    qctnR[1,au] <- qctnR[1,1]
    mqctn[1,au] <- mqctn[1,1]		
    qctnW[1,au] <- (qctnW[1,1] / qctnR[1,1] ) * qctnR[1,au]
    
    qcgnR[1,au] <- qcgnR[1,1]
    mqcgn[1,au] <- mqcgn[1,1]			
    qcgnW[1,au] <- (qcgnW[1,1] / qcgnR[1,1] ) * qcgnR[1,au]
  }
  
  # Grilse (j=2) and 2SW salmon	(j=3)
  
  
  #qlW[2] ~ dlnorm(mql[2], tauql)	# Catchability of reproductive salmon in the longline fishery
  #qlR[2] ~ dlnorm(mql[2], tauql)	# Catchability of non-reproductive salmon in the longline 
  #mql[2] <- log(-log(1-HRL[2])/El[1,1])	# Mean catchability coefficient longline fishery	
  
  # Autocorrelation model for 2SW (==MSW) catchability in longline fishery
  # =======================================================================
  # Note that autocorrelation takes place between years but in qlW[i,j] index i is smolt cohort 
  for(i in 1:(m+3)){
    
    # i+1=i+2-1: cohort i age 2 transformed to year i+1
    logit(qlW[i,2])<-logit_qlW[i+1]
    logit(qlR[i,2])<-logit_qlW[i+1]#logit_qlR[i+1]
    
    # mean reverting AR(1) (index i is year)
    logit_qlW[i+1]~dnorm(mu_qlW[i+1], tau_ql)
    mu_qlW[i+1]<-phi_ql*logit_qlW[i]+(1-phi_ql)*mean_qlW # mean reverting AR(1)
    
    #logit_qlR[i+1]~dnorm(mu_qlR[i+1], tau_ql)
    #mu_qlR[i+1]<-phi_ql*logit_qlR[i]+(1-phi_ql)*mean_qlR # mean reverting AR(1)
    
  }
  logit_qlW[1]~dnorm(0,0.39) # implies uniform[0,1] prior for initial catchability
  #logit_qlR[1]~dnorm(0,0.39) # implies uniform[0,1] prior for initial catchability
  
  for(i in 1:(m)){
    for(j in 3:5){
      # Same catchability of MSW salmon in the longline fishery
      qlW[i,j] <- qlW[i+j-2,2]	
      qlR[i,j] <-qlW[i,j]# qlR[i+j-2,2]
    }
  }
  
  mean_qlW~dnorm(0,0.39) # implies uniform[0,1] prior for mean catchability

  #mean_qlR~dnorm(0,0.39) # implies uniform[0,1] prior for mean catchability
  tau_ql<-1/((1-pow(phi_ql,2))*(sd_ql*sd_ql))  
  # Marginal variance chosen to give uniform when mean=0, otherwise unimodal
  sd_ql~dunif(0.01,1.6)
  phi_ql~dunif(0.01,0.99) #positive autocorrelation
  
  
  # Mean reverting AR(1)-model for catchability in driftnet fisheries (both offshore and coastal)
  # =======================================================================
  
  for(j in 2:3){
    for(i in 1:(m+3)){
      # i+j-1: cohort i age j transformed to year
      logit(qdW[i,j])<-logit_qdW[i+j-1,j]
      logit(qdR[i,j])<-logit_qdW[i+j-1,j]#logit_qdR[i+j-1,j]
    }
    
    for(i in 1:(m+4)){
      # mean reverting AR(1) (index i is year)
      logit_qdW[i+1,j]~dnorm(mu_qdW[i+1,j], tau_qd)
      mu_qdW[i+1,j]<-phi_qd*logit_qdW[i,j]+(1-phi_qd)*mean_qdW[j]
      
      #logit_qdR[i+1,j]~dnorm(mu_qdR[i+1,j], tau_qd)
      #mu_qdR[i+1,j]<-phi_qd*logit_qdR[i,j]+(1-phi_qd)*mean_qdR[j]
    }
    
    # implies uniform[0,1] prior for initial catchability
    logit_qdW[1,j]~dnorm(0,0.39)

    #logit_qdR[1,j]~dnorm(0,0.39)
  }
  mean_qdW[2]~dnorm(0,0.39) # 2SW

  #mean_qdW[1]<-mean_qdW[2]*eff_qd[1] # 1SW # this is dealt with in old way!
  mean_qdW[3]<-mean_qdW[2]*eff_qd # MSW
  ##mean_qdR[3]~dnorm(0,0.39) # 2SW
  ##mean_qdR[2]<-mean_qdR[3]*eff_qd[1] # 1SW
  ##mean_qdR[4]<-mean_qdR[3]*eff_qd[2] # MSW
  
  eff_qd~dbeta(10,5)
  # for(i in 1:2){
  #   eff_qd[i]~dbeta(10,5)
  # }
  
  tau_qd<-1/((1-pow(phi_qd,2))*(sd_qd*sd_qd))  
  # Marginal variance chosen to give uniform when mean=0, otherwise unimodal
  sd_qd~dunif(0.01,1.6)
  phi_qd~dunif(0.01,0.99) #positive autocorrelation
  
  for(i in 1:(m)){
    for(j in 4:6){
      # Same catchability of MSW salmon in the longline fishery
      qdW[i,j] <- qdW[i+j-3,3]	
      qdR[i,j] <-qdW[i,j]# qdR[i+j-4,4]	
      
    }
  }
  
  # River fishery
  # ================
  
  qrR[2] <- -log(1-HRR[2])/Er[1,1]	# Catchability of non-reproductive salmon in the river	
  qrR[3] <-qrR[2]
  
  for (j in 2:3){ 
    
    HRR[j] ~dbeta (5,1)	# Harvest rate in terminal river fishery	  
    
    #    	qdW[j] ~ dlnorm(mqd[j], tauqd)	# Catchability of reproductive salmon in the driftnet 	
    # 		qdR[j] ~ dlnorm(mqd[j], tauqd)	# Catchability of non-reproductive salmon in the driftnet 
    # 		mqd[j] <- log(-log(1-HRD[j])/Edo[1,1])	# Mean catchability coefficient driftnet fishery	
    #     
    # Catchability coefficient of salmon from area 1 by coastal fishery 
    # This data is available only for AU1
    qctnW[j,1] ~ dlnorm(mqctn[j,1], tauqctn)	# Catchability wild salmon coastal trapnet
    qctnR[j,1] ~ dlnorm(mqctn[j,1], tauqctn)	# Catchability reared salmon coastal trapnet
    mqctn[j,1] <- log(-log(1-HRCTN[j,1])/Ectn[1,1,1])	# Mean catchability coastal trapnet 
    
    qcgnW[j,1] ~ dlnorm(mqcgn[j,1], tauqcgn)	# Catchability wild salmon coastal gillnet
    qcgnR[j,1] ~ dlnorm(mqcgn[j,1], tauqcgn)	# Catchability reared salmon coastal gillnet 
    mqcgn[j,1] <- log(-log(1-HRCGN[j,1])/Ecgn[1,1,1])# Mean catchability coastal gillnet 
    
    
    # For AU2 and AU3 we assume the ratio between wild and reared to be the same as in AU1
    for (au in 2:3){
      
      qctnR[j,au] ~ dlnorm(mqctn[j,au], tauqctn)	# Catchability reared salmon coastal trapnet
      mqctn[j,au] <- log(-log(1-HRCTN[j,au])/Ectn[1,1,1])	# Mean catchability coastal trapnet 
      qctnW[j,au] <- (qctnW[j,1] / qctnR[j,1] ) * qctnR[j,au]	
      
      qcgnR[j,au] ~ dlnorm(mqcgn[j,au], tauqcgn)	# Catchability reared salmon coastal gillnet 
      mqcgn[j,au] <- log(-log(1-HRCGN[j,au])/Ecgn[1,1,1])# Mean catchability coastal gillnet 
      qcgnW[j,au] <- (qcgnW[j,1] / qcgnR[j,1] ) * qcgnR[j,au]
      
    }
    
    #  HRD[j] ~dbeta (2,5)	# Harvest rate driftnet fishery	
    #  HRL[j] ~dbeta (2,5)	# Harvest rate longline fishery	
  }
  for (au in 1:3){
    HRCGN[2,au]  ~dbeta (2,2.5) 	# Harvest rate coastal gillnet fishery of fish from area k
    HRCTN[2,au]  ~dbeta (2,2.5) 	# Harvest rate coastal trapnet fishery of fish from area k	
    
    HRCGN[3,au]  ~dbeta (2,5) 	# Harvest rate coastal gillnet fishery of fish from area k
    HRCTN[3,au]  ~dbeta (2,5) 	# Harvest rate coastal trapnet fishery of fish from area k	
  }
  for (j in 4:5){ 
    
    HRR[j] ~dbeta (5,1)	# Harvest rate in terminal river fishery	
    qrR[j] <- qrR[3]	# Same catchability of MSW salmon in the river
    
    #qlW[j] <- qlW[2]	# Same catchability of MSW salmon in the longline fishery
    #qlR[j] <- qlR[2]	# Same catchability of MSW salmon in the longline fishery
    
    #  	qdW[j] ~ dlnorm(mqd[j], tauqd)	# Catchability of reproductive salmon in the driftnet 	
    # 	qdR[j] ~ dlnorm(mqd[j], tauqd)	# Catchability of non-reproductive salmon in the driftnet 
    # 	mqd[j] <- log(-log(1-HRD[j])/Edo[1,1])	# Mean catchability coefficient driftnet fishery	
    
    qctnW[j,1] <- qctnW[3,1]# Same catchability of MSW salmon in the coastal trapnet fishery
    qcgnW[j,1] <- qcgnW[3,1]# Same catchability of MSW salmon in the coastal gillnet fishery
    qctnR[j,1] <- qctnR[3,1]# Same catchability of MSW salmon in the coastal trapnet fishery
    qcgnR[j,1] <- qcgnR[3,1]# Same catchability of MSW salmon in the coastal gillnet fishery
    
    HRCGN[j,1]  ~dbeta (2,5) 	# Harvest rate coastal gillnet fishery of fish from area k
    HRCTN[j,1]  ~dbeta (2,5) 	# Harvest rate coastal trapnet fishery of fish from area k	
    
    for (au in 2:3){
      qctnW[j,au] <- (qctnW[j,1] / qctnR[j,1] ) * qctnR[j,au]	
      qcgnW[j,au] <- (qcgnW[j,1] / qcgnR[j,1] ) * qcgnR[j,au]	
      qctnR[j,au] <- qctnR[3,au]# Same catchability of MSW salmon in the coastal trapnet fishery
      qcgnR[j,au] <- qcgnR[3,au]# Same catchability of MSW salmon in the coastal gillnet fishery
      
      HRCGN[j,au]  ~dbeta (2,5) 	# Harvest rate coastal gillnet fishery of fish from area k
      HRCTN[j,au]  ~dbeta (2,5) 	# Harvest rate coastal trapnet fishery of fish from area k	
      
    }
    #  HRD[j] ~dbeta (2,5)	# Harvest rate driftnet fishery	
    #	HRL[j] ~dbeta (2,5)	# Harvest rate longline fishery	
  }
  
  for (j in 6:6){ 
    
    HRR[j] ~dbeta (5,1)	# Harvest rate in terminal river fishery
    qrR[j] <- qrR[3]	# Catchability of non-reproductive salmon in the river	
    
    for(i in 1:(m+3)){
      qlW[i,j] <- 0 #not defined in bugs model
      qlR[i,j] <- 0 
    }
    
    #	qdW[j] ~ dlnorm(mqd[j], tauqd)	# Catchability of reproductive salmon in the driftnet 	
    #	qdR[j] ~ dlnorm(mqd[j], tauqd)	# Catchability of non-reproductive salmon in the driftnet 
    #	mqd[j] <- log(-log(1-HRD[j])/Edo[1,1])	# Mean catchability coefficient driftnet fishery	
    
    # Catchability coefficient of salmon from area l by coastal fishery		
    qctnW[j,1] ~ dlnorm(mqctn[j,1], tauqctn)	# Catchability wild salmon coastal trapnet
    qctnR[j,1] ~ dlnorm(mqctn[j,1], tauqctn)	# Catchability reared salmon coastal trapnet
    mqctn[j,1] <- log(-log(1-HRCTN[j,1])/Ectn[1,1,1])# Mean catchability coastal trapnet fishery	
    
    qcgnW[j,1] ~ dlnorm(mqcgn[j,1], tauqcgn)	# Catchability wild salmon coastal gillnet		
    qcgnR[j,1] ~ dlnorm(mqcgn[j,1], tauqcgn)	# Catchability reared salmon coastal gillnet 
    mqcgn[j,1] <- log(-log(1-HRCGN[j,1])/Ecgn[1,1,1])	# Mean catchability coastal gillnet fishery	
    
    HRCGN[j,1]  ~dbeta (2,5) 	# Harvest rate coastal gillnet fishery of fish from area k
    HRCTN[j,1]  ~dbeta (2,5) 	# Harvest rate coastal trapnet fishery of fish from area k	
    
    for (au in 2:3){
      qctnW[j,au] <- (qctnW[j,1] / qctnR[j,1] ) * qctnR[j,au]	
      qcgnW[j,au] <- (qcgnW[j,1] / qcgnR[j,1] ) * qcgnR[j,au]
      qctnR[j,au] ~ dlnorm(mqctn[j,au], tauqctn)	# Catchability reared salmon coastal trapnet
      mqctn[j,au] <- log(-log(1-HRCTN[j,au])/Ectn[1,1,1])# Mean catchability coastal trapnet fishery	
      qcgnR[j,au] ~ dlnorm(mqcgn[j,au], tauqcgn)	# Catchability reared salmon coastal gillnet 
      mqcgn[j,au] <- log(-log(1-HRCGN[j,au])/Ecgn[1,1,1])	# Mean catchability coastal gillnet fishery	
      
      HRCGN[j,au]  ~dbeta (2,5) 	# Harvest rate coastal gillnet fishery of fish from area k
      HRCTN[j,au]  ~dbeta (2,5) 	# Harvest rate coastal trapnet fishery of fish from area k	
    }
    
    #	HRD[j] ~dbeta (2,5)	# Harvest rate driftnet fishery	
    #	HRL[j] ~dbeta (2,5)	# Harvest rate longline fishery	
  }
  
  tauqd ~ T(dgamma(20,1),10,)
  tauql ~ T(dgamma(50,1),10,)
  tauqr ~ T(dgamma(20,1),10,)
  tauqcgn ~ T(dgamma(20,1),10,)			
  tauqctn ~ T(dgamma(20,1),10,)
  
  
  #Negbin ~1/overdispersion parameters
  for (j in 1:6){
    
    rrW[j] ~ dunif(1,200)
    rrRsp[j] ~  dunif(1,200)			
    rrR[j] ~  dunif(1,200)   
    rcW[j] ~  dunif(1,200)
    rcR[j] ~  dunif(1,200)		
    rdW[j]~ dunif(1,200)
    rdR[j] ~ dunif(1,200)		
    rlW[j] ~ dunif(1,200)									
    rlR[j] ~ dunif(1,200)								
    
  }	
  
  # Tag retention rate
  Tretain ~ T(dbeta(20,8),0.5,1)
  
  # Reporting rates in the different fisheries
  reportc ~ T(dbeta(11,9),0.2,0.8)	# coastal fishery
  reportrR ~ T(dbeta(16,6),0.3,0.95)	# river fishery
  reportrW ~ T(dbeta(16,6),0.3,0.95)	# river fishery	
  reportd ~ T(dbeta(8,4),0.2,0.95)	# driftnet fishery
  reportl ~ T(dbeta(10,4),0.3,0.95)	# longline fishery
  
  #Catch likelihoods	
  tauCRW<-1/(SCRW*SCRW)	
  tauCR<-1/(SCR*SCR)
  tauCC<-1/(SCC*SCC)
  tauCO<-1/(SCO*SCO)
  tauCT<-1/(SCT*SCT)				  
  
  SCRW~T(dbeta(1.45,8.19),,0.50)	
  SCR~T(dbeta(1.45,8.19),,0.50)
  SCC~T(dbeta(1.45,8.19),,0.50) 
  SCO~T(dbeta(1.45,8.19),,0.50) 
  SCT~T(dbeta(1.45,8.19),,0.50)					 
  
  # These CVs should be at the original scale, likely to be close to the values of SCR-SCO
  cvCRW<-sqrt(exp(SCRW*SCRW)-1)
  cvCR<-sqrt(exp(SCR*SCR)-1)
  cvCC<-sqrt(exp(SCC*SCC)-1)
  cvCO<-sqrt(exp(SCO*SCO)-1)
  cvCT<-sqrt(exp(SCT*SCT)-1)						  
  
  #Maturation
  
  for (i in 1:(m+5)){ # calendar years 1987-present+5
  p.rel[i]~dbeta(alpha_rel[i],beta_rel[i])   #i is calendar year
    for (j in 1:4){ # sea ages (1=1SW, 2=2SW etc.)
      
      muLW[i,j]<-cL[i]+bL[j]+delta[j]*Temp[i] # temperature effect differs on age groups
      #lw[i,j]~dnorm(muLW[i,j],tauL[j])
      lw[i,j]<-muLW[i,j]+(1/sqrt(tauL[j]))*eLW[i]
      logit(LW[i,j])<-lw[i,j]
      
      muLR[i,j]<-cL[i]+bL[j]+LReffect[j]+delta[j]*Temp[i]
      #lr[i,j]~dnorm(muLR[i,j],tauL[j])
      lr[i,j]<-muLR[i,j]+(1/sqrt(tauL[j]))*eLR[i]
      logit(LR[i,j])<-lr[i,j]
      
    }	
    eLW[i]~dnorm(0,1)
    eLR[i]~dnorm(0,1)
    Temp[i]~dnorm(muTemp[i],tauTemp[i])     #muTemp and tauTemp from Temperature data file
    cL[i]~dnorm(mucL,taucL)
    
    LW[i,5]<-1
    LR[i,5]<-1
    
    LW[i,6]<-1   #j=6 not  used as last year of spawners (j=6) use mat[,5] 
    LR[i,6]<-1
  }
  for (j in 1:3){ # sea ages (1=1SW, 2=2SW etc.)
    delta[j]~dunif(0.001,1)
    LReffect[j]~dlnorm(-1,2)
  }
  # MSW effects
  delta[4]<-delta[3]
  LReffect[4]<-LReffect[3]
  
  mucL~dnorm(0.4,5.4)
  taucL~dlnorm(4.6,2.8)
  
  bL[1]~dnorm(-2.9,5.4)
  bL[2]~dnorm(-0.84,5.4)
  bL[3]~dnorm(0.047,5.4)
  bL[4]~dnorm(1.40,5.4)
  
  tauL[1]~dlnorm(0.42,49)
  tauL[2]~dlnorm(1.7,44)
  tauL[3]~dlnorm(2.3,41)
  tauL[4]~dlnorm(1.4,46)
  
  #priors for spawner counting
  for(s in 1:stocks){
    
    
     logit_mu_spawn[s]~T(dnorm(mu_mu_sp[s],tau_mu_sp[s]),pdl[s],pdh[s]) 
    logit_CV_spawn[s]~T(dnorm(mu_CV_sp[s],tau_CV_sp[s]),pcvl[s],pcvh[s])  
    
    mu_spawn[s]<-ilogit(logit_mu_spawn[s]) 
    CV_spawn[s]<-ilogit(logit_CV_spawn[s]) 
    
    tau_spawn[s]<-(1-mu_spawn[s])^2/CV_spawn[s]^2

    
    for(i in 1:(m+5)){
      logit_pdetect[i,s]~T(dnorm(logit_mu_spawn[s],tau_spawn[s]),-7,7)
      p.detect[i,s]<-ilogit(logit_pdetect[i,s])
      
      p.ladder[i,s]~dbeta(alpha_ladder[i,s],beta_ladder[i,s])
      surv_migr[i,s]~dbeta(alpha_migr[i,s],beta_migr[i,s])   
  
    }
    
    eta_msw[s]<-(1/corr_msw[s])-1
    corr_msw[s]~dunif(0.0001,0.5)
    
    
  }
  
  tauDS<-1/(log(cvDS*cvDS+1))
  cvDS~dlnorm(-2.37,8)
  
  coefDS<-coefDS_tmp+1   # To ensure Simo count is always overestimation  
  coefDS_tmp~dlnorm(-2.029014,1.1211)
  #coefDS<-1.05 # assume that Simojoki Didson count is underestimate 
  
  for(rs in 1:rstocks){
    aTrap[rs]<-muTrap[rs]*etaTrap[rs]+1
    bTrap[rs]<-(1-muTrap[rs])*etaTrap[rs]+1
  }
  
  muTrap[1]~T(dbeta(2,2),0.01,0.99)        #Lule?lven AU2
  muTrap[2]~dbeta(72.7,197)   #Dal?lven AU3
  
  #etaTrap[1]~dlnorm(10,0.1)
  etaTrap[1]~dlnorm(10,1)
  etaTrap[2]~dlnorm(3.7,15.7)
  
  Usmolt~T(dlnorm(0.01, 85),0.8,1.5) 			# Uncertainty factor
  p.mort~dbeta(4,12)
  
})

get_tau <- nimble::nimbleFunction(
  run = function(N = double(0),p=double(0)) {
    returnType(double(0))
    mean<-N*p
    var<-N*p*(1-p)
    cv<-(sqrt(var))/mean
    s<-sqrt(log(cv*cv+1))
    tau<-min(max(1/(s*s),4),400)
    #tau<-max(1/(s*s),20) 
    return(tau)
  })
get_mu <- nimble::nimbleFunction(
  run = function(N = double(0),p=double(0)) {
    returnType(double(0))
    mean<-N*p
    var<-N*p*(1-p)
    cv<-(sqrt(var))/mean
    s<-sqrt(log(cv*cv+1))
    mu<-log(mean/exp(0.5*s^2))
    
    return(mu)
  })
#
#get_mu_tau <- nimble::nimbleFunction(           #LogN
#  run = function(N = double(0),p=double(0)) {
#    returnType(double(1))
#    mean<-N*p
#    var<-N*p*(1-p)
#    cv<-(sqrt(var))/mean
#    s<-sqrt(log(cv*cv+1))
#    tau<-min(max(1/(s*s),4),400)
#    #tau<-max(1/(s*s),20)
#    snew<-sqrt(1/tau)
#    mu<-log(mean/exp(0.5*snew^2))
#    return(c(mu,tau))
#  })

 
 #   HrWAU[i,j,1] <-  getmean(HrW[i,j,1:stocks],au1_stocks[1:nAU1])
getmean<- nimble::nimbleFunction(          
  run = function(hr = double(1),ind=double(1)) {
  returnType(double(0))
  return(mean(hr[ind]))
  })

get_mu_tau <- nimble::nimbleFunction(           #Normal
  run = function(N = double(0),p=double(0)) {
    returnType(double(1))
    bmean<-N*p
    bvar<-N*p*(1-p)      
    btau<-1/bvar               
  
    return(c(bmean,btau))
  })

get.err_par <- nimble::nimbleFunction(   #function to return bounds of Uniform dist for process errors
  
  #survtmp vector of reared wild values for a given year, age, fleet
  run = function(survtmp=double(1),zzv=double(0),nst=integer(0),monE=double(0)) {    
    returnType(double(1)) 
    kk<-numeric(length = nst)
    err_par<-numeric(length = 2)
    
    for(st in 1:nst){
      if(survtmp[st]>=0.5){
        kk[st]<-1
      }else{
        kk[st]<-0
      }
    }
    
    maxe<-((1/survtmp)*kk) + (2 * (1 - kk))
    mine<- (1 - (maxe - 1)) 
    vve<- (maxe - mine)^2 /12   # variance depending on survival (vector) 
    minv<-min(vve)      #compare stocks/find min
    
    ve<-monE*zzv/12    #variance depending on time step
    
    err_par[2] <- 1 + sqrt(3 * min(ve, minv))   #variance depending on time step
    err_par[1] <- 1 - sqrt(3 * min(ve, minv)) 
    
    return(err_par)
  })



get.err_par_riv <- nimble::nimbleFunction(   #function to return bounds of Uniform dist for process errors
  
  #survtmp vector of reared wild values for a given year, age, fleet
  run = function(survw=double(1),survr=double(0),zzv=double(0),nstw=integer(0),monE=double(0)) {    
    returnType(double(1)) 
    kkw<-numeric(length = nstw)
    err_par<-numeric(length = 2)
    
    for(st in 1:nstw){          #ifelse?
      if(survw[st]>=0.5){
        kkw[st]<-1
      }else{
        kkw[st]<-0
      }
    }
    
    
    if(survr>=0.5){
      kkr<-1
    }else{
      kkr<-0
    }
    
    
    maxew<-((1/survw)*kkw) + (2 * (1 - kkw))
    minew<- (1 - (maxew - 1)) 
    vvew<- (maxew - minew)^2 /12   # variance depending on survival (vector) 
    
    maxer<-((1/survr)*kkr) + (2 * (1 - kkr))
    miner<- (1 - (maxer - 1)) 
    vver<- (maxer - miner)^2 /12   # variance depending on survival (vector) 
    
    minv<-min(c(vvew,vver))      #compare stocks/find min
    
    ve<-monE*zzv/12    #variance depending on time step
    
    err_par[2] <- 1 + sqrt(3 * min(ve, minv))   #variance depending on time step
    err_par[1] <- 1 - sqrt(3 * min(ve, minv)) 
    
    return(err_par)
  })


get.err_par_coast <- nimble::nimbleFunction(   #function to return bounds of Uniform dist for process errors
  
  #survtmp vector of reared wild values for a given year, age, fleet
  run = function(survtmp=double(2),zzv=double(0),nau=integer(0),monE=double(0)) {    
    returnType(double(1)) 
    kk<-array(0,dim=c(2,nau))
    maxe<-array(0,dim=c(2,nau))
    mine<-array(0,dim=c(2,nau))
    vve<-array(0,dim=c(2,nau))
    
    err_par<-numeric(length = 2)
    
    for(ij in 1:2){ 
      for(au in 1:nau){                   
        if(survtmp[ij,au]>=0.5){
          kk[ij,au]<-1
        }else{
          kk[ij,au]<-0
        }
      }
    }
    
    for(ij in 1:2){
      for(au in 1:nau){   
        maxe[ij,au]<-((1/survtmp[ij,au])*kk[ij,au]) + (2 * (1 - kk[ij,au]))
        mine[ij,au]<- (1 - (maxe[ij,au] - 1)) 
        vve[ij,au]<- (maxe[ij,au] - mine[ij,au])^2 /12   # variance depending on survival (vector) 
      }
    }
    
    minv<-min(vve)      #compare stocks/find min
    
    ve<-monE*zzv/12    #variance depending on time step
    
    err_par[2] <- 1 + sqrt(3 * min(ve, minv))   #variance depending on time step
    err_par[1] <- 1 - sqrt(3 * min(ve, minv)) 
    
    return(err_par)
  })
  
get_EPR <- nimble::nimbleFunction(   
  
  run = function(mat=double(2),Msmolt=double(1),Mad=double(0),nages=integer(0),nyr=integer(0),Mseal=double(2),Mmigr=double(1),pl=double(1),fec=double(1),pfem=double(2)){ 
  
  returnType(double(1))
  pmat<-nimArray(0,dim=c(nyr,nages))
  pimm<-nimArray(0,dim=c(nyr,nages))
  smat<-nimArray(0,dim=c(nyr,nages))
  simm<-nimArray(0,dim=c(nyr,nages))
  epr<-nimNumeric(nyr)
  
  pmat[1,1]<-0                                                 
  pimm[1,1]<-1
  
  for(j in 2:nages){
    pmat[1,j]<-mat[1,(j-1)]*pimm[1,(j-1)]      
    pimm[1,j]<-(1-mat[1,(j-1)])*pimm[1,(j-1)]
    
  }
  
  for(i in 2:nyr){
    pmat[i,1]<-0
    pimm[i,1]<-1
    
    for(j in 2:6){
      pmat[i,j]<-mat[i,(j-1)]*pimm[(i-1),(j-1)]     
      pimm[i,j]<-(1-mat[i,(j-1)])*pimm[(i-1),(j-1)]
      
    }
  }
  
    simm[1,1]<-exp(-(11*Msmolt[1]/12))*exp(-(Msmolt[1]*Mseal[1,1]/12))   
    smat[1,1]<-1    #not used anywhere
    
    for(j in 2:nages){
      simm[1,j]<-exp(-(12*Mad/12))
      smat[1,j]<-exp(-(3*Mad/12))*exp(-(2*Mad*Mseal[1,j]/12))*pl[(1+j-1)]*Mmigr[(1+j-1)]
    }
    
    for(i in 2:nyr){
      
      simm[i,1]<-exp(-(11*Msmolt[i]/12))*exp(-(Msmolt[i]*Mseal[i,1]/12))   
      smat[i,1]<-1    #not used anywhere
      
      for(j in 2:nages){
        simm[i,j]<-exp(-(12*Mad/12))
        smat[i,j]<-exp(-(3*Mad/12))*exp(-(2*Mad*Mseal[i,j]/12))*pl[(i+j-1)]*Mmigr[(i+j-1)]
      }
    }
    
    for(i in 1:5){
      epr[i]<-0 
    }
    for(i in 6:nyr){ 
    epr[i]<-pmat[i,2]*pfem[(i-1),1]*fec[1]*simm[(i-1),1]*smat[(i-1),2]+pmat[i,3]*pfem[(i-2),2]*fec[2]*simm[(i-2),1]*simm[(i-1),2]*smat[(i-2),3]+pmat[i,4]*pfem[(i-3),3]*fec[3]*simm[(i-3),1]*simm[(i-2),2]*simm[(i-1),3]*smat[(i-3),4]+pmat[i,5]*pfem[(i-4),4]*fec[4]*simm[(i-4),1]*simm[(i-3),2]*simm[(i-2),3]*simm[(i-1),4]*smat[(i-4),5]+pmat[i,6]*pfem[(i-5),5]*fec[5]*simm[(i-5),1]*simm[(i-4),2]*simm[(i-3),3]*simm[(i-2),4]*simm[(i-1),5]*smat[(i-5),6]     #5SW
  }
  
  return(epr)

})


assign('getmean', getmean, envir = .GlobalEnv)
assign('get_mu', get_mu, envir = .GlobalEnv)
assign('get_tau', get_tau, envir = .GlobalEnv)
assign('get_mu_tau', get_mu_tau, envir = .GlobalEnv)
assign('get.err_par', get.err_par, envir = .GlobalEnv)
assign('get.err_par_riv', get.err_par_riv, envir = .GlobalEnv)
assign('get.err_par_coast', get.err_par_coast, envir = .GlobalEnv)
assign('get_EPR', get_EPR, envir = .GlobalEnv)

