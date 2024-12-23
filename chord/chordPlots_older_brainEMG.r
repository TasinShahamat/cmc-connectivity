### Created using Steve Peterson's code as a template
library(R.matlab)

pathname <- "C:\\Users\\tasin\\OneDrive - University of Central Florida\\bckup-asus\\BRaIN Lab\\CMC\\Data\\connMatChord\\"
# Define conditions and frequency bands
conds <- c("LEI", "LME")  
freqBands <- c("theta", "alpha")
n_components <- 9  # young - 11, old - 9

# Loop through conditions and frequency bands
for (cond in conds) {
  for (plotBand_cort  in freqBands) {
    netData_sbjs = readMat(paste0(pathname, "netVals_old_", cond, "_early_epoch_PPT_aveTime_sbjs.mat"))
    
    head(netData_sbjs$netVals[1])
    # df = as.data.frame(netData_sbjs$netVals)
    
    #Determine significant connections above 0
    pvals_theta = matrix(0L, nrow = 1, ncol = n_components*n_components-n_components)
    pvals_alpha = matrix(0L, nrow = 1, ncol = n_components*n_components-n_components)
    counter = 0
    counter_all = 0
    for(i in 1:n_components){
      for(j in 1:n_components){
        counter_all = counter_all+1
        if(i!=j){
          counter=counter+1
          pvals_theta[counter] = t.test(netData_sbjs$netVals[[1]][[counter_all]][[1]])$p.value
          pvals_alpha[counter] = t.test(netData_sbjs$netVals[[2]][[counter_all]][[1]])$p.value
        }
      }
    }
    
    #Reshape into n_components x n_components
    pvals_theta_box = matrix(1L, nrow = n_components, ncol = n_components)
    pvals_alpha_box = matrix(1L, nrow = n_components, ncol = n_components)
    counter = 0
    for(i in 1:n_components){
      for(j in 1:n_components){
        if(i!=j){
          counter=counter+1
          pvals_theta_box[j,i] = pvals_theta[counter]
          pvals_alpha_box[j,i] = pvals_alpha[counter]
        }
      }
    }
    
    #FDR correction
    # for(i in 1:n_components){
    #   pvals_theta_box[,i]=p.adjust(pvals_theta_box[,i],method='fdr')
    #   pvals_alpha_box[,i]=p.adjust(pvals_alpha_box[,i],method='fdr')
    # }
    pvals_alpha_box[3,4] = 1  # temporary fix. there was a nan value due to all array values being 0 in that connection
    
    #Check for significance
    isSig_theta_box = matrix(0L, nrow = n_components, ncol = n_components)
    isSig_alpha_box = matrix(0L, nrow = n_components, ncol = n_components)
    counter = 0
    for(i in 1:n_components){
      for(j in 1:n_components){
        if(i!=j){
          counter=counter+1
          if(pvals_theta_box[i,j]<0.05){
            isSig_theta_box[i,j]=1
          }
          if(pvals_alpha_box[i,j]<0.05){
            isSig_alpha_box[i,j]=1
          }
        }
      }
    }
    
    #Load average data and mask it by significance
    netData_ave = readMat(paste0(pathname, "netVals_old_", cond, "_early_epoch_PPT_aveTime.mat"))
    thetaAve = netData_ave$netVals.ave[[1]]
    alphaAve = netData_ave$netVals.ave[[2]]
    thetaAve[isSig_theta_box==0]=0
    alphaAve[isSig_alpha_box==0]=0
    
    #Chord plot (for cortico-cortical data)
    library(circlize)
    circos.clear()
    
    #Rearrange order in matrix
    exchangeInds1 = c(1,2,3) #c(5,2,6,1,3,8,4,7)
    exchangeInds2 = c(4,5,6,7,8,9)
    exchangeInds=c(exchangeInds1,exchangeInds2)
    thetaAve_orig = thetaAve
    alphaAve_orig = alphaAve
    for(i in 1:n_components){
      for(j in 1:n_components){
        thetaAve[i,j] = thetaAve_orig[exchangeInds[i],exchangeInds[j]]
        alphaAve[i,j] = alphaAve_orig[exchangeInds[i],exchangeInds[j]]
      }
    }
    
    cortTheta = 0*thetaAve
    cortAlpha = 0*alphaAve
    cortTheta[4:n_components,1:3] = thetaAve[4:n_components,1:3]
    cortAlpha[4:n_components,1:3] = alphaAve[4:n_components,1:3]
    cortTheta[1:3,4:n_components] = thetaAve[1:3,4:n_components]
    cortAlpha[1:3,4:n_components] = alphaAve[1:3,4:n_components]
    
    
    if(plotBand_cort=='theta'){
      dat_in = cortTheta
    }else{
      dat_in = cortAlpha
    }
    png(paste(pathname,cond,"_",plotBand_cort,"_brainEMG_old.png",sep=""),bg="transparent",width = 4160, height = 4160)
    # png(paste(pathname,cond,"_",plotBand_cort,"_brainEMG_young.svg",sep=""),bg="transparent",width = 4160, height = 4160)
    color_vals = c('gold','cyan','green','purple','orange',
                   'cyan','gold','purple','orange')
    color_vals_links = adjustcolor(color_vals, alpha.f = 0.6)
    factors = c("Motor","PPC_r","PPC_l","LTA","LSO","LRF","LST","LAD","LPD")
    circos.par(cell.padding = c(0, 0, 0, 0),
                start.degree = -35+180,
                canvas.xlim = c(-1.5, 1.5),  # Expand plot horizontally
                canvas.ylim = c(-1.5, 1.5),  # Expand plot vertically
                gap.degree=c(0,0,20,0,0,0,0,0,20)) 
    sign_vals = sign(dat_in)
    A = abs(dat_in)*10^4
    circos.initialize(factors, xlim = cbind(c(0,0,0,0,0,0,0,0,0), 2*c(1,1,1,1,1,1,1,1,1))) 
    
    
    circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05,
                           #panel.fun for each sector
                           panel.fun = function(x, y) {
                             offset = 0
                             #select details of current sector
                             name = get.cell.meta.data("sector.index")
                             i = get.cell.meta.data("sector.numeric.index")
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             
                             #text direction (dd) and adjusmtents (aa)
                             theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                             dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                             aa = c(1, 0.5)
                             if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                             
                             #plot labels
                             circos.text(x=mean(xlim), y = ylim[2] + 3, labels=name, facing = "bending.inside", cex=10, col="black",  adj = aa)
                             
                             #plot main sector
                             circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                         col = color_vals[i], border=color_vals[i])
                             
                             #plot axis
                            #  circos.axis(labels.cex=2.8, direction = "outside", major.at=c(0,1,2.1),#seq(from=0,to=floor(2),by=1), 
                            #              minor.ticks=1, labels.away.percentage = 0.15, lwd=2,labels.niceFacing=FALSE)
                             circos.axis(labels.cex=6, direction = "outside", major.at=c(0,1,2.1),#seq(from=0,to=floor(2),by=1), 
                                         minor.ticks=1, lwd=2,labels.niceFacing=FALSE)
                           }) 
    
    #Add nonsignificant links first
    curr_starts_gray=matrix(0L, nrow = 1, ncol = n_components)
    for(i in 1:n_components) {
      for(j in 1:n_components){
        if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
          if(sign_vals[i,j]==0){
            if(plotBand_cort=='theta'){
              data_temp = abs(netData_ave$netVals.ave[[1]])*10^4
              circos.link(factors[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), factors[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.4)) 
             }else{
              data_temp = abs(netData_ave$netVals.ave[[2]])*10^4
              circos.link(factors[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), factors[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.4)) 
            }
            curr_starts_gray[i] = curr_starts_gray[i]+data_temp[i,j]
            curr_starts_gray[j] = curr_starts_gray[j]+data_temp[j,i]
          }
        }
      }
    }
    
    curr_starts=matrix(0L, nrow = 1, ncol = n_components)
    for(i in 1:n_components) {
      for(j in 1:n_components){
        if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
          if(sign_vals[i,j]>0 & sign_vals[j,i]>0){
            circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('red', alpha.f = 0.6))#color_vals_links[8]) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }else if(sign_vals[i,j]<0 & sign_vals[j,i]<0){
            circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('blue', alpha.f = 0.6)) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }else if(sign_vals[i,j]!=0 & sign_vals[j,i]!=0){
            circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('magenta', alpha.f = 0.6)) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }
        }
      }
    }
    dev.off()
  }
}

