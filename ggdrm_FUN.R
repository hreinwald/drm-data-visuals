# Author:  Hannes Reinwald
# Contact: hannes.reinwald@bayer.com


### ggdrm() ### ----------------------------------------------------------------
#' ggdrm()    A function to plot and visualize DR-models created with drm() from the drc-package.
#'            Further the function can display EC values and computes model averaged EC values using drc's maED() function
#'            based on the fct.ls provided as input. ggdrm() also allows the export of these tables as a data frame object.
#' mod:       A drc::drm() model object
#' fct.ls:    List of dr-model functions which should be plotted for which EC values should be computed.
#' CIlevel:     Conf.Intervals to display. Default = 0.95 (= 95%)
#' respLev:   Numeric vector describing the relative effect concentrations levels which should be returned. Default c(20,50)
#' conc.unit: Character string describing the concentration unit. If not provided extracted from the dose argument from mod 
#' vcovFun:   Function used to Calculate Variance-Covariance Matrix for a Fitted Model Object. Default: vcov()
#' title:     drm-Plot title
#' subt:      drm-Plot subtitle 
#' ggtheme:   ggplot theme used for drm-plot visuals

#' F0: Numeric value by which to divide the lowest tested conc. which will replace any 0 values as they create problems with the log scale!!!
#' If not specified F0 value is computed as: F0 = ( sort(unique(df$conc))[3] / sort(unique(df$conc))[2] )*10
#' Which is the fold change between lowest test concentration (C1) and 2nd lowest (C2) times 10:
#' Default F0 = (C2/C1) * 10

## Required Packages ##
require(dplyr)
require(ggplot2)
require(ggpubr)
require(scales)
require(drc)

## ggdrm() functions ##

# sub-functions
getMspec = function(mod){ mod$fct$names %>% paste0(.,collapse = "-") %>% paste0(mod$fct$name,":",.) } # get model specifications
buildECtable = function(mod, unit = NULL, respLev = c(10,20,50), CIlevel = .9, covFUN=vcov){
  EC = ED(mod, respLev = respLev, "delta", level = CIlevel, display = F, vcov. = covFUN) %>% as.data.frame
  colnames(EC)[2:4] = colnames(EC)[3:4] %>% paste0(CIlevel*100,"CI.",.) %>% c("SE",.)
  EC = round(EC,2)
  EC$unit = if(is.null(unit)){ colnames(mod$data)[1] }else{ unit }
  EC$mtype = getMspec(mod)
  EC$ECXX = rownames(EC) %>% sub("^e:1:","EC",.)
  EC$AIC = AIC(mod) %>% round(.,1)
  
  # Append model summary stats
  df = summary(mod)[["rseMat"]] %>% round(.,2) %>% as.data.frame
  if(mod$type == "binomial"){ df = df[,2] }
  # if continous data, extract resVar
  if(mod$type == "continuous"){
    tmp = summary(mod)[["resVar"]] %>% round(.,2) %>% as.data.frame
    colnames(tmp) = "resVar"
    df = cbind(df,tmp)
  }
  EC = cbind(EC,df)
  
  # Add model fit (LOF) pvalue as well 
  lof = try( drc::modelFit(mod, method = if(mod$type == "continuous"){"gof"}else{"cum"})[2,5] %>% round(.,3) ,silent = T)
  if(inherits(lof, "try-error")){ lof = NA }
  EC$LOF.pval = lof
  return( EC[,c(7,1:6,8:ncol(EC))] )
}
buildMECtable = function(mod, fct.ls, unit = NULL, respLev = c(10,20,50), CIlevel = .9, verbose=F){
  EC = drc::maED(mod, fct.ls, respLev,"buckland",level = CIlevel, display = verbose, na.rm = T) %>% as.data.frame(.)
  colnames(EC)[2:4] = colnames(EC)[3:4] %>% paste0(CIlevel*100,"CI.",.) %>% c("SE",.)
  EC = round(EC,2)
  EC$unit = if(is.null(unit)){ colnames(mod$data)[1] }else{ unit }
  EC$mtype = c("maED")
  EC$ECXX = rownames(EC) %>% sub("^e:1:","EC",.)
  return( EC[,c(7,1:6)] )
}

# main plotting function
ggdrm = function(mod, fct.ls = NULL, respLev = c(10,20,50), CIlevel = .9, conc.unit = NULL, 
                 showEC = F, avEC = F, return.list = F, vcovFun=vcov, F0 = NULL,
                 title = NULL, subt = NULL, y.lim = NULL,
                 ggtheme = theme_light()+theme(legend.position="bottom")+rotate_x_text(30) ){
  
  # Extract df for plotting from drm() output object -------------------------------------------------
  df = mod$data
  D = colnames(df)[1] # name of Dose
  R = colnames(df)[2] # name of Response
  colnames(df)[1:2] = c("conc","resp")
  
  stopifnot(length(unique(df$conc)) > 2) # make sure that at least 3 different conc are provided!
  maxc = max(df$conc)
  #minc = sort(unique(df$conc))[2] #.[. != 0]
  minc = df$conc[df$conc != 0] %>% sort %>% .[1] # <-- this is better in case dataset does not contain a 0 test conc!
  minc2 = df$conc[df$conc != 0] %>% unique %>% sort %>% .[2]
  
  # plotting concentrations on a log scale will result in problems!!! We help by shifting the 0 concentration a bit up
  df$conc0 <- df$conc
  if(is.null(F0)){ F0 = (minc2/minc) * 10 } # specify F0 value by which to divide minc
  df$conc0[df$conc0 == 0] <- minc/F0
  
  # Predicting values with confidence intervals for all generated dose levels ------------------------
  
  # Generating new dose levels on the log scale as support for the line.
  f = (minc2/minc)*.5 # Extend the model plot to 1/2 of the dilution step space at the upper dose range
  if(f < 1){f = 1.1}
  #f = 2
  newdata = expand.grid(conc = exp(seq(log((minc/F0)/1.1),log(maxc*f),length=250)) ) # log transform for more even distributin of values!
  
  # function to compute model predictions with upper and lower CI
  my_PM = function(mod,tmpDf){
    pm = predict(mod, newdata=tmpDf, interval="confidence", level = CIlevel) # see: ?predict.drc for details
    tmpDf$p = pm[,1]
    tmpDf$pmin = pm[,2]
    tmpDf$pmax = pm[,3]
    return(tmpDf) }
  vis1 = my_PM(mod, newdata) # data points for model visualization 1
  
  # append model specifications
  vis1$Model = getMspec(mod)
  
  #' if !is.null(fct.ls) : 
  #' Create an individual vis object for each model from the fct.ls. 
  #' Then combine all in a single visual object. 
  #' Remember to use try(..., silent=T) when running update() on the initial model, to avoid breaking the function. 
  #' example code from: https://stackoverflow.com/questions/21192002/how-to-combine-2-plots-ggplot-into-one-plot 
  if(!is.null(fct.ls)){
    # update model - make sure to use try()!
    m.ls = lapply(fct.ls, function(fct){
      mod2 = try( update( mod, fct = fct, data = mod[["origData"]], control=drc::drmc(errorm = T)), silent = T)
      if(!inherits(mod2, "try-error")){
        # if model fitting worked predict values based on the newdata grid computed earlier & append the model name as info to vis2 
        vis2 = my_PM(mod2, newdata)
        vis2$Model = getMspec(mod2)
        return( list(model = mod2, vis = vis2) )
      } else { warning("Unable to find a suitable model fit for:\t",fct$name) }
    })
    names(m.ls) = lapply(m.ls, function(x){getMspec(x[["model"]])}) %>% unlist()
    
    # combine vis1 and vis2 to new df
    vis2 = lapply(m.ls, function(x) x[["vis"]]) %>% rlist::list.rbind(.) %>% na.omit(.)
    row.names(vis2) = NULL
    if(nrow(vis2) > 0){ vis1 = rbind(vis1, vis2)}
  }
  
  # Plotting the fitted DR-model(s) -----------------------------------
  br = unique(df$conc0) %>% sort(.) %>% .[-1] 
  gg.drg = ggplot(df, aes(x = conc0, y = resp)) +
    # Extend plot with model lines and CI
    geom_ribbon(data=vis1, aes(x=conc, y=p, ymin=pmin, ymax=pmax,
                               group=Model,col=Model, fill=Model), alpha=.25, linetype = "blank") +
    geom_line(data=vis1, aes(x=conc, y=p,  group=Model,col=Model), alpha=.75, linetype = "solid", linewidth=.8) +
    # add data points and rest
    geom_point(alpha=0.7) + #geom_point(position = position_jitter(width = 0, height = .02)) + 
    coord_trans(x="log") + xlab(D) + ylab(R) +
    labs(title = title, subtitle = paste("Data type:",mod$type,subt)) +
    # Custom breaks 
    scale_x_continuous(limits = c(NA, max(br)*f), expand = c(0,0.5), breaks = br, 
                       labels = scales::scientific(br, digits = 3)) + ggtheme
  
  # set scales for y axis
  if(mod$type == "binomial" & is.null(y.lim)){ y.lim = c(0,1) } # If not specified otherwise, for binomial data y.lim = c(0,1)
  if(!is.null(y.lim)){ gg.drg = gg.drg + ylim(y.lim) }
  
  ## OUTPUT LIST OBJECT ## 
  OUT.ls = list(gg.drg = gg.drg)
  
  # If showEC == T, add stats info to the plot -------------------------------------------------
  if(showEC|avEC){
    
    if( is.null(conc.unit) ){ conc.unit = D } # by default conc.unit will be replaced by the dose (D) argument. If you don't wish this set conc.unit = NA
    tb1 = buildECtable(mod, conc.unit, respLev, CIlevel, vcovFun)
    
    if(!is.null(fct.ls)){
      if(!is.null(m.ls)){ # if a m.ls was created generated EC/Mstats tables for other models as well. 
        tb2 = lapply(m.ls, function(x){ buildECtable( x[["model"]], conc.unit, respLev, CIlevel, vcovFun) }) %>% 
          rlist::list.rbind(.)
        if( nrow(tb2)>0 ){ tb1 = rbind(tb1,tb2) }
      } 
    }
    # Store in output list
    drop.rownames = function(x){
      row.names(x) = NULL
      return(x)}
    OUT.ls[["EC.all"]] = tb1[order(tb1$ECXX),] %>% drop.rownames
    OUT.ls[["gg.drt"]] = ggpubr::ggtexttable(tb1, theme = ttheme("mGreen"), rows = NULL)
  }
  
  # If wanted - compute and display the maED values ---------------------------------------------
  if(!is.null(fct.ls)){
    if(avEC & !is.null(m.ls)){
      
      if( is.null(conc.unit) ){ conc.unit = D }
      maEC = try( buildMECtable(mod, fct.ls, conc.unit, respLev, CIlevel=CIlevel), silent = F)
      if(!inherits(maEC, "try-error")){
        
        OUT.ls[["gg.maEC"]] = ggpubr::ggtexttable(maEC, theme = ttheme("mOrange"), rows = NULL)
        OUT.ls[["maEC"]] = maEC
        
      } else { warning("Unable to compute averaged model fit for provided dataset or model list. Check your input.") }
    } 
  }
  
  # Organize output of function -----------------------------------------------------------------
  if(return.list){ return(OUT.ls) }else{
    gg = grep("gg[.]", names(OUT.ls) ) %>% OUT.ls[.]
    # Print final plots
    if(length(gg) == 1){ g = gg[[1]] }
    if(length(gg) == 2){ g = ggpubr::ggarrange(plotlist=gg, ncol = 1, nrow = 2, heights = c(1,.6)) } 
    if(length(gg) == 3){ 
      if(!showEC){ 
        g = ggpubr::ggarrange(plotlist=gg[-2], ncol = 1, nrow = 2, heights = c(1,.6)) 
      } else { 
        g = ggpubr::ggarrange(plotlist=gg[c(1,3,2)], ncol = 1, nrow = 3, heights = c(1,.6,.6))
      }
    }
    return(g)
  } 
}
#######################


## See example code below ## ---------------------------------------------------------------------------
run = F
if(run){
  ### SOME DATA SETS to play around with ### ----------------------------------------------------------
  
  # Quantal / binomial data
  mod = drm(number/total ~ dose, data = earthworms, type = "binomial", weights = total, fct = LL.4()) # <- note the weights argument I used here!
  mod = drm(numdead/total ~ dose, data = H.virescens, type = "binomial", weights = total, fct = LL.4())
  
  # Continous data
  mod = drm(weight ~ conc, data = O.mykiss, fct = LL.4())
  mod = drm(rootl ~ conc, data = ryegrass, fct = LL.4())
  
  ## Choose most suitable models ## ------------------------------------------------------------------
  
  # List of functions you want to visualize in your drm-plot and from which to compute maED values
  fct.ls = list(LL.3(), W1.3(), W2.3())
  
  # simple way of comparing those models 
  modSel = mselect(mod, fct.ls, nested = T) # Models with the lowest IC value are usually best!
  modSel # inspect your output 
  
  # You could consider only using the top 3 models for computing averaged EC values. 
  topM = row.names(modSel)[1:3]
  names(fct.ls) = lapply(fct.ls, function(x) x$name) %>% unlist # name your input fct.ls after model names to select via topM
  fct.ls2 = fct.ls[ names(fct.ls) %in% topM ]
  
  
  ## Plot topM models, respective EC values and average EC values ## ---------------------------------
  # Below you find different examples how to use the ggdrm() plotting function. 
  ggdrm(mod)          # <- simplest way to display your drm-plot
  ggdrm(mod, fct.ls)  # <- with your list of model functions to compare to.
  ggdrm(mod, fct.ls2) + theme_light() # <- our best fit model collection
  
  # If you don't like the look you can make your own using ggplot's theme options
  my_theme = theme_classic2() + theme(legend.position="bottom") # <- specify your own ggplot theme for the drm-plot
  ggdrm(mod, fct.ls2, ggtheme = my_theme)
  
  ggdrm(mod, fct.ls2, showEC = T) # <- To display your EC values of interest simply set showEC = T
  ggdrm(mod, fct.ls2, avEC = T)   # <- To display average EC values from input mod and given fct.ls set avEC = T
  
  # To display all tables with the drm-plot run the following:
  ggdrm(mod, fct.ls2, avEC = T, showEC = T, title = "Test", subt = "subtitle", ggtheme =  my_theme)
  
  # If you wish to change the obtained EC values or CI you can do so by simply setting those parameters as: 
  ggdrm(mod, fct.ls2, avEC = T, showEC = T, title = "Test", subt = "subtitle", CIlevel = .8,
        respLev = c(10,25,35) )  # <- Ups! At some point it get's tricky to display all values in a single plot
  # But you can display each object of the plot indivudally by using return.list = T
  
  ## To export every single component of the ggdrmplot in a list object set return.ist = T ---------------
  g.drm = ggdrm(mod, fct.ls2, avEC = T, showEC = T, respLev = c(10,20,30,40,50), return.list = T)
  
  # Export tables from plots as df
  g.drm$maEC
  g.drm$EC.all
  
  # display only the tables
  g.drm$gg.drt
  g.drm$gg.maEC
  
  # display only the ggplot
  g = g.drm$gg.drg 
  g # <- this is a ggplot object, hence we can simply extend layers of the plot with ggplot expressions
  g + theme_classic2() + labs(title = "Fancy plot! :)", subtitle = "So cool!") 
}
############ END ###########