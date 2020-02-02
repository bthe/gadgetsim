##' Calculate yield per recruit of a stock in a Gadget model
##'
##' \code{gadget.ypr} returns the results from a yield per recruit simulation. 
##' The function sets up a new Gadget run based on a given model 
##' and parameter settings where the initial conditions and recruiment 
##' parameters are all set to zero appart from one year class which is 
##' set to a nominal value (10000 individuals). The model is the run with 
##' a range of fishing effort settings. The fishing effort is formulated in
##' relation to linearfleets whose selection is the same as the one in the 
##' fleets in the model. The user can specify what fleets will be used in
##' in the yield per recruit simulation.
##' @title Gadget Yield per Recruit
##' @name gadget.ypr
##' @param params.file Parameter file for the gagdet model
##' @param main.file Main file for the gagdet model
##' @param effort The range of fishing mortality
##' @param begin Start year of the simulation
##' @param end End year of the simulation
##' @param age.range at what age range should the YPR be calculated
##' @param fleets Data frame comtaining the fleet names and ratio in
##' future catches
##' @param ssb.stock name of the spawning stock
##' @param mat.par logistic curve parameters for the maturation
##' @param check.previous check if the analysis have been done before
##' @param save.results should the results be saved?
##' @param gd gadget directory object
##' @return a list containing the yield per recruit by F, estimate of
##' Fmax and F0.1
##' @author Bjarki Thor Elvarsson
##' @export
##' @examples \dontrun{
##' ypr <- gadget.ypr(ssb.stock = 'codmat',
##'                   params.file = 'WGTS/params.final')
##' plot(ypr)
##' }
gadget.ypr <- function(params.file = 'params.in',
                       main.file = 'main',
                       effort = seq(0, 1, by=0.01),
                       begin=1990,end=2020,
                       age.range=NULL,
                       fleets = data.frame(fleet='comm',ratio=1),
                       ssb.stock=NULL,
                       mat.par = NULL,
                       check.previous = FALSE,
                       save.results = TRUE,
                       gd=list(dir='.',rel.dir='YPR')){
  
  ypr <- paste(gd$dir,gd$rel.dir,sep='/')
  
  ## ensure that files exist
  if(!file.exists(params.file)) {
    stop('Parameter file not found')
  }
  
  if(!file.exists(main.file)) {
    stop('Main file not found')
  }
  
  ## model setup
  if(check.previous){
    if(file.exists(sprintf('%s/ypr.Rdata',ypr))){
      load(sprintf('%s/ypr.Rdata',ypr))
      return(res)
    }
  }
  
  ## File I/O
  dir.create(ypr,showWarnings = FALSE, recursive = TRUE)
  dir.create(sprintf('%s/Aggfiles',ypr),
             showWarnings = FALSE, recursive = TRUE)
  main <- read.gadget.main(main.file)
  stocks <- read.gadget.stockfiles(main$stockfiles)
  fleet <- read.gadget.fleet(main$fleetfiles)
  params <- read.gadget.parameters(params.file)
  time <- read.gadget.time(main$timefile)
  area <- read.gadget.area(main$areafile)
  
  ## basic setup
  time$lastyear <-  end
  time$firstyear <- begin
  time$laststep <- length(time$notimesteps)
  time$firststep <- 1
  
  time.grid <- expand.grid(year = time$firstyear:time$lastyear,
                           step = 1:length(time$notimesteps),
                           area = area$areas)
  
  area$temperature <- dplyr::mutate(time.grid,
                                    temperature = 5)
  
  main$areafile <- sprintf('%s/area',ypr)
  write.gadget.area(area,file=sprintf('%s/area',ypr))
  write.unix(sprintf('allareas %s',paste(area$areas,collapse=' ')),
             f=sprintf('%s/Aggfiles/allareas.agg',ypr))
  
  fleet <- plyr::llply(fleet,
                       function(x){
                         tmp <- subset(x,fleet %in% fleets$fleet)
                       })
  fleet$fleet <- dplyr::mutate(fleet$fleet,
                               multiplicative = '1#effort',
                               amount = sprintf('%s/fleet.ypr', ypr),
                               type = 'linearfleet')
  
  fleet.predict <- plyr::ddply(fleets,'fleet',function(x){
    tmp <- dplyr::mutate(time.grid,
                         ratio = x$ratio)
    return(tmp)
  })
  
  
  write.gadget.table(fleet.predict[c('year','step','area','fleet','ratio')],
                     file=sprintf('%s/fleet.ypr',ypr),
                     col.names=FALSE,row.names=FALSE,
                     quote = FALSE)
  
  main$fleetfiles <- sprintf('%s/fleet', ypr)
  write.gadget.fleet(fleet,file=sprintf('%s/fleet', ypr))
  
  write.gadget.time(time,file=sprintf('%s/time.ypr',ypr))
  main$timefile <- sprintf('%s/time.ypr',ypr)
  
  ## basic printfile
  
  print.txt <-
    paste('[component]',
          'type\tpredatorpreyprinter',
          sprintf('predatornames\t%s',
                  paste(fleets$fleet,collapse=' ')),
          'preynames\t%1$s',
          'areaaggfile\t%2$s/Aggfiles/allareas.agg',
          'ageaggfile\t%2$s/Aggfiles/%1$s.allages.agg',
          'lenaggfile\t%2$s/Aggfiles/%1$s.alllen.agg',
          'printfile\t%2$s/out/%1$s.prey',
          'yearsandsteps\tall all',
          sep = '\n')
  
  print.ssb <- NULL
  if(!is.null(ssb.stock)){
    if(sum(ssb.stock %in% names(stocks)) == length(ssb.stock)){
      print.ssb <-
        paste('[component]',
              'type\tstockprinter',
              'stocknames\t%1$s',
              'areaaggfile\t%2$s/Aggfiles/allareas.agg',
              'ageaggfile\t%2$s/Aggfiles/%1$s.allages.agg',
              'lenaggfile\t%2$s/Aggfiles/%1$s.alllen.agg',
              'printfile\t%2$s/out/%1$s.ssb',
              'yearsandsteps\tall 1',
              sep = '\n')
      
    } else {
      warning(sprintf('warning in gadget.ypr, SSB stock %s not found',
                      ssb.stock))
      
    }
  }
  print.mat <- NULL
  if(!is.null(mat.par)){
    print.mat <-
      paste('[component]',
            'type\tstockprinter',
            'stocknames\t%1$s',
            'areaaggfile\t%2$s/Aggfiles/allareas.agg',
            'ageaggfile\t%2$s/Aggfiles/%1$s.allages.agg',
            'lenaggfile\t%2$s/Aggfiles/%1$s.len.agg',
            'printfile\t%2$s/out/%1$s.mat',
            'yearsandsteps\tall 1',
            sep = '\n')
  }
  
  
  
  printfile <- paste(sprintf(print.txt,unique(fleet$prey$stock), ypr),
                     collapse='\n;\n')
  if(!is.null(ssb.stock)){
    printfile <-
      paste(printfile,
            paste(sprintf(print.ssb,ssb.stock,ypr),
                  collapse='\n;\n'),
            sep='\n;\n')
  }
  
  if(!is.null(mat.par)){
    printfile <-
      paste(printfile,
            paste(sprintf(print.mat,unique(fleet$prey$stock),ypr),
                  collapse='\n;\n'),
            sep='\n;\n')
  }
  
  
  dir.create(sprintf('%s/out',ypr),showWarnings = FALSE, recursive = TRUE)
  main$printfiles <- sprintf('%s/printfile.ypr',ypr)
  write.unix(printfile,f=sprintf('%s/printfile.ypr',ypr))
  
  
  ## remove recruitment and initialdata from the stockfiles
  
  plyr::l_ply(stocks,function(x){
    x@initialdata[,3] <- 0 ## nothing in the beginning
    if(x@doesrenew==1){
      x@renewal.data <- 
        time.grid %>% 
        dplyr::filter(.data$step == 1) %>% 
        dplyr::bind_cols(dplyr::slice(x@renewal.data %>% 
                                        dplyr::select(-c(.data$year,.data$step,.data$area)),rep(1,nrow(.)))) %>% 
        dplyr::mutate(number = ifelse(.data$year==begin,1,0))
      x@doesspawn <- 0
    }
    gadget_dir_write(gd,x)
  })
  
  main$stockfiles <- sprintf('%s/%s',ypr,
                             plyr::laply(stocks,function(x) x@stockname))
  
  main$likelihoodfiles <- ';'
  
  write.gadget.main(main,file=sprintf('%s/main.ypr',ypr))
  
  ## model parameters
  if(sum(names(params) %in% c('switch','value',
                              'lower','upper','optimise'))==5){
    tmp <- as.data.frame(t(params$value))
    names(tmp) <- params$switch
    params <- tmp
  }
  
  
  
  params.aug <- plyr::ldply(effort,
                            function(x){
                              tmp <- params
                              tmp$effort <- x
                              return(tmp)
                            })
  
  write.gadget.parameters(structure(params.aug,file_format = 'wide'),file=sprintf('%s/params.ypr',ypr))
  
  callGadget(s=1,i=sprintf('%s/params.ypr',ypr),main=sprintf('%s/main.ypr',ypr))
  
  ## read output
  
  out <- plyr::ddply(data.frame(stock = unique(fleet$prey$stock),tmp=1),
                     'stock',
                     function(x){
                       stock.prey <- utils::read.table(file = sprintf("%1$s/out/%2$s.prey",
                                                                      ypr,x$stock),
                                                       comment.char = ';')
                       
                       names(stock.prey) <-
                         c('year', 'step','area','age','length','number.consumed',
                           'biomass.consumed','fishing.mortality')
                       
                       stock.prey$trial <-
                         rep(1:c(nrow(stock.prey)/(
                           length(unique(stock.prey$area))*
                             length(unique(stock.prey$step))*
                             length(unique(stock.prey$year)))),
                           each=length(unique(stock.prey$area))*
                             length(unique(stock.prey$year))*
                             length(unique(stock.prey$step)))
                       stock.prey <- merge(stock.prey,
                                           data.frame(trial=1:length(effort),
                                                      effort=effort),
                                           all.x=TRUE)
                       stock.prey$age <- stock.prey$year - begin +
                         getMinage(stocks[[x$stock]])
                       ## clean up
                       file.remove(sprintf('%s/out/%s.prey',ypr,x$stock))
                       
                       return(stock.prey)
                     })
  if(!is.null(ssb.stock)){
    ssb.out <- plyr::ldply(ssb.stock,function(x){
      ssb.out <- utils::read.table(file = sprintf("%1$s/out/%2$s.ssb",
                                                  ypr,x), comment.char = ';')
      file.remove(sprintf('%s/out/%s.ssb',ypr,x))
      names(ssb.out) <-
        c('year', 'step','area','age','length','number',
          'biomass')
      ssb.out$trial <-
        rep(1:c(nrow(ssb.out)/(
          length(unique(ssb.out$area))*
            length(unique(ssb.out$step))*
            length(unique(ssb.out$year)))),
          each=length(unique(ssb.out$area))*
            length(unique(ssb.out$year))*
            length(unique(ssb.out$step)))
      ssb.out <- merge(ssb.out,
                       data.frame(trial=1:length(effort),
                                  effort=effort),
                       all.x=TRUE)
      ssb.out %>% 
        dplyr::group_by(.data$effort) %>% 
        dplyr::summarise(ssb = max(.data$number*.data$biomass)) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(ssb.ratio=.data$ssb/max(.data$ssb))
      
      
    })
    
  } else {
    ssb.out <- NULL
  }
  
  if(!is.null(mat.par)){
    mat.out <- plyr::ldply(unique(fleet$prey$stock),function(x){
      mat.out <- 
        utils::read.table(file = sprintf("%1$s/out/%2$s.mat",
                                         ypr,x), comment.char = ';')
      file.remove(sprintf('%s/out/%s.mat',ypr,x))
      names(mat.out) <-
        c('year', 'step','area','age','length','number',
          'biomass')
      mat.out %>% 
        dplyr::mutate(length = as.numeric(gsub('len','',.data$length)),
                      trial = cut(1:length(.data$year),c(0,which(diff(.data$year)<0),1e9),labels = FALSE)) %>% 
        dplyr::left_join(data.frame(trial=1:length(effort),
                                    effort=effort)) %>% 
        dplyr::group_by(.data$effort,.data$year) %>% 
        dplyr::summarise(ssb=sum(.data$number*.data$biomass*logit(mat.par[1],mat.par[2],.data$length))) %>% 
        dplyr::group_by(.data$effort) %>% 
        dplyr::mutate(ssb.ratio = .data$ssb/max(.data$ssb))
      
    })
    
  } else {
    mat.out <- NULL
  }
  
  
  
  if(!is.null(age.range)){
    out <- subset(out,out$age >= min(age.range) & out$age <= max(age.range))
  }
  res <- plyr::ddply(out,'effort',
                     function(x) c(num=sum(x$number.consumed)/1e6,
                                   bio=sum(x$biomass.consumed)/1e6))
  secant <- diff(res$bio)/diff(res$effort)
  f0.1 <- res$effort[min(which(secant<0.1*secant[1]))]
  fmax <- min(res$effort[which(res$bio==max(res$bio,na.rm=TRUE))])
  res <- list(params=params,out=out,ypr=res,fmax=fmax,
              f0.1=f0.1,ssb=ssb.out,mat.out=mat.out)
  
  
  class(res) <- c('gadget.ypr',class(res))
  if(save.results){
    save(res, file = sprintf('%s/ypr.Rdata',ypr))
  }
  
  return(res)
}

##' @rdname gadget.ypr
##' @name plot.gadget.ypr
##' @param x Output from \code{gadget.ypr}
##' @param ... Unused
##' @export
plot.gadget.ypr <- function(x, ...) {
  ypr <- x
  if(!is.null(ypr$ssb)){
    tmp <- merge(dplyr::mutate(ypr$ypr,bio=.data$bio/max(.data$bio)),
                 ypr$ssb)
  } else {
    tmp <- dplyr::mutate(ypr$ypr,bio=.data$bio/max(.data$bio))
  }
  plo <- 
    ggplot2::ggplot(tmp,ggplot2::aes(.data$effort,.data$bio)) +
    ggplot2::geom_line() +
    ggplot2::geom_segment(ggplot2::aes(x = .data$effort,xend=.data$effort,y=-Inf,yend=.data$bio),
                          data=subset(tmp, tmp$effort == ypr$fmax)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$effort,xend=.data$effort,y=-Inf,yend=.data$bio),
                          data=subset(tmp, tmp$effort == ypr$f0.1)) +
    ggplot2::geom_text(data=subset(tmp,tmp$effort == ypr$fmax),
                       ggplot2::aes(label = sprintf('Fmax = %s',.data$effort),
                                    x = .data$effort+0.04,y=0.2,angle=90)) +
    ggplot2::geom_text(data=subset(tmp,tmp$effort == ypr$f0.1),
                       ggplot2::aes(label = sprintf('F0.1 = %s',.data$effort),
                                    x = .data$effort+0.04,y=0.2,angle=90)) +
    ggplot2::labs(x='Fishing mortality', y='Yield per recruit') +
    ggplot2::theme(legend.position='none',plot.margin = ggplot2::unit(c(0,0,0,0),'cm'))
  if(!is.null(ypr$ssb)){
    plo <- plo + ggplot2::geom_line(ggplot2::aes(.data$effort,.data$ssb.ratio),
                                    lty=2,col='gray')
  }
  return(plo)
}
