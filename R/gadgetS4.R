

##' @rdname gadgetFitHelpers
##' @description \code{get.gadget.recruitment} gets gadget recruitment
##' @param stocks gadget.stock
##' @param params gadget.params
##' @param collapse Logical, group output by stock/year/area
##' @return recruitment by year
get.gadget.recruitment <- function(stocks,params,collapse=TRUE){
  plyr::ldply(stocks, function(x){
    if(x@doesrenew == 1){
      tmp <- 
        x@renewal.data %>% 
        dplyr::mutate(stock = x@stockname,
                      recruitment = 1e4*unlist(eval.gadget.formula(.data$number,params))) %>% 
        dplyr::select(.data$stock,.data$year,.data$step,.data$area,.data$recruitment) %>% 
        stats::na.omit() 
      if(collapse){
        tmp %>% 
          # area added to accomodate merge with res.by.year in gadget.fit() - pfrater
          dplyr::group_by(.data$stock,.data$year, .data$area) %>% 
          dplyr::summarise(recruitment=sum(.data$recruitment)) %>% 
          as.data.frame()
      } else{
        tmp
      }
    } else {
      data.frame(stock = x@stockname,year=NA,area=NA,recruitment=NA)
    }
  })
}


##' @rdname gadgetFileIO
##' @description \code{read.gadget.stockfiles} reads in Gadget stockfiles DEPRECATED
##' @param stock.files location of stock file
##' @return list of gadget-stock objects
read.gadget.stockfiles <- function(stock.files){
  tmp.func <- function(sf){
    stock <- strip.comments(sf)
    
    growth.loc <- grep('doesgrow', stock, ignore.case = TRUE)
    mort.loc <- grep('naturalmortality', stock, ignore.case = TRUE)
    init.loc <- grep('initialconditions', stock, ignore.case = TRUE)
    pred.loc <- init.loc - 1
    if(tolower(stock[init.loc+1])=='numbers'){
      init.loc <- init.loc + 1
    }
    #    initfile.loc <- #grep('normalcondfile',stock, ignore.case =TRUE)
    #      c(grep('normalcondfile',stock, ignore.case =TRUE),
    #        grep('normalparamfile',stock, ignore.case =TRUE),
    #        grep('numberfile',stock, ignore.case =TRUE))
    eat.loc <- grep('doeseat', stock, ignore.case = TRUE)
    migrate.loc <- grep('doesmigrate', stock, ignore.case = TRUE)
    initfile.loc <- migrate.loc -1
    mature.loc <- grep('doesmature', stock, ignore.case = TRUE)
    move.loc <- grep('doesmove', stock, ignore.case = TRUE)
    renew.loc <- grep('doesrenew', stock, ignore.case = TRUE)
    spawn.loc <- grep('doesspawn', stock, ignore.case = TRUE)
    stray.loc <- grep('doesstray', stock, ignore.case = TRUE)
    
    growth.info <- function(tmp){
      if(length(tmp)==1)
        tmp <- methods::new('gadget-growth')
      else {
        names.tmp <- sapply(tmp,function(x) x[1])
        tmp <- plyr::llply(tmp,function(x) paste(x[-1],collapse=' '))
        names(tmp) <- names.tmp
        
        if(is.null(tmp$growthparameters))
          tmp$growthparameters <- vector()
        if(is.null(tmp$growthfunction))
          tmp$growthfunction <- vector()
        if(is.null(tmp$wgrowthfunction))
          tmp$wgrowthfunction <- vector()
        if(is.null(tmp$lgrowthfunction))
          tmp$lgrowthfunction <- vector()
        if(is.null(tmp$yeareffect))
          tmp$yeareffect <- vector()
        if(is.null(tmp$stepeffect))
          tmp$stepeffect <- vector()
        if(is.null(tmp$areaeffect))
          tmp$areaeffect <- vector()
        if(tolower(tmp$growthfunction) == 'lengthvbsimple'){
          tmp <- methods::new('gadget-growth',
                              growthfunction = tmp$growthfunction,
                              ## growthfunction parameters
                              growthparameters = tmp$growthparameters,
                              #                   wgrowthparameters = tmp$wgrowthparameters,
                              #                   lgrowthparameters = tmp$lgrowthparameters,
                              #                   yeareffect = tmp$yeareffect,
                              #                   stepeffect = tmp$stepeffect,
                              #                   areaeffect = tmp$areaeffect,
                              ## growth implementation
                              beta = tmp$beta,
                              maxlengthgroupgrowth = tmp$maxlengthgroupgrowth)
        } else if(tolower(tmp$growthfunction) == 'lengthpower'){
          tmp <- methods::new('gadget-growth',
                              growthfunction = tmp$growthfunction,
                              ## growthfunction parameters
                              growthparameters = tmp$growthparameters,
                              weightgrowthfile = tmp$weightgrowthfile,
                              #                   lgrowthparameters = tmp$lgrowthparameters,
                              #                   yeareffect = tmp$yeareffect,
                              #                   stepeffect = tmp$stepeffect,
                              #                   areaeffect = tmp$areaeffect,
                              ## growth implementation
                              beta = tmp$beta,
                              maxlengthgroupgrowth = tmp$maxlengthgroupgrowth)
        } else {
          tmp <- methods::new('gadget-growth',
                              growthfunction = tmp$growthfunction,
                              ## growthfunction parameters
                              growthparameters = tmp$growthparameters,
                              wgrowthparameters = tmp$wgrowthparameters,
                              lgrowthparameters = tmp$lgrowthparameters,
                              ##                   yeareffect = tmp$yeareffect,
                              ##                   stepeffect = tmp$stepeffect,
                              ##                   areaeffect = tmp$areaeffect,
                              ## growth implementation
                              beta = tmp$beta,
                              maxlengthgroupgrowth = tmp$maxlengthgroupgrowth)
        }
      }
      return(tmp)
    }
    
    prey.info <- function(tmp){
      if(length(tmp)==1){
        tmp <- methods::new('gadget-prey')
      } else {
        tmp <- methods::new('gadget-prey',
                            name = stock[[1]][2],
                            preylengths = utils::read.table(tmp[[2]][2],comment.char=';'),
                            energycontent = ifelse(length(tmp)==3,as.numeric(tmp[[3]][2]),
                                                   1))
      }
      return(tmp)
    }
    
    pred.info <- function(tmp){
      if(length(tmp)==1){
        tmp <- methods::new('gadget-predator')
      } else {
        suit.loc <- grep('suitability',tmp)
        pref.loc <- grep('preference',tmp)
        maxcon.loc <- grep('maxconsumption',tmp)
        half.loc <- grep('halffeedingvalue',tmp)
        suit <- plyr::ldply((suit.loc+1):(pref.loc-1),
                            function(x){
                              c(stock = tmp[[x]][1],
                                suitability = paste(tmp[[x]][-1],collapse=' '))
                            })
        pref <- plyr::ldply((pref.loc+1):(maxcon.loc-1),
                            function(x){
                              c(stock = tmp[[x]][1],
                                preference = paste(tmp[[x]][-1],collapse=' '))
                            })
        tmp <- methods::new('gadget-predator',
                            suitability = suit,
                            preference = pref,
                            maxconsumption = paste(tmp[[maxcon.loc]][-1], collapse = ' '),
                            halffeedingvalue = paste(tmp[[half.loc]][2], collapse = ' '))
      }
      return(tmp)
    }
    
    initialdata <- read.gadget.table(stock[[initfile.loc]][2])
    if(length(names(initialdata)) == 7){
      names(initialdata) <- c('age','area','age.factor','area.factor',
                              'mean','stddev','relcond')
    } else if (length(names(initialdata)) == 8){
      names(initialdata) <- c('age','area','age.factor','area.factor',
                              'mean','stddev','alpha','beta')
    } else if (length(names(initialdata)) == 5){
      names(initialdata) <- c('area','age','length','number','weight')
    }
    renewal.data <-
      tryCatch(read.gadget.table(stock[[renew.loc+3]][2]),
               error=function(x){
                 tryCatch(read.gadget.table(stock[[renew.loc+4]][2]),
                          error=function(x) data.frame(text='No renewal data'))
               })
    if(length(names(renewal.data)) == 8){
      names(renewal.data) <- c('year','step','area','age','number',
                               'mean','stddev','relcond')
    } else if (length(names(renewal.data)) == 9){
      names(renewal.data) <- c('year','step','area','age','number',
                               'mean','stddev','alpha','beta')
    } else if (length(names(renewal.data)) == 6){
      names(renewal.data) <- c('year','area','age','length','number','weight')
    }
    doesmature <- as.numeric(stock[[mature.loc]][2])
    if(doesmature == 1){
      maturity.function <- tolower(stock[[mature.loc+1]][2])
      maturity.file <- strip.comments(stock[[mature.loc+2]][2])
      maturestocksandratios <- (maturity.file[[1]][-1])
      coefficients <- (maturity.file[[2]][-1])
      maturitylengths <- ''
      maturitysteps <- ''
      if(tolower(maturity.function) == 'fixedlength'){
        maturitylengths <- (maturity.file[[3]][-1])
      } else if(tolower(maturity.function) != 'continuous'){
        maturitysteps <- (maturity.file[[3]][-1])
      } 
      
    } else {
      maturity.function <- ''
      maturestocksandratios <- ''
      coefficients <- ''
      maturitylengths <- ''
      maturitysteps <- ''
    }
    doesmove <- as.numeric(stock[[move.loc]][2])
    if(doesmove == 1){
      transitionstocksandratios <- (stock[[move.loc+1]][-1])
      transitionstep <- as.numeric(stock[[move.loc+2]][-1])
    } else {
      transitionstocksandratios <- ''
      transitionstep <- 0
    }
    
    doesmigrate <- as.numeric(stock[[migrate.loc]][2])
    if(doesmigrate == 1){
      yearstep <- utils::read.table(stock[[migrate.loc+1]][-1],
                                    comment.char = ';',
                                    stringsAsFactors=FALSE)
      tmp <- strip.comments(stock[[migrate.loc+2]][-1])
      mat.loc <- grep('[migrationmatrix]',tmp,fixed=TRUE)
      tyler <- diff(mat.loc)[1]
      if(length(tyler) == 0){
        tyler <- length(tmp)-1
      }
      migrationratio <-
        plyr::llply(mat.loc,  ## all puns intended;)
                    function(x){
                      plyr::laply((x+2):(x+tyler-1),
                                  function(y){merge.formula(tmp[[y]])})
                    })
      names(migrationratio) <- plyr::laply(mat.loc,function(x){ tmp[[x+1]][2]})
    } else {
      migrationratio <- list()
      #      transitionstocksandratios <- ''
      yearstep <- data.frame()
    }
    
    doesspawn <- as.numeric(stock[[spawn.loc]][2])
    if(doesspawn == 1){
      tmp <- strip.comments(stock[[spawn.loc+1]][2])
      ssar <- tmp[[5]][-1]
      nssar <- 2*(1:(length(ssar)/2))-1
      ssar <- data.frame(stock=ssar[nssar],ratio=ssar[1+nssar])
      stockparameters <- as.data.frame(t(merge.formula(tmp[[10]][-1])))
      names(stockparameters) <- c('mean','std.dev','alpha','beta')
      spawning <-
        methods::new('gadget-spawning',
                     spawnsteps = as.numeric(tmp[[1]][-1]),
                     spawnareas = as.numeric(tmp[[2]][-1]),
                     firstspawnyear = as.numeric(tmp[[3]][-1]),
                     lastspawnyear = as.numeric(tmp[[4]][-1]),
                     spawnstocksandratio = ssar,
                     proportionfunction = merge.formula(tmp[[6]][-1]),
                     mortalityfunction = merge.formula(tmp[[7]][-1]),
                     weightlossfunction = merge.formula(tmp[[8]][-1]),
                     recruitment = merge.formula(tmp[[9]][-1]),
                     stockparameters = stockparameters)
      
      
    } else {
      spawning <- methods::new('gadget-spawning')
    }
    
    
    
    
    st <-
      methods::new('gadget-stock',
                   stockname = stock[[1]][2],
                   livesonareas = as.numeric(stock[[2]][-1]),
                   minage = as.numeric(stock[[3]][2]),
                   maxage = as.numeric(stock[[4]][2]),
                   minlength = as.numeric(stock[[5]][2]),
                   maxlength = as.numeric(stock[[6]][2]),
                   dl = as.numeric(stock[[7]][2]),
                   refweight = utils::read.table(stock[[8]][2],comment.char=';'),
                   growthandeatlengths = utils::read.table(stock[[9]][2],comment.char=';'),
                   doesgrow = as.numeric(stock[[growth.loc]][2]),
                   growth = growth.info(stock[growth.loc:(mort.loc-1)]),
                   naturalmortality = merge.formula(stock[[mort.loc]][-1]),
                   iseaten = as.numeric(stock[[mort.loc+1]][2]),
                   preyinfo = prey.info(stock[(mort.loc+1):(eat.loc-1)]),
                   doeseat = as.numeric(stock[[eat.loc]][2]),
                   predator = pred.info(stock[eat.loc:(pred.loc)]),
                   initialconditions = list(minage = stock[[init.loc + 1]][2],
                                            maxage = stock[[init.loc + 2]][2],
                                            minlength = stock[[init.loc + 3]][2],
                                            maxlength = stock[[init.loc + 4]][2],
                                            dl = ifelse(tolower(stock[[init.loc + 5]][1])=='dl',
                                                        stock[[init.loc + 5]][2],as.numeric(stock[[7]][2])),
                                            sdev = ifelse(tolower(stock[[init.loc + 6]][1])=='sdev',
                                                          stock[[init.loc + 6]][2], 1)),
                   initialdata = initialdata,
                   doesmigrate = as.numeric(stock[[migrate.loc]][2]),
                   yearstep = yearstep,
                   migrationratio = migrationratio,
                   doesmature =  as.numeric(stock[[mature.loc]][2]),
                   maturityfunction = maturity.function,
                   maturestocksandratios = maturestocksandratios,
                   coefficients = coefficients,
                   maturitysteps = maturitysteps,
                   maturitylengths = maturitylengths,
                   doesmove = as.numeric(stock[[move.loc]][2]),
                   transitionstocksandratios = transitionstocksandratios,
                   transitionstep = transitionstep,
                   doesrenew =  as.numeric(stock[[renew.loc]][2]),
                   renewal = list(
                     minlength = ifelse(as.numeric(stock[[renew.loc]][2]) == 0,
                                        '',
                                        as.numeric(stock[[renew.loc + 1]][2])),
                     maxlength = ifelse(as.numeric(stock[[renew.loc]][2]) == 0,
                                        '',
                                        as.numeric(stock[[renew.loc + 2]][2]))),
                   renewal.data = renewal.data,
                   doesspawn = as.numeric(stock[[spawn.loc]][2]),
                   spawning = spawning,
                   doesstray = ifelse(length(stray.loc)==0,
                                      0,as.numeric(stock[[stray.loc]][2]))
      )
    
    return(st)
  }
  stocks <- plyr::llply(stock.files,tmp.func)
  names(stocks) <- plyr::laply(stocks,getStockNames)
  return(stocks)
}
