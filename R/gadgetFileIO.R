##' @title gadget.fit helper functions
##' @description \code{eval.gadget.formula} Evaluate gadget formulas, which are in reverse polish notation, ie
##' '(* x y)' which is equivalent to 'x*y'. The evaluation supports the following
##' symbols '*','/','+','-','exp','log','sqrt'. The evaluation uses a gadget
##' parameter object for its evaluation.
##' @rdname gadgetFitHelpers
##' @param gad.for gadget formula
##' @param par gadget parameters object
##' @return a vector of evaluated gadget formulas
##' @author Bjarki Thor Elvarsson
eval.gadget.formula <- function(gad.for,par=data.frame()){
  tmp <- strsplit(gsub(')',' )',gsub('(','',gad.for,fixed=TRUE)),' ')
  plyr::ldply(tmp,
              function(x){
                x <- x[!x=='']
                x[x=='-'] <- "'-'("
                x[x=='+'] <- "'+'("
                par.ind <- grep('#',x,fixed=TRUE)
                x <- gsub("*","prod(",x,fixed=TRUE)
                x <- gsub("/","'/'(",x,fixed=TRUE)
                x <- gsub("+ ","sum(",x,fixed=TRUE)
                x <- gsub("- ","'-'(",x,fixed=TRUE)
                x <- gsub('exp','exp(',x,fixed = TRUE)
                x <- gsub('log','log(',x,fixed = TRUE)
                x <- gsub('sqrt','sqrt(',x,fixed = TRUE)
                ## remove initial values
                x <- gsub('[0-9]+.[0-9]+#|[0-9]+#','#',x)
                x[par.ind] <- par[gsub('#','',x[par.ind],fixed=TRUE),'value']
                x <- gsub(',)',')',gsub('(,','(',paste(x,collapse=','),fixed=TRUE),
                          fixed=TRUE)
                return(eval(parse(text=x)))
              })
}

##' @rdname gadgetFitHelpers
##' @description \code{merge.formula} merges txt to gadget.formula
##' @param txt gadget formula text
##' @return txt where the formula has been properly compiled
merge.formula <- function(txt){
  openP <- grep('(',txt,fixed=TRUE)
  closeP <- grep(')',txt,fixed=TRUE)
  if(length(openP) + length(closeP) == 0)
    return(txt)
  
  if(length(openP) != length(closeP))
    stop('numbers of paranthesis dont match in gadget formula')
  
  braces <- data.frame(begin=openP,end=closeP,group=openP)
  for(i in 1:length(openP)){
    braces$end[i] <- closeP[i]
    braces$begin[i] <- openP[max(which(openP < closeP[i]))]
    openP[max(which(openP < closeP[i]))] <- length(txt)
  }
  braces <- dplyr::arrange(braces, .data$begin)
  for(i in 1:length(openP)){
    braces$group[braces$end<braces$end[i] & braces$end>braces$begin[i]] <-
      braces$group[i]
  }
  
  braces <- plyr::ddply(braces,'group',function(x) utils::head(x,1))
  for(i in length(braces$group):1){
    txt[braces$begin[i]] <- paste(txt[braces$begin[i]:braces$end[i]],
                                  collapse=' ')
    txt <- txt[-c((braces$begin[i]+1):braces$end[i])]
  }
  return(txt)
}

##' \code{read.gadget.table} reads gadget tables
##' @rdname gadgetFileIO
##' @param file path to file
##' @param header logical, should the header be read from the file
##' @return data.frame 
read.gadget.table <- function(file,header=FALSE){
  dat <- strip.comments(file)
  if(class(dat) == 'list')
    gad.tab <- plyr::ldply(dat,merge.formula)
  else {
    gad.tab <- plyr::adply(dat,2,merge.formula)
    gad.tab$X1 <- NULL
  }
  if(header){
    comments <- attr(dat,'comments')
    header <- utils::tail(comments,1)
    ## unfinised business
  }
  
  return(gad.tab)
}



##' @rdname gadgetFileIO
##' @description \code{read.gadget.fleet} read gadget fleet
##' @param fleet.file location of the fleet file
##' @return fleet file object
read.gadget.fleet <- function(fleet.file='fleet'){
  fleet <- strip.comments(fleet.file)
  comp.loc <- grep('fleetcomponent|component',fleet)
  suit.loc <- grep('suitability',fleet)
  fleet.dat <-
    data.frame(fleet = plyr::laply(fleet[comp.loc+1],function(x) x[2]),
               type = plyr::laply(fleet[comp.loc+1],function(x) x[1]),
               livesonareas = plyr::laply(fleet[comp.loc+2],
                                          function(x) paste(x[-1],collapse=' ')),
               multiplicative = plyr::laply(fleet[comp.loc+3],
                                            function(x) as.numeric(x[2])),
               amount =  plyr::laply(fleet[c(comp.loc[-1]-1,
                                             length(fleet))],
                                     function(x) x[2]),
               stringsAsFactors=FALSE
    )
  diff.suit <- data.frame(fleet=plyr::laply(fleet[comp.loc+1],function(x) x[2]),
                          begin=suit.loc+1,
                          end=c(comp.loc[-1]-2,length(fleet)-1))
  prey <- plyr::ddply(diff.suit,'fleet',
                      function(x){
                        plyr::ldply(fleet[x$begin:x$end],
                                    function(x)
                                      c(stock=x[1],suitability=x[3],
                                        params=paste(utils::tail(x,-3),collapse=' ')))
                        
                      })
  return(list(fleet=fleet.dat,prey=prey))
}

##' \code{write.gadget.fleet} writes gadget fleet DEPRECATED
##' @rdname gadgetFileIO
##' @param fleet a gadget.fleet 
##' @param file location of the file
##' @return nothing
write.gadget.fleet <- function(fleet,file='fleet'){
  base.text <-
    paste('[fleetcomponent]',
          '%s\t%s',
          'livesonareas\t%s',
          'multiplicative\t%s',
          'suitability',
          '%s',
          'amount\t%s',
          sep='\n')
  
  suit.text <- plyr::ddply(fleet$prey,'fleet',
                           function(x){
                             c(suitability=paste(x$stock,'function',
                                                 x$suitability,x$params,
                                                 sep='\t', collapse='\n'))
                           })
  tmp <- merge(fleet$fleet,suit.text,by='fleet')
  tmp$suitability <- ifelse(tmp$type=='quotafleet',
                            paste(tmp$suitability,
                                  sprintf('quotafunction\t%s\nbiomasslevel\t%s\nquotalevel\t%s\nselectstocks\t%s',
                                          tmp$quotafunction,tmp$biomasslevel,
                                          tmp$quotalevel,
                                          tmp$selectstocks),
                                  sep='\n'),
                            tmp$suitability)
  
  write.unix(sprintf(base.text,tmp$type,tmp$fleet,tmp$livesonareas,
                     tmp$multiplicative,tmp$suitability, tmp$amount),
             f=file)
  
}

##' @rdname gadgetFileIO
##' @description \code{read.gadget.area} reads the gadget area file
##' @param area.file location of the area file
##' @return areafile object
##' @author Bjarki Thor Elvarsson
read.gadget.area <- function(area.file='area'){
  area <- strip.comments(area.file)
  areas <- area[[1]][-1]
  size <- area[[2]][-1]
  temperature <-
    as.data.frame(t(sapply(area[-c(1:3)],function(x) as.numeric(x))))
  names(temperature) <- c('year','step','area','temperature')
  area <- list(areas=areas,size=size,temperature=temperature)
  class(area) <- c('gadget.area',class(area))
  return(area)
}



##' @rdname gadgetFileIO
##' @description \code{write.gadget.area} writes a gadget area file
##' @param area data frame with area and temperature
##' @param file location of the area file
##' @return nothing
##' @author Bjarki Thor Elvarsson
write.gadget.area <- function(area,file='area'){
  header <- sprintf('; time file created in Rgadget\n; %s - %s',file,date())
  area.file <-
    paste(header,
          paste('areas',paste(area$areas,collapse=' '),sep='\t'),
          paste('size',paste(area$size,collapse=' '),sep='\t'),
          'temperature',
          '; year - step - area - temperature',
          sep='\n')
  write.unix(area.file,f=file)
  write.gadget.table(area$temperature,file=file,col.names=FALSE,append=TRUE,
                     quote=FALSE,sep='\t',row.names=FALSE)
}

##' \code{read.gadget.time} reads gadget time file
##' @rdname gadgetFileIO
##' @param time.file location of the file
##' @return timefile objects
read.gadget.time <- function(time.file='time'){
  time <- strip.comments(time.file)
  time.names <- sapply(time,function(x) x[1])
  time <- sapply(time,function(x) as.numeric(x[-1]))
  names(time) <- time.names
  if(sum(time$notimesteps[-1])!=12)
    warning('Error in timefile - notimesteps does not sum to 12')
  if(length(time$notimesteps[-1])!=time$notimesteps[1])
    warning('Error in timefile - notimesteps does not contain the right number of timesteps')
  time$notimesteps <- time$notimesteps[-1]
  class(time) <- c('gadget.time',class(time))
  return(time)
}

##' \code{write.gadget.time} writes gadget time file
##' @rdname gadgetFileIO
##' @param time data.frame with time
##' @param file location of the time file
##' @return nothing
write.gadget.time <- function(time,file='time'){
  header <- sprintf('; time file created in Rgadget\n; %s - %s',file,date())
  time.file <-
    paste(header,
          paste('firstyear',time$firstyear,sep='\t'),
          paste('firststep',time$firststep,sep='\t'),
          paste('lastyear',time$lastyear,sep='\t'),
          paste('laststep',time$laststep,sep='\t'),
          paste('notimesteps',
                paste(length(time$notimesteps),
                      paste(time$notimesteps,collapse=' ')),
                sep='\t'),
          sep='\n')
  write.unix(time.file,f=file)
}



