Modarea <- function(x){
  z=x[1]
  d <- as.numeric(as.character(x[2]))
  L <- as.numeric(as.character(x[3]))
    if(tolower(z)=='a'){ if(d<20) {return(4) } else { return (5) }  }
  if(tolower(z)=='b'){ if(abs(L)<28) {if(d<20) {return(6) } else { return (7) }  } else { if(d<20) {return(2) } else { return (3) } } } 
  if(tolower(z)=='c'){ if(d<20) {return(0) } else { return (1) }  }
}

