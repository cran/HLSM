adjust.my.tune=function(tune, acc, dim){
    n=length(acc)
    if(dim==1){
        for(i in 1:n){
            if(acc[i] <= 0.48 & acc[i] >= 0.4){tune[i] = tune[i]
            } else if(acc[i]>0.9){tune[i]=tune[i]*3
            } else if(acc[i]>0.8){tune[i]=tune[i]*2
            } else if(acc[i]>0.7){tune[i]=tune[i]*1.5
            } else if(acc[i]>0.6){tune[i]=tune[i]*1.25
            } else if(acc[i]>0.5){tune[i]=tune[i]*1.1
            } else if(acc[i]>0.48){tune[i]=tune[i]*1.05
            } else if(acc[i]<0.1){tune[i]=tune[i]*0.5
            } else if(acc[i]<0.25){tune[i]=tune[i]*0.75
            } else if(acc[i]<0.35){tune[i]=tune[i]*0.9
            } else if(acc[i]<0.4){tune[i]=tune[i]*.95
            }# else{}
        }}
    
    if(dim>1){
        for(i in 1:2){
        if(acc[i] <= 0.48 & acc[i] >= 0.4){tune[i] = tune[i]
        } else if(acc[i]>0.9){tune[i]=tune[i]*3
        } else if(acc[i]>0.8){tune[i]=tune[i]*2
        } else if(acc[i]>0.7){tune[i]=tune[i]*1.5
        } else if(acc[i]>0.6){tune[i]=tune[i]*1.25
        } else if(acc[i]>0.5){tune[i]=tune[i]*1.1
        } else if(acc[i]>0.48){tune[i]=tune[i]*1.05
        } else if(acc[i]<0.1){tune[i]=tune[i]*0.5
         } else if(acc[i]<0.25){tune[i]=tune[i]*0.75
        } else if(acc[i]<0.35){tune[i]=tune[i]*0.9
        } else if(acc[i]<0.4){tune[i]=tune[i]*.95
        }#else{}
        #tune[i] = min(tune[i],20)
        }
        for(i in 1:n){
            if(acc[i] <= 0.2 & acc[i] >= 0.26){tune[i] = tune[i]
            } else if(acc[i]>0.9){tune[i]=tune[i]*3
            } else if(acc[i]>0.8){tune[i]=tune[i]*2.5
            } else if(acc[i]>0.7){tune[i]=tune[i]*2
            } else if(acc[i]>0.55){tune[i]=tune[i]*1.5
            } else if(acc[i]>0.4){tune[i]=tune[i]*1.25
            } else if(acc[i]>0.3){tune[i]=tune[i]*1.1
            } else if(acc[i]>0.26){tune[i]=tune[i]*1.05
            } else if(acc[i]<0.05){tune[i]=tune[i]*0.75
            } else if(acc[i]<0.15){tune[i]=tune[i]*0.8
            } else if(acc[i]<0.18){tune[i]=tune[i]*0.85
            } else if(acc[i]<0.2){tune[i]=tune[i]*.9
            } # else{}
        }
    }
    return(tune)
}



poor.acceptance=function(vector, dim){
    if(dim>1){which(vector<0.10 | vector>0.4)}
    else{which(vector<0.30 | vector>0.55)}
}


which.suck=function(vectors, dim){
 if(dim>1){which(vectors<0.10 | vectors>0.4)
	}else{which(vectors<0.30 | vectors>0.55)}

}
