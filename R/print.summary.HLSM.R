print.summary.HLSM = function(x,...){
    message("Call:\n")
    print(x$call)
    #message("\n Estimated Intercept:\n")
    #print(x$est.intercept)
    if(all(!is.na(x$est.slopes))){
    message("\n Estimated Slopes:\n")
    print(x$est.slopes) }
}

