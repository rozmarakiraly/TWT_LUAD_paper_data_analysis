two_group_test <- function(x,y,log_data) {
  normality <- shapiro.test(x)[[2]] > 0.05 & shapiro.test(y)[[2]] > 0.05
  variance_equality <- lawstat::levene.test(c(x,y),c(rep("x",length(x)),rep("y",length(y))))[[2]] > 0.05
  if (normality & variance_equality) {
    temp_tt = list(p_value=t.test(x,y,var.equal=T)[[3]],
                   test_type="StudentT",
                   fold_change=if(log_data){mean(x,na.rm=T)-mean(y,na.rm=T)}else{mean(x,na.rm=T)/mean(y,na.rm=T)})
    return(temp_tt)
  } else if (normality) {
    temp_tt = list(p_value=t.test(x,y,var.equal=F)[[3]],
                   test_type="WelchT",
                   fold_change=if(log_data){mean(x,na.rm=T)-mean(y,na.rm=T)}else{mean(x,na.rm=T)/mean(y,na.rm=T)})
    return(temp_tt)
  } else {
    temp_tt = list(p_value=wilcox.test(x,y)[[3]],
                   test_type="Wilcoxon",
                   fold_change=if(log_data){mean(x,na.rm=T)-mean(y,na.rm=T)}else{mean(x,na.rm=T)/mean(y,na.rm=T)})
    return(temp_tt)
  }
}
