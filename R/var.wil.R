#a is one number, b is vector
var.wil = function(a, b, flip = FALSE){
  if(!flip) return(mean(sapply(b, my.wilcox, x= a)))
  if(flip) return(mean(sapply(b, my.wilcox, y= a)))
}