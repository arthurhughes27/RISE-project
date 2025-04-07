test.surrogate = function (full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, 
          szero = NULL, epsilon = NULL, power.want.s = 0.7, u.y.hyp = NULL) 
{
  if (!is.null(full.data)) {
    yone = full.data[full.data[, 3] == 1, 1]
    yzero = full.data[full.data[, 3] == 0, 1]
    sone = full.data[full.data[, 3] == 1, 2]
    szero = full.data[full.data[, 3] == 0, 2]
  }
  dd = delta.calculate(yone = yone, yzero = yzero, sone = sone, 
                       szero = szero)
  z.alpha = qnorm(0.95)
  ci.delta = c(-1, dd$delta.estimate + z.alpha * dd$sd.delta)
  if (is.null(epsilon)) {
    n1 = length(yone)
    n0 = length(yzero)
    sd.null = sqrt((n1 + n0 + 1)/(12 * n1 * n0))
    z.alpha.2 = qnorm(0.975)
    u.s.power = 1/2 - (qnorm(1 - power.want.s) - z.alpha.2) * 
      (sd.null)
    if (is.null(u.y.hyp)) {
      epsilon = dd$u.y - u.s.power
    }
    else epsilon = u.y.hyp - u.s.power
  }
  if (ci.delta[2] < epsilon) {
    is.surrogate = TRUE
  }
  if (ci.delta[2] >= epsilon) {
    is.surrogate = FALSE
  }
  return(list(u.y = as.numeric(dd$u.y), u.s = as.numeric(dd$u.s), 
              delta.estimate = as.numeric(dd$delta.estimate), sd.u.y = dd$sd.u.y, 
              sd.u.s = dd$sd.u.s, sd.delta = as.numeric(dd$sd.delta), 
              ci.delta = ci.delta, epsilon.used = epsilon, is.surrogate = is.surrogate))
}