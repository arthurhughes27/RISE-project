delta.calculate = function (full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, 
          szero = NULL) 
{
  if (!is.null(full.data)) {
    yone = full.data[full.data[, 3] == 1, 1]
    yzero = full.data[full.data[, 3] == 0, 1]
    sone = full.data[full.data[, 3] == 1, 2]
    szero = full.data[full.data[, 3] == 0, 2]
  }
  test.y = wilcox.test(yone, yzero, exact = F)
  test.s = wilcox.test(sone, szero, exact = F)
  n1.f = length(yone)
  n0.f = length(yzero)
  u.y = (n1.f * n0.f)^(-1) * test.y$statistic
  u.s = (n1.f * n0.f)^(-1) * test.s$statistic
  delta.estimate = u.y - u.s
  m.count = length(yone)
  n.count = length(yzero)
  V10.Xi.Y = sapply(yone, var.wil, b = yzero)
  V10.Xi.S = sapply(sone, var.wil, b = szero)
  V01.Yj.Y = sapply(yzero, var.wil, b = yone, flip = TRUE)
  V01.Yj.S = sapply(szero, var.wil, b = sone, flip = TRUE)
  s10.11.YY = 1/(m.count - 1) * sum((V10.Xi.Y - u.y) * (V10.Xi.Y - 
                                                          u.y))
  s10.12.YS = 1/(m.count - 1) * sum((V10.Xi.Y - u.y) * (V10.Xi.S - 
                                                          u.s))
  s10.22.SS = 1/(m.count - 1) * sum((V10.Xi.S - u.s) * (V10.Xi.S - 
                                                          u.s))
  s10.21.SY = 1/(m.count - 1) * sum((V10.Xi.Y - u.y) * (V10.Xi.S - 
                                                          u.s))
  s01.11.YY = 1/(n.count - 1) * sum((V01.Yj.Y - u.y) * (V01.Yj.Y - 
                                                          u.y))
  s01.12.YS = 1/(n.count - 1) * sum((V01.Yj.Y - u.y) * (V01.Yj.S - 
                                                          u.s))
  s01.22.SS = 1/(n.count - 1) * sum((V01.Yj.S - u.s) * (V01.Yj.S - 
                                                          u.s))
  s01.21.SY = 1/(n.count - 1) * sum((V01.Yj.Y - u.y) * (V01.Yj.S - 
                                                          u.s))
  S10 = matrix(c(s10.11.YY, s10.12.YS, s10.21.SY, s10.22.SS), 
               nrow = 2, ncol = 2, byrow = TRUE)
  S01 = matrix(c(s01.11.YY, s01.12.YS, s01.21.SY, s01.22.SS), 
               nrow = 2, ncol = 2, byrow = TRUE)
  S.mat = (1/m.count) * S10 + (1/n.count) * S01
  sd.y = sqrt(S.mat[1, 1])
  sd.s = sqrt(S.mat[2, 2])
  L = t(as.matrix(c(1, -1)))
  sd.est = (L %*% ((1/m.count) * S10 + (1/n.count) * S01) %*% 
              t(L))^(1/2)
  return(list(u.y = as.numeric(u.y), u.s = as.numeric(u.s), 
              delta.estimate = as.numeric(delta.estimate), sd.u.y = sd.y, 
              sd.u.s = sd.s, sd.delta = as.numeric(sd.est)))
}