QRLM <- function (x, y, case.weights = rep(1, nrow(x)), var.weights = rep(1, nrow(x)), ..., w = rep(1, nrow(x)), init = "ls", psi = psi.huber, scale.est = c("MAD", "Huber", "proposal 2"), k2 = 1.345, method = c("M", "MM"), maxit = 20, acc = 1e-04, test.vec = "resid", q = 0.5){

  irls.delta <- function(old, new)
    sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r*w,1,length(r)) %*% x)/sqrt(matrix(w,1,length(r))%*%(x^2))))/sqrt(sum(w*r^2))
  }
  #IRLS digunakan untuk menyelesaikan influence function yang umumnya tidak linear

  method <- match.arg(method)
  nmx <- deparse(substitute(x))

  if(is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)

  if (is.null(colnames(x)))
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x))
    stop("x is singular: singular fits are not implemented in rlm")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec)))
    stop("invalid testvec")
  if (length(var.weights) != nrow(x))
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0))
    stop("Negative var.weights value")
  if (length(case.weights) != nrow(x))
    stop("Length of case.weights must equal number of observations")

  w <- (w * case.weights)/var.weights

  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi))
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0))
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }

    if (is.character(init)) {
      if (init == "ls")
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts")
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init))
        coef <- init$coef
      else coef <- init
      resid <- y - x %*% coef
    }}
  else if (method == "MM") {
    scale.est <- "MM"
    temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...)))
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else stop("method is unknown")
  done <- FALSE
  conv <- NULL
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM")
    scale <- mad(resid/sqrt(var.weights), 0)
  theta <- 2 * pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec))
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD")
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else
          scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights))) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec))
        convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done)
        break
    }
    if (!done)
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
  }
  list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt,
       coefficients = qest)
}


"zerovalinter"<-function(y, x)
{
  if(min(y) > 0) {
    xmin <- x[y == min(y)]
    if(length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  else {
    if(max(y) < 0) {
      xmin <- x[y == max(y)]
      if(length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if(length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if(length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if(length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if(length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
      xmin <- x1
      if(abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <- xzero
  resu
}


"gridfitinter"<-function(y,expectile,Q)
  # computing of the expectile-order of each observation of y by interpolation
{
  nq<-length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile
  vectordest <- apply(diff, 1, zerovalinter,Q)

  #print(vectordest)
  #qord<-list(ord=c(vectordest))
  #qord
}

#' @export
mq=function(y,x,regioncode.s,x.outs,regioncode.r,m,p,SE=0,ydir=0,tol.value=0.0001,maxit.value=100,k.value=1.345) {

  # mengurutkan kode area tersampel
  id.area<- sort(unique(regioncode.s))
  # jumlah area tersampel
  m<-length(id.area)
  # mengurutkan kode area tidak tersampel
  id.area.r=sort(unique(regioncode.r))
  # jumlah area tidak tersampel
  m.r<-length(id.area.r)
  # membuat vektor berisi angka 0 sebanyak area tidak tersampel
  tmp.cont=rep(0,m.r)
  # sintkas dibawah akan menghasilkan nilai 0 untuk area level yang tidak ada pada data susenas
  # dan menghasilkan nilai 1 untuk area level yang ada pada data susenas
  for (i in 1:m.r) {
    for (j in 1:m) {
      if (id.area.r[i]==id.area[j])
        tmp.cont[i]=1
    }
  }
  # which digunakan untuk melihat indeks ke berapa
  # hasil tmp0 berupa indeks2 yang mempunyai nilai 0 pada tmp.cont
  tmp0=which(tmp.cont==0)
  # mengambil kode area dengan nilai tmp0 = 0
  # id.area.out : berisi vektor kode area untuk wilayah tidak tersampel (kretek saja)
  id.area.out=id.area.r[tmp0]
  # id.area.in : berisi vektor kode area untuk wilayah tersampel
  id.area.in=id.area


  # membuat kar.kode.s : berisi kode area format karakter
  # karena klo tidak dibuat karakter pada var wadah.s kode berubah menjadi indeks
  kar.kode.s <- as.character(regioncode.s)
  wadah.s <- as.data.frame(cbind(kar.kode.s,rep(1,length(regioncode.s))))
  wadah.s[,1] <- as.factor(wadah.s[,1])
  names(wadah.s) <- c("kode","jml")
  kar.kode.r <- as.character(regioncode.r)
  wadah.r <- as.data.frame(cbind(kar.kode.r,rep(1,length(regioncode.r))))
  wadah.r[,1] <- as.factor(wadah.r[,1])
  names(wadah.r) <- c("kode","jml")
  # left join digunakan untuk mengambil kode area yang tidak ada pada data susenas
  wadah.baru <- left_join(wadah.r, wadah.s, by = "kode")
  wadah.baru[,2] <- as.numeric(wadah.baru[,2])
  wadah.baru[,3] <- as.numeric(wadah.baru[,3])
  indek <- which(wadah.baru$jml.y==1,)
  # mengambil kode kecamatan yang tidak ada pada data susenas
  kode.daerah.luar <- wadah.baru$kode[-indek]
  kode.daerah.luar <- unique(kode.daerah.luar)
  kode.daerah.luar <- as.character(kode.daerah.luar)
  logika <- regioncode.r == kode.daerah.luar
  # kode.r menghasilkan kode area seperti regioncode.r tapi tanpa kecamatan kretek
  kode.r <- regioncode.r[logika == FALSE]
  # mengatur ulang jumlah faktor
  kode.r <- as.character(kode.r)
  # kode.r berisi kecamatan pada regioncode.s (kecamatan tdk tersampel) tapi tanpa kretek
  #kode.r <- as.factor(kode.r)


  if(m>m.r){
    ember.s <- as.data.frame(cbind(as.character(regioncode.s),rep(1,length(regioncode.s))))
    ember.s[,2] <- as.numeric(ember.s[,2])
    names(ember.s) <- c("kode","jml")
    ember.r <- as.data.frame(cbind(as.character(kode.r),rep(1,length(kode.r))))
    ember.r[,2] <- as.numeric(ember.r[,2])
    names(ember.r) <- c("kode","jml")
    # ember.baru menghasilkan kode area yang jumlah nonsampel = 0
    ember.baru <- anti_join(ember.s, ember.r, by="kode")
    tes <- ember.baru$kode
    tes <- unique(tes)
    cc <- NULL
    # menghasilkan indek2 regioncode.s untuk area yang nonsampel = 0
    for(i in 1:length(tes)){
      log <- regioncode.s == tes[i]
      indek <- which(log==TRUE)
      # cc tidak unik (level unit)
      cc <- c(cc,indek)
    }

    # regioncode.s.baru : region tanpa 20 kecamatan
    regioncode.s.baru <- regioncode.s[-cc]
    #y <- y[-cc]
    #x <- x[-cc,]

    # cc2 : indeks untuk data tersampel yang ada nonsampelnya (unik/ level area)
    cc2 <- as.numeric(table(regioncode.s.baru))
    cc2 <- which(cc2!=0)
  }

  # regioncode.s adalah kode area untuk unit yang tersampel
  #datanew=cbind(y,x,regioncode.s)
  # kode wilayahnya sudah urut
  # ni : berisi area level dan jumlah unitnya untuk area tersampel
  ni=as.numeric(table(regioncode.s))
  # regioncode.r adalah kode area untuk unit tdk tersampel
  # kode wilayahnya sudah urut
  # sample.sizer : berisi area level dan jumlah unitnya untuk area tidak tersampel
  #kodkod <- as.character(kode.r)
  #levels(kodkod) <- levels(regioncode.s)
  #sample.sizer<-as.numeric(table(kodkod))

  if(m!=m.r){
    t <- c(as.character(kode.r),as.character(regioncode.s))
  }else{
    t <- c(as.character(regioncode.r),as.character(regioncode.s))
    kode.r <- regioncode.r
    #t <- c(as.character(kode.r),as.character(regioncode.s))
  }
  datasam <- aggregate(t, list(t), FUN=length)
  aa <- aggregate(regioncode.s, list(regioncode.s), FUN=length)
  wadah <- full_join(datasam, aa, by = "Group.1")
  datasam[,3] <- wadah[,3]
  datasam[,4] <- datasam[,2] - datasam[,3]
  datasam[,1] <- as.numeric(datasam[,1])
  names(datasam) <- c("Kec","pop","sampel","nonsampel")
  datasam <- datasam[order(datasam$Kec),]

  sample.sizer <- as.numeric(datasam$nonsampel)

  # total keseluruhan unit desa dalam bentuk tabel per wilayah kecamatan
  Ni = sample.sizer + ni

  # total seluruh unit desa
  N<-sum(Ni)
  # total seluruh unit desa code yang masuk sampel
  n<-sum(ni)

  # p adalah ukuran atas x+1 (termasuk intersep)
  # membentuk matrix dari variabel penyerta yang tersampel
  x=matrix(x,n,p-1)
  # x.out adalah informasi kovariat atas unit yang tidak tersampel
  # membentuk matrix dari variabel penyerta yang tidak tersampel tanpa kretek (FALSE)
  if(m!=m.r){
    x.outs.r <- x.outs[which(logika==FALSE),]
  }else{
    x.outs.r <- x.outs
  }
  x.r=matrix(x.outs.r,(N-n),p-1)
  # penggabungan baris dari matriks var penyerta tersampel dan tidak tersampel
  x.t=rbind(x,x.r)
  x.c<-rep(1,n)
  # membuat data koefisien (intersep, var penyerta..,)
  x.design<-cbind(x.c,x)
  # jumlah p akan menjadi jmlh var penyerta + 1
  p=ncol(x.design)

  ob<-QRLM(x.design, y,q=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98)),k=k.value,maxit=maxit.value,acc=tol.value)

  qo<-matrix(c(gridfitinter(y,ob$fitted.values,ob$q.values)),nrow=n,ncol=1)

  qmat<-matrix(c(qo,regioncode.s),nrow=sum(ni),ncol=2)

  mqo<-aggregate(qmat[,1],list(d2=qmat[,2]),mean)[,2]

  # mengurutkan dan mengunikan kode wilayah tersampel
  uar<-sort(unique(regioncode.s))
  saq<-matrix(c(mqo,uar),nrow=m,ncol=2)

  saq<-rbind(saq,c(0.5,9999))

  beta.stored=matrix(0,m,2)
  res.s=NULL
  tttmp1<-NULL
  ci<-array(rep(0,n*m),dim=c(n,m,1))
  ci1<-array(rep(0,n*m),dim=c(n,m,1))
  prs<-NULL
  prr<-NULL
  wi<-matrix(0,n,m)

  # m adalah jumlah area kecilnya
  MQNAIVE<-c(rep(0,m))

  if(m>m.r){
    for(i in cc2){
      ob1<-QRLM(x.design,y,q=mqo[i],psi=psi.huber,k = k.value,maxit=maxit.value,acc=tol.value)
      wd<-diag(c(ob1$q.weights))
      # Regional parameters from multiquantile model
      # need to be ordered by area
      coef<-matrix(c(t(ob1$coefficients)),nrow=1,ncol=p)
      coef<-t(coef)
      meat<-wd%*%x.design%*%solve(t(x.design)%*%wd%*%x.design)
      x1<-c(rep(1,(Ni[i]-ni[i])))
      ir<-rep(0,n)
      ir[regioncode.s==uar[i]]<-1
      rj1<-sample.sizer[i]
      # berbentuk vektor, berisi jumlah var penyerta tdk tersampel di setiap wilayah kecamatan
      r=NULL
      for (kk in 1:(p-1)) {
        r<-c(r,sum(x.r[,kk][kode.r==uar[i]]))
      }
      r=c(rj1,r)
      sj1<-sum(rep(1,ni[i]))
      tss=NULL
      # berbentuk vektor, berisi jumlah var penyerta tersampel di setiap wilayah kecamatan
      for (kk in 1:(p-1)) {
        tss<-c(tss,sum(x[,kk][regioncode.s==uar[i]]))
      }
      tss<-c(sj1,tss)
      w.welsh<-((Ni[i])/(ni[i]))*ir+meat%*%(r-((Ni[i]-ni[i])/ni[i])*tss)
      y.i<-y[regioncode.s==uar[i]]
      y.pred.s<-cbind(1,x[regioncode.s==uar[i],,drop = FALSE])%*%coef
      residual<-y.i-y.pred.s
      tttmp1[i]<-(sample.sizer[i]/ni[i])*sum(residual)
      prs<-c(prs,y.pred.s)
      res.s<-c(res.s,residual)
      y.pred<-cbind(1,x.r[kode.r==uar[i],,drop = FALSE])%*%coef
      prr<-c(prr,y.pred)
      data<-cbind(regioncode.s,w.welsh)
      for (kk in 1:n){
        if (data[kk,1]==uar[i])
          ci[kk,i,1]<-data[kk,2]-1
        else if (data[kk,1]!=uar[i])
          ci[kk,i,1]<-data[kk,2]
      }
      MQNAIVE[i]<-(1/Ni[i])*as.numeric(sum(y.i)+sum(y.pred))
      ai<-r
      bi<-ai%*%solve(t(x.design)%*%wd%*%x.design)%*%t(x.design)%*%wd
      bi<-c(bi)
      wi[,i]<-c(ir+bi)
      datanaive<-cbind(regioncode.s,wi[,i])
      for (kk in 1:n) {
        if (datanaive[kk,1]==uar[i])
          ci1[kk,i,1]<-datanaive[kk,2]-1
        else if (datanaive[kk,1]!=uar[i])
          ci1[kk,i,1]<-datanaive[kk,2]
      }
    }
    mm <- c(1:m)
    mm <- mm[-cc2]
    for(i in mm){
      MQNAIVE[i]<-ydir[i]
    }
  }else{
    for(i in 1:m){
      ob1<-QRLM(x.design,y,q=mqo[i],psi=psi.huber,k = k.value,maxit=maxit.value,acc=tol.value)
      wd<-diag(c(ob1$q.weights))
      # Regional parameters from multiquantile model
      # need to be ordered by area
      coef<-matrix(c(t(ob1$coefficients)),nrow=1,ncol=p)
      coef<-t(coef)
      meat<-wd%*%x.design%*%solve(t(x.design)%*%wd%*%x.design)
      x1<-c(rep(1,(Ni[i]-ni[i])))
      ir<-rep(0,n)
      ir[regioncode.s==uar[i]]<-1
      rj1<-sample.sizer[i]
      # berbentuk vektor, berisi jumlah var penyerta tdk tersampel di setiap wilayah kecamatan
      r=NULL
      for (kk in 1:(p-1)) {
        r<-c(r,sum(x.r[,kk][kode.r==uar[i]]))
      }
      r=c(rj1,r)
      sj1<-sum(rep(1,ni[i]))
      tss=NULL
      # berbentuk vektor, berisi jumlah var penyerta tersampel di setiap wilayah kecamatan
      for (kk in 1:(p-1)) {
        tss<-c(tss,sum(x[,kk][regioncode.s==uar[i]]))
      }
      tss<-c(sj1,tss)
      w.welsh<-((Ni[i])/(ni[i]))*ir+meat%*%(r-((Ni[i]-ni[i])/ni[i])*tss)
      y.i<-y[regioncode.s==uar[i]]
      y.pred.s<-cbind(1,x[regioncode.s==uar[i],,drop = FALSE])%*%coef
      residual<-y.i-y.pred.s
      tttmp1[i]<-(sample.sizer[i]/ni[i])*sum(residual)
      prs<-c(prs,y.pred.s)
      res.s<-c(res.s,residual)
      y.pred<-cbind(1,x.r[kode.r==uar[i],,drop = FALSE])%*%coef
      prr<-c(prr,y.pred)
      data<-cbind(regioncode.s,w.welsh)
      for (kk in 1:n){
        if (data[kk,1]==uar[i])
          ci[kk,i,1]<-data[kk,2]-1
        else if (data[kk,1]!=uar[i])
          ci[kk,i,1]<-data[kk,2]
      }
      MQNAIVE[i]<-(1/Ni[i])*as.numeric(sum(y.i)+sum(y.pred))
      ai<-r
      bi<-ai%*%solve(t(x.design)%*%wd%*%x.design)%*%t(x.design)%*%wd
      bi<-c(bi)
      wi[,i]<-c(ir+bi)
      datanaive<-cbind(regioncode.s,wi[,i])
      for (kk in 1:n) {
        if (datanaive[kk,1]==uar[i])
          ci1[kk,i,1]<-datanaive[kk,2]-1
        else if (datanaive[kk,1]!=uar[i])
          ci1[kk,i,1]<-datanaive[kk,2]
      }
    }
  }

  cek <- m-m.r
  # Sintetik MQuantile
  # out of sample area (1 kecamatan, kretek)
  if (length(tmp.cont)!=sum(tmp.cont)){
    # rr.sample adalah jumlah yang tersampel pada data susenas
    rr.sample=length(id.area.out)
    MQ.Mean.out <- c()
    for (i in 1:rr.sample)
    {
      # Ri : jumlah area yang tidak ada pada data
      Ri=sum(unique(regioncode.r)==id.area.out[i])
      x.r.area<-matrix(cbind(1,x.outs)[regioncode.r==id.area.out[i]],Ri,p)
      pred.medr=0
      for (j in 1:Ri){
        ob1=QRLM(x.design, y,maxit = 100,q=0.5)
        coef<-matrix(c(t(ob1$coef)),nrow=1,ncol=p)
        # need to be ordered by area
        coef<-t(coef)
        pred.medr[j]<-(x.r.area[j,]%*%coef[,1])
      }
      MQ.Mean.out[i]<-1/(Ri)*(sum(pred.medr))}
  }

  #Evaluasi statistik
  #res.d1<-cbind(res.s^2,regioncode.s,ci1[,,1])
  v1<-NULL
  bias<-NULL
  mse<-NULL
  if(m>m.r){
    qq <- c(1:256)
    qq <- qq[-cc]
    apple <- c(1:256)
    prs2 <- c(1:256)
    wi[cc,] <- 0
    ci1[cc,,1] <- 0
    for(i in 1:(length(res.s))){
      apple[qq[i]] <- res.s[i]
      prs2[qq[i]] <- prs[i]
    }
    for(i in 1:length(cc)){
      apple[cc[i]] <- 0
      prs2[cc[i]] <- 0
    }
    res.d1<-cbind(apple^2,regioncode.s,ci1[,,1])
    for(oo in cc2){
      v1[oo]<-1/Ni[oo]^2*(sum((res.d1[,(oo+2)][res.d1[,2]==uar[oo]]^2+(sample.sizer[oo])/n)*res.d1[,1][res.d1[,2]==uar[oo]])+sum(res.d1[,(oo+2)][res.d1[,2]!=uar[oo]]^2*res.d1[,1][res.d1[,2]!=uar[oo]]))

      bias[oo]<-(1/Ni[oo])*(sum(wi[,oo]*prs2)-sum(c(prs2[regioncode.s==uar[oo]],prr[kode.r==uar[oo]])))
      mse[oo]<-v1[oo]+(bias[oo])^2
    }
    for(oo in mm){
      mse[oo]<-(SE[oo]/MQNAIVE[oo]*100)^2
    }
  }else{
    res.d1<-cbind(res.s^2,regioncode.s,ci1[,,1])
    for (oo in 1:m){
      v1[oo]<-1/Ni[oo]^2*(sum((res.d1[,(oo+2)][res.d1[,2]==uar[oo]]^2+(sample.sizer[oo])/n)*res.d1[,1][res.d1[,2]==uar[oo]])+sum(res.d1[,(oo+2)][res.d1[,2]!=uar[oo]]^2*res.d1[,1][res.d1[,2]!=uar[oo]]))

      bias[oo]<-(1/Ni[oo])*(sum(wi[,oo]*prs)-sum(c(prs[regioncode.s==uar[oo]],prr[kode.r==uar[oo]])))
      mse[oo]<-v1[oo]+(bias[oo])^2
    }
  }


  if(m!=m.r){
    list(mq.naive = MQNAIVE,
         area.tersampel = uar,
         mq.nirsampel = MQ.Mean.out,
         area.nirsampel = id.area.out,
         mse.tersampel = mse
         #mse.naive.nirsampel = median(mse)
    )
  }else{
    list(mq.naive = MQNAIVE,
         area.tersampel = uar,
         mse.tersampel = mse)
  }


}
