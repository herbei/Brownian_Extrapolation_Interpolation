

myalpha <- function(t,a){
   out = pi^2*t / 2 / (a^2)
   return(out)
}

mybeta <- function(t,a){
   out = a^2 / 2 / t
   return(out)
}

myKh <- function(t,a){
   out = sqrt(2) * pi^(3/2) * t^(3/2) / a^2
   return(out)
}


func_n0 <- function(alpha,beta,rho=1/2){
   
   T1 = (1/alpha)*log(4/rho) - 1
   T1 = (1/2)*T1
   
   T2 = (1/beta)*log(1/rho) - 1
   T2 = (1/2)*T2
   
   N = max(T1,T2)
   
   N0 = max( ceiling(N)+1 , 2 ) + 1
   
   return(N0)
   
}


#######################################################
func_g0 = function(y,x){
      
      # Initialize output vector with same length as y
      term1 <- numeric(length(y))
      term2 <- numeric(length(y))
      
      for(i in seq_along(y)){
            term2[i] = 2*x^2*exp(-(y[i]-x)^2 / 2)
            
            if(y[i] >= x){
                  term1[i]=2*x*(y[i]-x) * exp( - (y[i]-x)^2 / 2)
            }else{
                  term1[i]=0
            }
            
            
      }
      out = term1+term2
      return(out)
}
#######################################################

#######################################################
func_f_n = function(y,x,a,n){
      
      
      # Initialize output vector with same length as y
      out <- numeric(length(y))
      
      for(i in seq_along(y)){
            if( (y[i]<0) || (y[i]>a) ){
                  out[i]=0
            }else{
                  term1 = exp(-(y[i] + 2*n*a - x)^2 / 2)
                  term2 = exp(-( y[i]+ 2*n*a + x)^2 / 2)
                  out[i] = (1/sqrt(2*pi)) * (term1 - term2)
            }
            
            
      }
      return(out)
}
#######################################################

#######################################################
simulate_f0 = function(x,a){
      flag = (x>0) & (a > 0) & (2*x <= a)
      if (!flag) stop("conditions for f_0 are not met.")
      
      
      # Region 1: a <= 1
      #------------------------
      if(a <= 1){
            
            sw=0
            while(!sw){
                  U=runif(1,0,1)
                  V=runif(1,0,1)
                  temp = a*sqrt(U)
                  if( V * 2 * x * temp < exp( - (temp-x)^2 / 2 )  - exp( - (temp+x)^2 / 2 ) ){
                        sw=1
                        out=temp
                  }
            }
            
            return(out)
            
      }
      #------------------------
      
      
      # Region 2: x >= 1 ==> a >= 2
      #------------------------
      if(x >= 1){
            sw=0
            while(!sw){
                  V = runif(1,0,1)
                  temp = x + rnorm(1,0,1)
                  if( (0 <= temp) & (temp <= a)  ){
                        if(V  >=  exp(-2*x*temp)){
                              sw=1
                              out=temp
                        }
                  }
            }
            
            return(out)
      }
      #-----------------------
      
      
      # Region 3: x < 1, a > 1 
      #-----------------------
      if( (x<1) & (a>1) ){
            sw=0
            while(!sw){
                  U = runif(1,0,1)
                  V = runif(1,0,1)
                  if(U*(2*x + x^2*sqrt(8*pi)) <= 2*x ){
                        temp = x + sqrt(2*rexp(1))
                  }else{
                        temp = x + rnorm(1,0,1)
                  }
                  Y = V * func_g0(temp,x)
                  
                  if( (0 <= temp) & (temp <= a)  ){
                        if(Y <= exp( -(temp-x)^2 / 2 ) - exp(- (temp+x)^2 / 2) ){
                              sw=1
                              out=temp
                              return(out)
                        }
                        
                  }
                  
            }
            
            
      }
      #--------------------
}
#######################################################

#######################################################
g_series = function(y,t,x,a, Nmax=500){
      
      n <- 1:Nmax
      
      # Initialize output vector with same length as y
      out <- numeric(length(y))
      
      # Pre-compute terms that don't depend on y
      sin_term2 <- sin(pi * n * x / a)
      exp_term <- exp(-n^2 * pi^2 * t / (2 * a^2))
      
      # Vectorized computation for each y value
      for (i in seq_along(y)) {
            if (y[i] < 0 || y[i] > a) {
                  out[i] <- 0
            } else {
                  sin_term1 <- sin(pi * n * y[i] / a)
                  out[i] <- (2 / a) * sum(sin_term1 * sin_term2 * exp_term)
            }
      }
      
      return(out)
      
}

lower_g_bound <- function(y,t,x,a,N){
   
   alpha = myalpha(t,a)
   
   out = numeric(length(y))
   for(i in seq_along(y)){
      T1 = g_series(y[i],t,x,a,N)
      
      T2 = 2*x*pi^2*y[i]/(a^3)
      XX = 2 * (N+1)^2 * exp(-alpha*(N+1)^2)
      out[i] = T1 - T2*XX
      
   }
   
   return(out)
}


upper_g_bound <- function(y,t,x,a,N){
   alpha = myalpha(t,a)
   
   out = numeric(length(y))
   for(i in seq_along(y)){
      T1 = g_series(y[i],t,x,a,N)
      
      T2 = 2*x*pi^2*y[i]/(a^3)
      XX = 2 * (N+1)^2 * exp(-alpha*(N+1)^2)
      out[i] = T1 + T2*XX
      
   }
   
   return(out)
}

#######################################################

#######################################################
# Implementation of Devroye (2009) algorithm
# 
extrapolate_bm_0 = function(x,a){
      flag = (x>0) & (a > 0) & (2*x <= a)
      if (!flag) stop("conditions for f_0 are not met.")
      
      
      # a >= 2
      #-----------------
      if(a >= 2){
            sw=0
            while(!sw){
                  temp = simulate_f0(x,a)
                  V = runif(1,0,1)
                  Y = V * func_f_n(temp,x,a,0)
                  n=1
                  S = func_f_n(temp,x,a,0)
                  decision = "undecided"
                  while( decision == "undecided" ){
                        
                        S = S + func_f_n(temp,x,a,-n)
                        if(Y <= S){ decision="accept" }
                        
                        S = S + func_f_n(temp,x,a,n)
                        if(Y>=S){ decision = "reject" }
                        
                        n = n+1
                  }
                  
                  if(decision=="accept"){
                        sw = 1
                        out = temp
                  }
            }
            
            return(out)
      }
      #-----------------
      
      
      # a<2
      #------------------
      if(a<2){
            rho_term = 4*exp(-3*pi^2/8)
            sw=0
            while(!sw){
                  
                  temp = a*sqrt(runif(1,0,1))
                  V = runif(1,0,1)
                  
                  gg = 2*pi^2*x*temp / (a^3*(1-rho_term))
                  gg = gg * exp(-pi^2 / (2*a^2))
                  Y = V * gg
                  n=1
                  S=0
                  
                  decision = "undecided"
                  while(decision == "undecided"){
                        
                        fn_term = (2/a) * sin(pi*n*temp /a) * sin(pi*n*x/a)* exp( - n^2 * pi^2 / (2*a^2) )
                        S = S + fn_term
                        
                        N = n+1
                        hh = 2*N^2*pi^2*x*temp / (a^3*(1-rho_term))
                        hh = hh * exp(-N^2*pi^2 / (2*a^2))
                        if(Y <= S - hh){ decision = "accept" }
                        if(Y >= S + hh){ decision = "reject" }
                        n = n+1
                  }
                  
                  if(decision == "accept"){
                        sw = 1
                        out = temp
                  }
                  
            }
            
            
            return(out)
            
            
      }
      #------------------         
} # end extrapolate
#######################################################

#######################################################
extrapolate_bm_I <- function(t,x,a){
      if(x <= a/2){
            out = sqrt(t)*extrapolate_bm_0(x/sqrt(t), a/sqrt(t))
      }
      else{
            out = a - sqrt(t)*extrapolate_bm_0((a-x)/sqrt(t), a/sqrt(t))
      }
      return(out)
}
#######################################################



#######################################################


#######################################################
myseries1 = function(alpha, Nmax=500){
      n = 1:Nmax
      out = sum(n^2*exp(-alpha*n^2))
      return(out)
}


myseries2 <- function(beta, Nmax=500){
      n=1:Nmax
      out = sum(exp(-n^2* beta))
      return(out)
}


myseries3 <- function(beta, Nmax = 500){
      n=1:Nmax
      out=sum( (-1)^n * exp(-beta*n^2))
      out = 1 + 2*out
      return(out)
}
#######################################################




#######################################################
q_R = function(y,t,a){
      out = numeric(length(y))
      CC = t*(1-exp(-a^2/(2*t)))
      
      for(i in seq_along(y)){
            if( (y[i]<0) || (y[i])>a  ){
                  out[i]=0
            }else{
                  out[i] =  y[i] * exp(-y[i]^2/(2*t)) / CC
            }
      }
      return(out)
}

q_root_U = function(y,a){
      out = numeric(length(y))
      CC = a^2/2
      
      for(i in seq_along(y)){
            if( (y[i]<0) || (y[i])>a  ){
                  out[i]=0
            }else{
                  out[i] = y[i] / CC
            }
      }
      return(out)
}

q_triangle = function(y,a){
      out=numeric(length(y))
      CC = a^2/4
      
      for(i in seq_along(y)){
            if( (y[i]<0) || (y[i]>a) ){
                  out[i]=0
            }else{
                  out[i] = min(y[i], a-y[i]) / CC
            }
      }
      return(out)
}    
#######################################################


#######################################################
h_series = function(y,t,a, Nmax=500){
      n = 1:Nmax
      Kh = myKh(t,a)
      out = numeric(length(y))
      for(i in seq_along(y)){
            if( (y[i] < 0 ) || (y[i]) >a ){
                  out[i] = 0
            }else{
                  out[i] = Kh * sum(n*sin(n*pi*y[i]/a)*exp(-pi^2 * n^2 * t / 2 /a^2))  
            }
            
      }
      return(out)
}

lower_h_bound <-function(y,t,a,N){
   
   Kh = myKh(t,a)
   alpha = myalpha(t,a)
   
   out = numeric(length(y))
   for(i in seq_along(y)){
      T1 = h_series(y[i],t,a,N)
      
      T2 = Kh*pi*y[i]/a
      XX = 2 * (N+1)^2 * exp(-alpha*(N+1)^2)
      out[i] = T1 - T2*XX
      
   }

   return(out)
}

upper_h_bound <-function(y,t,a,N){
   
   Kh = myKh(t,a)
   alpha = myalpha(t,a)
   
   out = numeric(length(y))
   for(i in seq_along(y)){
      T1 = h_series(y[i],t,a,N)
      
      T2 = Kh*pi*y[i]/a
      XX = 2 * (N+1)^2 * exp(-alpha*(N+1)^2)
      out[i] = T1 + T2*XX
      
   }
   
   return(out)
}
#######################################################



######################################################
lower_D3_bound <-function(t,a,N){
      alpha = myalpha(t,a)
      out = pi^(-1/2) * alpha^(3/2) * a^2 * myseries1(alpha,N+1)
      return(out)
}


upper_D3_bound <-function(t,a,N){
      alpha = myalpha(t,a)
      out = pi^(-1/2) * alpha^(3/2) * a^2 * ( myseries1(alpha,N+1) + 2*(N+1)^2 * exp(-alpha*(N+1)^2) )
      return(out)
}

######################################################




#######################################################
extrapolate_bm_II = function(t,a){
      
      alpha = myalpha(t,a)
      beta = mybeta(t,a)
      alpha0 = 1.086165
      N0 = func_n0(alpha, beta)      
      
      
      
      if(alpha <= alpha0){
            D1 = t*(1-exp(-a^2/(2*t)))
            D2 = (a^2/2) * ( myseries2(beta,N0 + 1) + exp(-beta * (N0+1)^2) )
            p = D1 / (D1+D2)
            
            sw=0
            while(!sw){
                  #--- propose
                  if( runif(1,0,1) < p){
                        U = runif(1,0,1)
                        Y = sqrt( - 2 * t * log(1-U*(1-exp(-a^2/(2*t))) ))
                  }else{
                        U = runif(1,0,1)
                        Y = a*sqrt(U)
                  }
                  
                  
                  #--- accept/reject
                  denominator = D1*q_R(Y,t,a) + D2 * q_root_U(Y,a)
                  U = runif(1,0,1)
                  #if( U* denominator < h_series(temp,t,a)){
                  #      sw=1
                  #      out=temp
                  #      return(out)
                  #}
                  
                  sw1=0
                  N=N0
                  while(!sw1){
                        N = N+50
                        LB = lower_h_bound(Y, t, a, N)
                        UB = upper_h_bound(Y, t, a, N)
                        if(U * denominator < LB){
                              # decision done, accept.
                              sw1 = 1
                              sw = 1
                              out = Y
                              return(out)
                        }
                        if(U*denominator > UB){
                              # decision done, reject
                              sw1 = 1
                        }
                        
                  }
                  
                  
                  
                  
            }
            
      }
      
      
      if(alpha > alpha0){
            
            #D3 = 4*pi^(-1/2)*alpha^(3/2)*myseries1(alpha)*(a^2/4)
            sw=0
            while(!sw){
                  
                  #-- propose
                  U = runif(1,0,1)
                  if(U<1/2){
                        Y = (a/sqrt(2))*sqrt(U)
                  }else{
                        Y = a- (a/sqrt(2))*sqrt(1-U)
                  }
                  #return(temp)
                  
                  #-- accept/reject step
                  U=runif(1,0,1)
                  
                  
                  
                  #denominator = D3 * q_triangle(Y,a)
                  #if(U*denominator < h_series(Y,t,a)){
                  #      sw=1 
                  #      out=Y
                  #      return(out)
                  #}

                  sw1 = 0
                  N = N0
                  while(!sw1){
                        N = N+50
                        LB = lower_h_bound(Y, t, a, N)
                        UB = upper_h_bound(Y, t, a, N)
                        
                        LB_D3 = lower_D3_bound(t,a,N)
                        UB_D3 = upper_D3_bound(t,a,N)
                        if( U * UB_D3 * q_triangle(Y,a) < LB){
                              # decision done, accept
                              sw1 = 1
                              sw = 1
                              out = Y
                              return(out)
                        }
                        if( U * LB_D3 * q_triangle(Y,a) > UB){
                              # decision done, reject
                              sw1 = 1
                        }
                        
                  }
                  
                  
                  
                  
                  
                  
            }
            
      }     
}
#######################################################





#######################################################
# Interpolation functions
#######################################################

func_C1_g <-function(t,x,a){
      delta_term = 4*exp(-3*pi^2/8)
      out = 2* pi^2 *x / a^2 / (1-delta_term)
      out = out * exp( - pi^2 * t / 2 / a^2)
      return(out)
}

func_C2_g <- function(t,x,a){
      out = 1/sqrt(2*pi*t)
      out = out * (1- exp( - (x+a)^2 / (2*t) ))
      return(out)
}


#######################################################
interpolate_BB_I = function(t,T,x,z,a){
      
   
      alpha = myalpha(t,a)
      beta = mybeta(t,a)
      alpha0 = 1.086165
      N0 = func_n0(alpha, beta)      
   
      if( a/sqrt(T-t) < 2 ){
            C <- func_C1_g(T-t,z,a)
      }else{
            C <- func_C2_g(T-t,z,a)
      }      
      
      sw=0
      while(!sw){
            
            # propose 
            Y = extrapolate_bm_I(t,x,a)
            U = runif(1,0,1)
            
            
            sw1 = 0
            N = N0
            while(!sw1){
                 N = N+50
                 LB = lower_g_bound(Y,T-t,z,a,N)
                 UB = upper_g_bound(Y,T-t,z,a,N)
                 if(U*C < LB){
                    # decision done, accept
                    sw1=1
                    sw=1
                    out=Y
                    return(out)
                 }
                 if( U*C > UB ){
                    # decision done, reject
                    sw1=1
                 }
            }
            
      }
}
#######################################################

#######################################################
interpolate_BB_II = function(t,T,z,a){
   
   #interpolate a BB from 0 to z constrained in [0,a]
   
   alpha = myalpha(t,a)
   beta = mybeta(t,a)
   alpha0 = 1.086165
   N0 = func_n0(alpha, beta)   
   
   if( a/sqrt(T-t) <2 ){
      C <- func_C1_g(T-t,z,a)
   }else{
      C <- func_C2_g(T-t,z,a)
   }      
   
   sw=0
   while(!sw){
      
      Y = extrapolate_bm_II(t,a)
      U = runif(1,0,1)
      
      
      sw1 = 0
      N=N0
      while(!sw1){
         N = N+50
         LB = lower_g_bound(Y,T-t,z,a,N)
         UB = upper_g_bound(Y,T-t,z,a,N)
         if(U*C < LB){
            # decision done, accept
            sw1=1
            sw=1
            out=Y
            return(out)
         }
         if( U*C > UB ){
            # decision done, reject
            sw1=1
         }
      }
      
   }
}
#######################################################


#######################################################
interpolate_BB_III = function(t,T,x,a){
   
   #interpolate a BB from x to a, constrained in [0,a]
   
   
   alpha = myalpha(t,a)
   beta = mybeta(t,a)
   alpha0 = 1.086165
   N0 = func_n0(alpha, beta) 
   
   if( a/sqrt(t) < 2 ){
      C <- func_C1_g(t,x,a)
   }else{
      C <- func_C2_g(t,x,a)
   }      
   
   sw=0
   while(!sw){
      
      Y <- extrapolate_bm_II(T-t,a)
      Y <- a-Y
      U <- runif(1,0,1)
      
      # accept / reject test
      sw1 = 0
      N = N0
      while(!sw1){
         N = N+50
         LB = lower_g_bound(Y,t,x,a,N)
         UB = upper_g_bound(Y,t,x,a,N)
         if(U*C < LB){
            # decision done, accept
            sw1=1
            sw=1
            out=Y
            return(out)
         }
         if( U*C > UB ){
            # decision done, reject
            sw1=1
         }
      }
      
   }
}
#######################################################






func_D1_h <- function(t,a){
   
   if(a<sqrt(t)){
      mysup = a * exp(-a^2/2/t)
   }else
   {
      mysup = sqrt(t) * exp(-1/2)
   }
   
   alpha = myalpha(t,a)
   beta = mybeta(t,a)
   N0 = func_n0(alpha,beta)
   
   T2 = myseries2(beta,N0+1)
   T3 = exp(-beta * (N0+1)^2)
   
   out = mysup + a * (T2 + T3)
   return(out)
}


func_D2_h <- function(t,a){
   
   alpha = myalpha(t,a)
   beta = mybeta(t,a)
   N0 = func_n0(alpha,beta)
   
   T1 = 4 * pi^(-1/2) * alpha^(3/2)
   
   T2 = myseries1(alpha, N0+1)
   
   T3 = (N0+1)^2 * exp( - alpha * (N0+1)^2 )
   
   out= T1*(T2 + T3) *(a/2) 
   return(out)
}


#######################################################
interpolate_BB_IV <- function(t,T,a){
   
   
   alpha0 = 1.086165
   alpha = myalpha(T-t,a)
   beta = mybeta(T-t,a)
   N0 = func_n0(alpha,beta)
   
   if(alpha<alpha0){
      D = func_D1_h(T-t,a)
   }else{
      D = func_D2_h(T-t,a)
   }
   
   sw = 0
   while( !sw ){
      
      #propose
      Y = extrapolate_bm_II(t,a)
      U = runif(1,0,1)
      
      
      sw1 = 0
      N=N0
      while(!sw1){
         N=N+50
         LB = lower_h_bound(a-Y,T-t,a,N)
         UB = upper_h_bound(a-Y,T-t,a,N)
         if(U*D < LB){
            sw1=1 # series comparison done
            sw=1 # value us accepted
            out=Y
            return(out)
         }
         if(UB < U*D){
            sw1=1 # series comparison done, proposed value is rejected.
         }
      }
      
      
   } 
}
#######################################################








#######################################################
# Constrained path simulation
#######################################################
sim_bb_path <- function(start_point, end_point, T, a, N=300){

   flag = (start_point>=0) & (start_point < a) & (end_point > 0) & (end_point <= a)
   if (!flag) stop("conditions for f_0 are not met.")
   
   
   
   t_vals <- seq(0,T, length=N)
   dt <- t_vals[2]-t_vals[1]
   PP = numeric(N)
   PP[1] = start_point
   PP[N] <- end_point 
   
   t_value = dt
   
   start_p=start_point

   for(i in 2:(N-1)){
      i1=(start_p==0)
      i2=(end_point==a)
      
      if(i1){
         if(i2){
            PP[i] = interpolate_BB_IV(t_value,T-(i-2)*dt,a)
         }else{
            PP[i] = interpolate_BB_II(t_value,T-(i-2)*dt,end_point,a)
         }
      }else{
         if(i2){
            PP[i] = interpolate_BB_III(t_value,T-(i-2)*dt,start_p,a)
         }else{
            PP[i] = interpolate_BB_I(t_value,T-(i-2)*dt, start_p, end_point,a)
         }
      }
      
      
      start_p = PP[i]   
   }
   

   
   mypath = data.frame(t_vals=t_vals, BB_vals = PP)
      
   return(mypath)
}


