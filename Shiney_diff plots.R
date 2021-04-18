## R shiny app to Simulate Disease Seasonality Using Temporally Forced Models
## and Test how seasonality changes with R0 and amplitude of seasonality
## Model: SEIR; Parameters based on childhood infections
## 4/15/18, by Wan Yang

library(shiny)
library(deSolve)


HIVTB = function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = caplambda - (beta*c*S*(I+J3)/N) - lambda*sigma*S*(Jas/R) - mu*S
    dL = beta*c * (S + Tb)*(I+J3)/N - lambda*sigma*L*(Jas/R) - ( mu+k+r1)*L
    dI = k*L - (mu+d+r2)*I
    dTb = r1*L + r2*I - beta*c*Tb*(I+J3)/N - lambda*sigma*Tb*(Jas/R) - mu*Tb
    dJ1 = lambda*sigma*(S+Tb)*(Jas/R) - beta*c*J1*(I+J3)/N - (alpha1 + mu)*J1 + ras*J2
    dJ2 = lambda*sigma*L*(Jas/R) - beta*c*J1*(I+J3)/N - (alpha1 + mu+kas+ras)*J2
    dJ3 = kas*J2 - (alpha3 + mu + das)*J3
    dA = alpha1*J1 + alpha2*J2 + alpha3*J3 - (mu+f)*A
    list(c(dS, dL, dI, dTb, dJ1, dJ2, dJ3, dA))
  })
}
N=1; # 100%

ui <- fluidPage(h4('Simulating TB-HIV Co-infections with changes in R0 and beta??'),
                h5('Note: it may take a while to run... wait... wait... wait...'),
                # outputs
                sidebarLayout(
                  sidebarPanel(width=4,
                               sliderInput(inputId = 'R0',h6('R0 for TB'),value=1.7,min=0,max=3,step=.1),
                               sliderInput(inputId = 'lambda',h6('Transmission rate for HIV (lambda)'),value=0.02,min=0,max=0.25,step=.01),
                               #sliderInput(inputId = 'lambda11',h6('Transmission rate for HIV (lambda) for Run 1'),value=0.02,min=0,max=0.25,step=.01),
                               #sliderInput(inputId = 'lambda12',h6('Transmission rate for HIV (lambda) for Run 2'),value=0.1,min=0,max=0.25,step=.01),
                               #sliderInput(inputId = 'lambda13',h6('Transmission rate for HIV (lambda) for Run 3'),value=0.23,min=0,max=0.25,step=.01),
                               h6('Note: other parameters based on childhood infections: \nlatent period = 8 days; infectious period = 5 days')
                  ),
                  
                  mainPanel(
                    plotOutput(outputId = 'plots',width = '100%', height = "550px")
                  )
                )
                
)

server <- function(input, output){
  
  output$plots=renderPlot({
    
    R0=input$R0         # R0 for all three runs
    lambda=input$lambda

    
    
    # other parameters and initial conditions
    TN = 8336817; N = .825 * TN; L = 0.16 * N; I = 6.9/100000 * N;
    Tb = 0.82 * I; J1 = 1458/100000 *N; J2 = .14*J1; J3 = 0.6 * I;
    A = 0.166 * J1; gamma = .95; beta = R0 * gamma; 
    #lambda = 4.025/100000
    #R0 = 1.7
    
    state = c(
      S = N - L - J1 - J2 - J3 , L = L, I = I, Tb = Tb, J1 = J1,
      J2 = J2, J3 = J3, A = A
    )
    
    
    params = c(
      N = N,
      R = R,
      Jas = J1 + J2 + J3,
      caplambda = (.175*TN)/15,
      beta = beta,
      lambda = lambda,
      c = (10 * I) / R,
      sigma = 0.009/100000 ,
      mu = 555.1/100000 ,
      k = 0.08,
      kas = 0.1,
      d = 0.9,
      das = 4.7/1000 ,
      f = 21.35/100000 ,
      r1 = 0.78 ,
      r2 = 0.82 ,
      ras = 1,
      alpha1 = 10/100000 ,
      alpha2 = 0.25 ,
      alpha3 = 0.25
    )

    tm_step=7; # run it by week
    times=seq(0,10,by=tm_step/365) # need to run for ~10 yrs to have the 3rd sim stable!
    
    sim=ode(y=state,times=times,func=HIVTB,parms=params)

    
    # plot results
    par(mfrow=c(3,1),cex=1,mar=c(3,3,1,1),oma=c(0,0,1,0),mgp=c(1.8,.5,0))
    plot(sim[,'time'],(sim[,'I']),ylab='TB infectious',xlab='Time (year)',
         type='l',lwd=2, ylim=c(-5,20000))
    mtext(bquote(lambda[1]==.(lambda)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    mtext(bquote(R[0]==.(R0)),side=3,line=0,outer=F,adj=.5,cex=1.5)
    abline(v=990:1000,col='grey',lty=2)
    plot(sim[,'time'],(sim[,'J1']),ylab='HIV infectious',xlab='Time (year)',
         type='l',lwd=2, ylim=c(0,150000))
    mtext(bquote(lambda[1]==.(lambda)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    abline(v=990:1000,col='grey',lty=2)
    plot(sim[,'time'],(sim[,'J3']),ylab='TB-HIV infectious',xlab='Time (year)',
         type='l',lwd=2, ylim=c(-5,2000))
    mtext(bquote(lambda[1]==.(lambda)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    abline(v=990:1000,col='grey',lty=2)
  })
  
}

shinyApp(ui=ui, server = server)

