pro fait_plot
  rd_1d_temp,file='log',out    
  dat=*out[1]

plot,dat.x,alog10(dat.er),psym=6

xx=findgen(100)/100.
chi=1.d10;/7.5d0
x0=0.5
E0=1.d5
tt=dat.t
ana2 = E0/(2.0d0*(chi*3.1415*tt)^.5)
ana2 = 1.+ana2*exp(-(xx-x0)^2./(4.0d0*tt*chi));^0.5)
print,tt
oplot,xx,alog10(ana2),color=2
end
