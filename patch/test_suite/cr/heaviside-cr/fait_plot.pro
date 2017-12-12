pro fait_plot
  
  rd_1d_temp,ngroup=1,file='ana.log',out    
  ana=out

  rd_1d_temp,ngroup=1,file='VN_4.log',out    
  vn=out

  rd_1d_temp,ngroup=1,file='DIRICHLET_4.log',out    
  dir=out

  rd_1d_temp,ngroup=1,file='ROBIN05_4.log',out    
  rob=out


nout=81

er_tot = fltarr(nout,3)
time   =  fltarr(nout,3)
for iout=0,nout-1 do begin
   out_vn = *vn[iout]
   out_dir=*dir[iout]
   out_rob=*rob[iout]

   er_tot(iout,0)=total(out_vn.er*0.5^out_vn.l)
   er_tot(iout,1)=total(out_dir.er*0.5^out_dir.l)
   er_tot(iout,2)=total(out_rob.er*0.5^out_rob.l)
  
   time(iout,0)=out_vn.t
   time(iout,1)=out_dir.t
   time(iout,2)=out_rob.t

endfor

nout=15

er_ana = fltarr(nout)
time_ana   =  fltarr(nout)

for iout=0,nout-1 do begin
   out_ana=*ana[iout]

   er_ana(iout)=total(out_ana.er*0.5^out_ana.l)
   time_ana(iout)=out_ana.t
endfor

set_plot,'x'
plot,time_ana,er_ana/er_ana,xr=[1.d-13,2.05d-12],yr=[0.95,1.05],line=1
oplot,time(*,0),er_tot(*,0)/1.d5,color=4
oplot,time(*,1),er_tot(*,1)/1.d5,color=1
oplot,time(*,2),er_tot(*,2)/1.d5,color=2




set_plot,'ps'
thick=3.2 ; thickness of lines and characters
csize=1.7
set_plot,'ps'
device,filename='cons_energy.eps',/encaps,yoffset=5.3,ysize=17,xsize=20,/color
!p.charsize=csize
!p.charthick=thick
!p.thick=thick
!x.thick=thick
!y.thick=thick
!p.position=[0.15,0.15,0.9,0.95]
iplot=0

plot,time_ana,er_ana/er_ana,xr=[1.d-13,2.0d-12],yr=[0.96,1.1],line=1,ytitle=textoidl("\SigmaE_r(t)/E_0"),xtitle='time (s)'
oplot,time(*,0),er_tot(*,0)/1.d5,color=4
oplot,time(*,1),er_tot(*,1)/1.d5,color=0
oplot,time(*,2),er_tot(*,2)/1.d5,color=2

xyouts,1.d-13,0.965,'!17(b)'

t1=textoidl("Analytic")
t2=textoidl("\alpha=1")
t3=textoidl("\alpha=0.5")
t4=textoidl("\alpha=0")
t5=textoidl("AMR level")

legend,[t2,t3,t4],/fill,lines=[0,0,0],color=[0,2,4],box=0,spacing=1., /bottom, /right, thick=[2.5,2.5,2.5],linsize=0.3

device,/close

device,filename='dirac.eps',/encaps,yoffset=5.3,ysize=17,xsize=20,/color
out_ana=*ana[2]
out_dir=*dir[8]
out_rob=*rob[8]
out_vn=*vn[8]

xx=findgen(100)/100.
chi=1.d10;/7.5d0
x0=0.5
E0=1.d5
tt=out_ana.t
anac = E0/(2.0d0*(chi*3.1415*tt)^.5)
ana2 = 1.+anac*exp(-(xx-x0)^2./(4.0d0*tt*chi));^0.5)

out_dir.p=1.+anac*exp(-(out_dir.x-x0)^2./(4.0d0*tt*chi));^0.5)
out_rob.p=1.+anac*exp(-(out_rob.x-x0)^2./(4.0d0*tt*chi));^0.5)
out_vn.p=1.+anac*exp(-(out_vn.x-x0)^2./(4.0d0*tt*chi));^0.5)

out_dir.p=abs(out_dir.p-out_dir.er)/out_dir.p
out_rob.p=abs(out_rob.p-out_rob.er)/out_rob.p
out_vn.p=abs(out_vn.p-out_vn.er)/out_vn.p


;oplot,xx,alog10(ana2),psym=6

;; plot,out_ana.x,alog10(out_ana.er),line=1,ytitle=textoidl("log(E_r)"),xtitle='x',/ynozero,ystyle=8

!p.position=[0.15,0.4,0.9,0.95]
iplot=0

plot,xx,alog10(ana2),line=2,ytitle=textoidl("log(E_r)"),/ynozero,ystyle=8,yr=[0,6],/noerase,xtickn=[' ',' ',' ',' ',' ',' ']

oplot,out_dir.x,alog10(out_dir.er),color=0
oplot,out_rob.x,alog10(out_rob.er),color=2
oplot,out_vn.x,alog10(out_vn.er),color=4

axis,yaxis=1,ytitle='!17 AMR level',yr=[4,9]
oplot,out_dir.x,1.2*(out_dir.l-4),line=1

xyouts,0.05,0.3,'!17(a)'

t1=textoidl("Analytic")
t2=textoidl("\alpha=1")
t3=textoidl("\alpha=0.5")
t4=textoidl("\alpha=0")
t5=textoidl("AMR level")
legend,[t5,t1,t2,t3,t4],/fill,lines=[1,2,0,0,0],color=[0,0,0,2,4],box=0,spacing=1., /top, /right, thick=[2.5,2.5,2.5,2.5,2.5],linsize=0.3

!p.position=[0.15,0.1,0.9,0.4]
iplot=1
plot_io,out_dir.x,out_dir.p,xtitle='x',ytitle='Relative error',/noerase,yr=[1.d-5,10.],ytickn=['10!u-5',' ','10!u-3',' ','10!u-1',' ','10'];,ytickn=[' ','10!u-6','10!u-4','10!u-2',' ']
oplot,out_rob.x,out_rob.p,color=2,line=0
oplot,out_vn.x,out_vn.p,color=4

device,/close
print,out_dir.l

print,time_ana,er_ana
print,time(*,0),er_tot(*,0)
print,time(*,1),er_tot(*,1)
print,time(*,2),er_tot(*,2)
end
