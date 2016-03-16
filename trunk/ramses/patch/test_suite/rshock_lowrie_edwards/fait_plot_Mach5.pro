pro fait_plot,nout=nout

 rd_1d_Temp,file='res_Mach5.log',dat;,ngroup=1
 out=*dat[nout]          

readcol,'Tg_Mach5.txt',xtg,tgi
readcol,'Tr_Mach5.txt',xtr,tri

tp=out.p/out.d*1.66d-24/1.38d-16
tr=(out.er/7.5657d-15)^0.25
u=out.u
aa=max(abs(tp),xa)

print,out.x(xa-1)
out.x=out.x-out.x(xa-1)
set_plot,'x'
window,0
plot,out.x,tp,psym=6,xr=[-2500,500]
oplot,out.x,tr,psym=4

oplot,out.x,out.l*10,line=1

a=sort(xtg)
xtg=xtg(a)*1.d5
tgi=tgi(a)
coef=8.57274781455567/855.7;tgi(0)/tp(0)
tgi=tgi/coef

;; coef=tri(0)/tr(0)
;; tri=tri/coef

;; coef=tri(0)/tr(0)

print,tr

;calcul error

sz=size(xtr)
print,sz
nsz=sz(1)

for i=0,out.nc-1 do begin
   for j=1,nsz-1 do begin
      if(xtr(j) gt out.x(i)) then begin
         aa=(tri(j)-tri(j-1))/(xtr(j)-xtr(j-1))
         b=tri(j)-aa*xtr(j)
         er_int = aa*out.x(i)+b
         out.p(i)=abs(tr(i)-er_int)/er_int
         break
      endif
   endfor
endfor

sz=size(xtg)
print,sz
nsz=sz(1)
print,xtg,tgi
for i=0,out.nc-1 do begin
   for j=1,nsz-1 do begin
      if(xtg((j)) gt out.x(i)) then begin
         aa=(tgi((j))-tgi((j-1)))/(xtg((j))-xtg((j-1)))
         b=tgi((j))-aa*xtg((j))
         er_int = aa*out.x(i)+b
         out.u(i)=abs(tp(i)-er_int)/er_int
         break
      endif
   endfor
endfor
;; print,a
;; for i=0,out.nc-1 do begin
;;    for j=1,nsz-1 do begin
;;       if(xtg(a(j)) gt out.x(i)) then begin
;;          aa=(tgi(a(j))-tgi(a(j-1)))/(xtg(a(j))-xtg(a(j-1)))
;;          b=tgi(a(j))-aa*xtg(a(j))
;;          er_int = aa*out.x(i)+b
;;          out.u(i)=abs(tp(i)-er_int)/er_int
;; print,out.u(i),xtg(a(j)),tgi(a(j))
;;          break
;;       endif
;;    endfor
;; endfor






window,1

plot,out.x,tp,psym=6,xr=[-50,50],yr=[840,920]
oplot,out.x,tr,psym=4

oplot,out.x,out.l*10,line=1


oplot,xtr,tri,color=2
;oplot,xtg(a)*1d5,tgi(a),color=4

;oplot,xtr,tpi,color=3




plot,out.x,tp,psym=6,xr=[-3000,500],yr=[75,975],/xstyle,ystyle=9,ytitle='Temperature (K)',xtitle='Distance (cm)'
axis,yaxis=1,ytitle='!17 AMR level',yr=[0,7],/ynozero
oplot,out.x,tr,psym=4
oplot,xtr,tri,color=2
oplot,xtg,tgi,color=4

oplot,out.x,(out.l-5)*975./7.+75,line=1


set_plot,'ps'
thick=3.2 ; thickness of lines and characters
csize=1.7
set_plot,'ps'
device,filename='supercritical.eps',/encaps,yoffset=5.3,ysize=17,xsize=20,/color
!p.charsize=csize
!p.charthick=thick
!p.thick=thick
!x.thick=thick
!y.thick=thick
!p.position=[0.15,0.4,0.9,0.95]
iplot=0


plot,out.x,tp,psym=6,xr=[-3000,400],yr=[75,975],/xstyle,ystyle=9,ytitle='Temperature (K)',/noerase,xtickn=[' ',' ',' ',' ',' ',' ',' ',' ',' '];,xtitle='Distance (cm)',xtitle='Distance (cm)'
axis,yaxis=1,ytitle='!17 AMR level',yr=[0,8],/ynozero,/ystyle
oplot,out.x,tr,psym=4
oplot,xtr,tri,color=2
oplot,xtg,tgi,color=4

;oplot,out.x,(out.l-4)*90+75,line=1
oplot,out.x,(out.l-5)*900./8.+75,line=1

t0=textoidl("AMR level")
t1=textoidl("T_{gas} ana")
t2=textoidl("T_{gas}")
t3=textoidl("T_{rad} ana")
t4=textoidl("T_{rad}")

legend,[t0,t1,t2,t3,t4],/fill,psym=[0,0,6,0,4],color=[0,4,0,2,0],box=0,spacing=2., /top, /left, thick=[2.5,2.5,2.5,2.5,2.5],linsize=0.4,lines=[1,0,0,0,0]
;xyouts,150,100,'!17(b)'

!p.position=[0.15,0.1,0.9,0.4]
iplot=1
plot_io,out.x,out.p,xr=[-3000,400],ytitle='Relative error',/noerase,/xstyle,yr=[1.d-6,1],xtitle='Distance (cm)',ytickn=['10!u-6',' ','10!u-4',' ','10!u-2',' ',' '];,ytickn=[' ','10!u-6','10!u-4','10!u-2',' ']
oplot,out.x,out.p,color=2,line=0
oplot,out.x,out.u,color=4,line=0
;legend,[t2,t3,t4],/fill,lines=[0,0,0],color=[0,2,4],box=0,spacing=1.5, /bottom, /left, thick=[2.5,2.5,2.5],linsize=0.3


device,/close
print,out.p
;; print,di


;; print,pgasi
;; print,tpi

end
