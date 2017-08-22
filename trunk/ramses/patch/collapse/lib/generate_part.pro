;; cette routine permet de generer des particules traceurs pour ramses
pro generate_part
;boxlen=50.
Npart = 100L
dir=''

print, 'write the file ', dir+'ic_part'
openw,lun,/get_lun,dir+'ic_part'

for i=0L,Npart-1L do begin

for j=0L,Npart-1L do begin


for k=0L,Npart-1L do begin

rho = (0.001+ float(i) / float(Npart-1))*0.28d0

phi = (float(j) / float(Npart))*!pi ; 180.0d0

theta =(float(k) / float(Npart))*2.0d0*!pi;360.0d0


x = 0.5d0 + rho*sin(phi)*cos(theta);(0.25 +  float(i) / (Npart -1L)*0.5) 

y = 0.5d0 + rho*sin(phi)*sin(theta); (0.25 +  float(j) / (Npart -1L)*0.5) 

z = 0.5d0 + rho*cos(phi); (0.25 +  float(k) / (Npart -1L)*0.5) 

printf,lun,x,y,z,0.,0.,0.,0.,format='(7f13.7)'

;print,x,y,z,rho,theta,phi,cos(phi)
endfor
endfor
endfor

close,lun

end
;
;        filename=TRIM(initfile(levelmin))//'/ic_part'
;        if(myid==1)then
;           open(10,file=filename,form='formatted')
;           indglob=0
;        end if
;        eof=.false.
;
;        do while (.not.eof)
;           xx=0.0
;           if(myid==1)then
;              jpart=0
;              do i=1,nvector
;                 read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
;                 jpart=jpart+1
;                 indglob=indglob+1
;                 xx(i,1)=xx1
;                 xx(i,2)=xx2
;                 xx(i,3)=xx3
;                 vv(i,1)=vv1
;                 vv(i,2)=vv2
;                 vv(i,3)=vv3
;                 mm(i  )=mm1
;                 ii(i  )=indglob
;              end do




