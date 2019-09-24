dans sink_particle : 

l'ancienne correction :
if(msink(isink)-msink_all(isink) .gt. 50.*msink_all(isink)) xsink(isink,1:ndim)=xsink(isink,1:ndim)/msink(isink)
n'a pas été implémentée



la routine "get_height_scale" qui se trouvait dans courant_fine a été bougée dans cloud_module

appliquer les changement à hydro_boundary.f90
get_cell_index3 qui se trouvait dans hydro_boundary a été bougée dans cloud_module

note : 
#if USE_M_1==1
                 ! Reflection for the radiative flux                                                                                              
                 ! [E_nener(1),E_nener(2),E_nener(3),E(1),E(2),Fx(1),Fx(2),Fy(1),Fy(2),Fz(1),Fz(2)]                                               
                 !                                                                                                                                
                 ! ivar=8+nener-ngrp+igrp+idim*ngrp so ivar-9-nener+ngrp=igrp-1 + idim*ngrp                                                       
                 ! so igrp-1=ivar-9-nener+ngrp [ngrp]                                                                                             
                 if(ivar>8+nener.and.ivar<=8+nener+ndim*ngrp) then
                    igrp=1+modulo(ivar-9-nener+ngrp,ngrp)
                    idim=(ivar-8-nener+ngrp-igrp)/ngrp
                    switch=gs(idim)
                 endif
#endif

qui est dans le hydro_boundary the ramses_ism n'a pas été reporté 




copie de rt_boundana


remplacement d eget_cell_index3 par get_cell_index4 car get_cell_index3 est defini dans hydro/cooling_frig_module

 
vérifier : 

rt_hydro_boundary 




refaire les tests

le garder 

le donner à matthias pour la test suite

ordre 2


