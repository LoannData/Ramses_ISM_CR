#!/usr/import/epd/python
# -*- coding: utf-8 -*-
#
# Script from M. Joos, adapted by M. Gonzalez
#
# run map2d.py -p . -o 89 -s 0.01,0.01 -L -k
import sys
import getopt
import pymses as pm
import numpy  as np
import pylab  as pl
import math   as ma
import os
import glob
from pymses.utils import constants as ct
from pymses.filters import CellsToPoints
from pymses.analysis.visualization import Camera, SliceMap, ScalarOperator, FractionOperator, raytracing
#From pyfrag.utils.files import *

bold  = "\033[1m"
reset = "\033[0;0m"

def usage():
    print(bold + "Usage: %%s [<options>]" + reset)
    print(bold + "-h --help:      " + reset + "print usage summary")
    print(bold + "-p --path:      " + reset + "path of directory containing the outputs")
    print(bold + "-o --output:    " + reset + "number of the output to plot")
    print(bold + "-k --sink:      " + reset + "if 'True', plot sinks position")
    print(bold + "-c --center:    " + reset + "center of the slice in [0,1]")
    print(bold + "-s --size:      " + reset + "size of the slice in [0,1]")
    print(bold + "-m --minimum:   " + reset + "minimum value for the colormap")
    print(bold + "-M --maximum:   " + reset + "maximum value for the colormap")
    print(bold + "-t --typlot:    " + reset + "variable to plot")
    print(bold + "-r --record:    " + reset + "save the plot")
    print(bold + "-i --igrp:      " + reset + "group number of Er to plot, if -1 Ertot")
    print(bold + "-l --los:       " + reset + "line of sight")

def main(argv):
    try:
        opts,args = getopt.getopt(argv, "hp:d:o:kc:s:m:M:t:ri:",
           ["help", "path=", "output=", "sink", "center=" \
            , "size=", "minimum=", "maximum=", "typlot=", "record", "igrp="])
    except getopt.GetoptError:
        usage()
        sys.exit(0)

    path      = "."
    output    = 1
    sink      = False
    center    = [0.5,0.5,0.5]
    size      = [0.5,0.5]
    minimum   = None
    maximum   = None
    typlot    = "rho"
    record    = False
    igrp      = -1
    los       = "x"

    for o, a in opts:
        if o in ["-h", "help"]:
            usage()
            sys.exit(0)
        elif o in ["-p", "path"]:
            path = a
        elif o in ["-o", "output"]:
            output = eval(a)
        elif o in ["-k", "sink"]:
            sink = True
        elif o in ["-c", "center"]:
            try:
                center = map(float, a.split(","))
                assert len(center) == 3
            except ValueError:
                print >>sys.stderr, "Invalid center coordinates"
                sys.exit(1)
            except AssertionError:
                print >>sys.stderr, "Invalid center coordinates"
                sys.exit(1)
        elif o in ["-s", "size"]:
            try:
                size = map(float, a.split(","))
                assert len(size) == 2
            except ValueError:
                print >>sys.stderr, "Invalid size"
                sys.exit(1)
            except AssertionError:
                print >>sys.stderr, "Invalid size"
                sys.exit(1)
        elif o in ["-m", "minimum"]:
            minimum = eval(a)
        elif o in ["-M", "maximum"]:
            maximum = eval(a)
        elif o in ["-t", "typlot"]:
            typlot = a
        elif o in ["-r", "record"]:
            record = True
        elif o in ["-i", "igrp"]:
            igrp = eval(a)
        elif o in ["-l", "los"]:
            los = a
    source = "".join(args)        

    myplot(path=path,output=output,sink=sink,center=center, \
               los=los,size=size,minimum=minimum,maximum=maximum, \
               typlot=typlot,record=record,igrp=igrp)


def myplot(path='.',output=-1,sink=False,center=[0.5,0.5,0.5],los='x',size=[0.5,0.5], \
           minimum=None,maximum=None,typlot='rho',record=False,igrp=-1,plot_velocity=False,slice=False):
    
    filerep       = path 
    if(record):
        savedir       = os.path.join(filerep,"results")
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        #mu,alpha,lmax = parseDirSink(directory)
        #savedir       = checkDir(savedir, "mu" + mu)
        #savedir       = checkDir(savedir, "lmax" + lmax)
            
    # if output=-1, take last output
    if output==-1:
        outputname = np.sort(glob.glob(os.path.join(filerep,"output_?????")))[-1]
        output=int(outputname.split('_')[-1])

    # Define the file
    ro = pm.RamsesOutput(filerep, output)

    if igrp != -1 and (igrp>ro.info["ngrp"] or igrp<0):
        print(bold+'PROBLEM:'+reset+' igrp='+str(igrp)+'  but simulation with ngrp='+str(ro.info["ngrp"]))
        sys.exit(1)

    # Define the output data structure (to use with pymses5)
    ro.define_amr_scalar_field("hydro", "rho", 0)
    ro.define_amr_vector_field("hydro", "vel", [1, 2, 3])
    ro.define_amr_vector_field("hydro", "Bl", [4, 5, 6])
    ro.define_amr_vector_field("hydro", "Br", [7, 8, 9])
    ro.define_amr_scalar_field("hydro", "P", 10)
    ro.define_amr_multivalued_field("hydro", "Er", 11, ro.info["ngrp"])
    ro.define_amr_vector_field("hydro", "J", [11+ro.info["ngrp"], 12+ro.info["ngrp"], 13+ro.info["ngrp"]])
    ro.define_amr_scalar_field("hydro", "eint", 14+ro.info["ngrp"])
    if ro.info["eos"]:
        ro.define_amr_scalar_field("hydro", "T", 15+ro.info["ngrp"])
    ro.define_amr_vector_field("grav", "g", [0, 1, 2])

    time   = ro.info["time"]*ro.info["unit_time"].express(ct.kyr)
    if(sink):
        # Read sink particles
        sink_filename=os.path.join(filerep,'output_%05d' %output +'/sink_%05d' %output +'.csv')
        print("Reading sinks     : " + reset + sink_filename)
        sinkp=np.loadtxt(sink_filename,dtype={'names': ('Id','mass','x','y','z','vx','vy','vz','rot_period','lx','ly','lz','acc_rate','acc_lum','age','int_lum','Teff'),'formats':('i','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f')},delimiter=',')
        sinkp=np.atleast_1d(sinkp) # to work even if only 1 sink
        nsink=sinkp.size

        # Define the accretion radius
        if("ncell_racc" in ro.info):
            ncell_racc=ro.info['ncell_racc']
        else:
            ncell_racc=4

        dxmin = 0.5**ro.info['levelmax']*ro.info["unit_length"].express(ct.au)
        r_acc = dxmin*ncell_racc
        #area  = np.pi*r_acc**2

        xc,yc,zc = sinkp['x']/ro.info['boxlen'],sinkp['y']/ro.info['boxlen'],sinkp['z']/ro.info['boxlen']
    else:
        nsink=0

    if nsink==0:
        nsink=1
        # Center
        xc,yc,zc = np.atleast_1d(center[0]),np.atleast_1d(center[1]),np.atleast_1d(center[2])
        
    cmap = pl.cm.get_cmap('seismic',100) 
    if typlot == 'rho':
        cmap = pl.cm.get_cmap('jet',100) 
    elif typlot == 'Tr':
        cmap = pl.cm.get_cmap('hot',100)
    elif typlot == 'B':
        cmap = pl.cm.get_cmap('PuOr',100)
    elif typlot == 'beta_plasma':
        cmap = pl.cm.get_cmap('seismic',100)
    else:
        cmap = pl.cm.get_cmap('jet',100) 


    # Width of the region to plot
    xr,yr = size[0],size[1]
    xbl =  (- xr/2.*ro.info["unit_length"]).express(ct.au)
    xbr =  (+ xr/2.*ro.info["unit_length"]).express(ct.au)
    ybl =  (- yr/2.*ro.info["unit_length"]).express(ct.au)
    ybr =  (+ yr/2.*ro.info["unit_length"]).express(ct.au)
    extent = [xbl,xbr,ybl,ybr]

    # Define camera
    if los == 'x':
        upvec = "z"
        vx = 1
        vy = 2
        labelx = 'Y (AU)'
        labely = 'Z (AU)'
        plane = 'yz'
    elif los == 'y':
        upvec = "x"
        vx = 2
        vy = 0
        labelx = 'Z (AU)'
        labely = 'X (AU)'
        plane = 'zx'
    elif los == 'z':
        upvec = "y"
        vx = 0
        vy = 1
        labelx = 'X (AU)'
        labely = 'Y (AU)'        
        plane = 'xy'

    def plot_func(dset):
        if typlot == 'rho':
            return dset['rho']*ro.info["unit_density"].express(ct.g/ct.cm**3)
        if typlot == 'eint':
            return dset['eint']*(ro.info["unit_density"]*ro.info["unit_velocity"]**2).express(ct.erg/ct.cm**3)
        if typlot == 'T':
            if ro.info["eos"]:
                return dset['T']
            else:
                return dset['P']/dset['rho']*ro.info["unit_temperature"].express(ct.K)*ro.info["mu_gas"]
        if typlot == 'B':
            B = 1./4*( (dset['Bl'][:,0]+dset['Br'][:,0])**2  \
                      +(dset['Bl'][:,1]+dset['Br'][:,1])**2  \
                      +(dset['Bl'][:,2]+dset['Br'][:,2])**2)
            return B
        # To use with pymses5
        if typlot == 'Er' or typlot == 'Tr':
            if igrp==-1:
                Er=np.sum(dset['Er'],axis=-1)*(ro.info["unit_density"]*ro.info["unit_velocity"]**2).express(ct.erg/ct.cm**3)
            else:
                Er=dset['Er'][:,igrp-1]*(ro.info["unit_density"]*ro.info["unit_velocity"]**2).express(ct.erg/ct.cm**3)

            if typlot=='Er':
                return Er
            else:
                return (Er/ct.a_R.express(ct.erg/ct.cm**3/ct.K**4))**0.25
                
        if typlot == 'entropy':
            return dset['P']/dset['rho']**ro.hydro_info["gamma"]

        if typlot == 'Pth_Pdyn':
            return dset['P']/(0.5*dset['rho']*np.sum(dset['vel']**2,axis=-1))

        if typlot == 'beta_plasma':
            B = 1./4.*( (dset['Bl'][:,0]+dset['Br'][:,0])**2  \
                      +(dset['Bl'][:,1]+dset['Br'][:,1])**2  \
                      +(dset['Bl'][:,2]+dset['Br'][:,2])**2)
            if igrp==-1:
                Er=np.sum(dset['Er'],axis=-1)#*(ro.info["unit_density"]*ro.info["unit_velocity"]**2).express(ct.erg/ct.cm**3)
            else:
                Er=dset['Er'][:,igrp-1]#*(ro.info["unit_density"]*ro.info["unit_velocity"]**2).express(ct.erg/ct.cm**3)

            return dset['P']/(B)#/dset['rho'])

        else:
            print(bold + 'Problem in typlot')
            sys.exit()
            
    # Fields to be read
    fields_to_read=[]
    if (not ro.info["eos"] and typlot=='T') or typlot=='entropy':
        fields_to_read=["rho","P"]
    elif typlot == 'Tr':
        fields_to_read=["Er"]
    elif typlot == 'Pth_Pdyn':
        fields_to_read=["P","rho","vel"]
    elif typlot == 'B':
        fields_to_read=["Bl","Br"]
    elif typlot == 'beta_plasma':
        fields_to_read=["Bl","Br","rho","P","Er"]
    else:
        fields_to_read=[typlot]
    if(plot_velocity):
        fields_to_read.append("vel")
    if not slice:
        fields_to_read.append("rho")
        
    fields_to_read=list(set(fields_to_read)) #to remove duplicates

    if slice:
        source = ro.amr_source(fields_to_read)
        varplot_op = ScalarOperator(plot_func)
    else:
        rt = raytracing.RayTracer(ro,fields_to_read)

        if typlot == 'rho':
            func = lambda dset: plot_func(dset)*ro.info["unit_length"].express(ct.cm)
            varplot_op = ScalarOperator(func)
        else:
            up_func = lambda dset: (dset["rho"]*ro.info["unit_density"].express(ct.g/ct.cm**3)*plot_func(dset))
            down_func = lambda dset: (dset["rho"]*ro.info["unit_density"].express(ct.g/ct.cm**3))
            varplot_op = FractionOperator(up_func,down_func)
        
    for i in range(nsink):
        fig = pl.figure()
        ax  = pl.subplot(111)

        pl.xlabel(labelx)
        pl.ylabel(labely)

        cam = Camera(center=[xc[i], yc[i], zc[i]],line_of_sight_axis=los,region_size=[xr, yr] \
                         , up_vector=upvec,map_max_size=512, log_sensitive=True)

        if slice:
            mapx = SliceMap(source, cam, varplot_op, z=0.)
        else:
            # WARNING! surf_qty=True to get an integral over dz and not dSdz
            # so that if typlot=='rho', get surface density map
            # and if not, get a density-weighted map (and not a mass-weighted map)
            mapx = rt.process(varplot_op,cam,surf_qty=True)

        if typlot=='Pth_Pdyn':
            imx    = pl.imshow(mapx.transpose(), extent = extent, origin='lower' \
                                   , vmin = minimum, vmax = maximum,cmap=cmap)
            #imx    = pl.contourf(mapx.transpose(), extent = extent, origin='lower' \
            #                       , vmin = minimum, vmax = maximum)
            cbar   = pl.colorbar()
            cbar.set_label(typlot)
        else:
            imx    = pl.imshow(np.log10(mapx.transpose()), extent = extent, origin='lower' \
                                   , vmin = minimum, vmax = maximum,cmap=cmap)
            #imx    = pl.contourf(np.log10(mapx.transpose()), extent = extent, origin='lower' \
            #                       , vmin = minimum, vmax = maximum)
            cbar   = pl.colorbar()

#            cmap = pl.cm.get_cmap('jet',100)

            if typlot=='rho':
               cbar.set_label( r'$\mathrm{log(}\rho)$', fontsize=15)
            elif typlot == 'Tr':
                cbar.set_label( r'$\mathrm{log(T_r)}$', fontsize=15)
            elif typlot == 'Pth_Pdyn':
                cbar.set_label('Pth/Pdyn')
            elif typlot == 'B':
                cbar.set_label('log(|B|)')
            elif typlot == 'beta_plasma':
                cbar.set_label( r'$\mathrm{log(}\beta)$', fontsize=15)
            else:
                cbar.set_label('Log('+typlot+')')


        if(sink):
            for isink in range(nsink):
                # Plot sinks position
                mass=sinkp[isink]["mass"]
                col  = cmap(int(mass*100))
                xp = (sinkp[isink]['x']/ro.info['boxlen']-xc[i])*ro.info["unit_length"].express(ct.au)
                yp = (sinkp[isink]['y']/ro.info['boxlen']-yc[i])*ro.info["unit_length"].express(ct.au)
                zp = (sinkp[isink]['z']/ro.info['boxlen']-zc[i])*ro.info["unit_length"].express(ct.au)

                if los == 'x':
                    ax.add_patch(pl.Circle((yp,zp),radius=r_acc,fc='white', alpha=0.65))
                    ax.add_patch(pl.Circle((yp,zp),radius=(r_acc/ncell_racc),fc=col,ec=col))
                #pl.plot(yp,zp,'k+',mew=2,ms=10)
                elif los == 'y':
                    ax.add_patch(pl.Circle((zp,xp),radius=r_acc,fc='white', alpha=0.65))
                    ax.add_patch(pl.Circle((zp,xp),radius=(r_acc/ncell_racc),fc=col,ec=col))
                #pl.plot(zp,xp,'k+',mew=2,ms=10)
                elif los == 'z':
                    ax.add_patch(pl.Circle((xp,yp),radius=r_acc,fc='white', alpha=0.65))
                    ax.add_patch(pl.Circle((xp,yp),radius=(r_acc/ncell_racc),fc=col,ec=col))
                #pl.plot(xp,yp,'k+',mew=2,ms=10)    

        # Plot velocity field
        if(slice and plot_velocity):
            p = cam.get_slice_points(0.0)
            nx, ny = cam.get_map_size()
            dset = pm.analysis.sample_points(source, p)
            vel  = dset["vel"]*ro.info["unit_velocity"].express(ct.km/ct.s)
            rs = 32
            x = np.linspace(xbl,xbr,nx)
            y = np.linspace(ybl,ybr,ny)
            u,v = np.zeros((nx,ny)),np.zeros((nx,ny))
            mask = np.zeros((nx,ny))
            for ii in range(nx):
                for jj in range(ny):
                    if(ii%rs == 0 and jj%rs == 0):
                        u[ii,jj] = vel[:,vx].reshape(nx,ny)[ii,jj]
                        v[ii,jj] = vel[:,vy].reshape(nx,ny)[ii,jj]
                    else:
                        u[ii,jj] = 'Inf'
                        v[ii,jj] = 'Inf'
                        mask[ii,jj] = 1

            u2=u[::rs,::rs]
            v2=v[::rs,::rs]
            x2=x[::rs]
            y2=y[::rs]

            vel_mean=np.mean(np.sqrt(u2**2+v2**2))
            vel_max =np.max (np.sqrt(u2**2+v2**2))

            #pl.quiver(x,y,u.transpose(),v.transpose(),scale=20*nx)
            #Q=pl.quiver(x,y,u.transpose(),v.transpose(),scale=100,pivot='mid')
            #u_masked=np.ma.masked_array(u,mask=mask)
            #v_masked=np.ma.masked_array(v,mask=mask)
            #vel_mean=np.mean(np.sqrt(u_masked**2+v_masked**2))

            Q=pl.quiver(x2,y2,u2.transpose(),v2.transpose(),scale=100,pivot='mid')
            #pl.quiverkey(Q,0.7,0.92,vel_mean,r'%.2f km/s'%vel_mean,coordinates='figure')
            pl.quiverkey(Q,0.7,0.92,vel_max,r'%.2f km/s'%vel_max,coordinates='figure')
            pl.axis(extent)
            del(u,v,x,y,u2,v2,x2,y2)

            if(sink):
                # Print the mass of the sinks
                mass=0
                for j in range(nsink):
                    mass = mass + sinkp["mass"][j]
                    
                pl.figtext(0.01,00.01, "$\mathrm{M}_*=$"+ "%4.1f " %mass + "$M_\odot$", fontsize=15)

        pl.suptitle("%.3f kyr" %time)


        if(record):
            for format in ("png","eps","pdf"):
                if slice:
                    namefig = os.path.join(savedir,"slice_" + typlot + "_" + plane + "_eps" +"_%.3f" %size[0]+  "_out"+"_%05d" %output + "_isink_%05d." %sinkp['Id'][i] + format)
                else:
                    namefig = os.path.join(savedir,"proj_" + typlot + "_" + plane + "_eps" +"_%.3f" %size[0]+  "_out"+ "_%05d" %output + "_isink_%05d." %sinkp['Id'][i] + format)
                print(bold + "Save figure: " + reset + namefig)
                pl.savefig(namefig)
            pl.close()
        else:
            pl.ion()
            pl.show()

    if(sink):
        # Print the mass of the sinks
        for i in range(nsink):
            mass = sinkp["mass"][i]
            print(bold + "Mass(" + str(i) + "): " + reset + "%.2e Msun" %mass)
 
if __name__ == "__main__":
    main(sys.argv[1:])

def allplot(path='./',output=1,size=[0.02,0.02]):
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='rho',los='x',slice=True,plot_velocity=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='rho',los='y',slice=True,plot_velocity=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='rho',los='z',slice=True,plot_velocity=True)

    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='rho',los='x')
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='rho',los='y')
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='rho',los='z')

    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='T',los='x')
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='T',los='y')
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='T',los='z')

    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='T',los='x',slice=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='T',los='y',slice=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='T',los='z',slice=True)

    #myplot(path=path,output=output,size=size,record=True,sink=True,typlot='Tr',los='x')
    #myplot(path=path,output=output,size=size,record=True,sink=True,typlot='Tr',los='y')
    #myplot(path=path,output=output,size=size,record=True,sink=True,typlot='Tr',los='z')

    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='Tr',los='x',slice=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='Tr',los='y',slice=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='Tr',los='z',slice=True)

    #myplot(path=path,output=output,size=size,record=True,sink=True,typlot='entropy',los='x')
    #myplot(path=path,output=output,size=size,record=True,sink=True,typlot='entropy',los='y')
    #myplot(path=path,output=output,size=size,record=True,sink=True,typlot='entropy',los='z')

    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='entropy',los='x',slice=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='entropy',los='y',slice=True)
    myplot(path=path,output=output,size=size,record=True,sink=True,typlot='entropy',los='z',slice=True)

# map2d.allplot(path='../Greylmax13beta01facc0/',output=1995)
# map2d.allplot(path='../Greylmax14beta01facc0/',output=1035)
# map2d.allplot(path='../SPlmax13beta01facc0/',output=578)
# map2d.allplot(path='../SPlmax14beta01facc0/',output=313)

# map2d.allplot(path='../Greylmax13beta01facc0/',output=300)
# map2d.allplot(path='../Greylmax14beta01facc0/',output=300)
# map2d.allplot(path='../SPlmax13beta01facc0/',output=300)
# map2d.allplot(path='../SPlmax14beta01facc0/',output=300)
