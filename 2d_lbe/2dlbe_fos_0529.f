c******************************************************
c  use equilibrium  condition at the open boundary
c  nest_incompres4.f
c  simple extrapolation in at i=nx
c  symmetry condition at j=1 and ny
c  multi-relaxation-time model in the nested loop
c  nest_loop.f
c  from Like Li 12/2013
c******************************************************
	     implicit double precision (a-h,o-z) !modified
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension af(0:8),un(nxd,nyd),n_march(6),sum_flux(nxd),
     &             sig_xx(nxd),sig_xy(nxd),sig_yy(nxd),
     &             tau_xx(nxd),tau_xy(nxd),tau_yy(nxd),
     &	           n_used(6)
	     common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
	     common /solution/ux(nxd,nyd),uy(nxd,nyd),rho(nxd,nyd),
     &			          fi(0:8,nxd,nyd)
         common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
         common /con/ s2,s3,s5,s7,s8,s9
	     common /bk_c/nx_c,ny_c,imid_c,jmid_c,tau_c,radius_c
         common /cxy/cix,ciy,cf(0:8)
         common /maxerr/err_max,i_max,j_max
         common /rev/kreverse
	     common /constants/c0,c1,c2,c3,c4,c5,c6,c9,c12,c36
	     INTEGER  kreverse(8)/5,6,7,8,1,2,3,4/
	     real*8 cix(0:8)/0.,1.,1.,0.,-1.,-1.,-1.,0.,1.0/
	     real*8 ciy(0:8)/0.,0.,1.,1.,1.,0.,-1.,-1.,-1.0/
		 
         open(9,file='datain.lbe',status='unknown',action='read')
         open(20,file='velocity.his',status='unknown')
         open(30,file='fi.his',status='unknown')
         open(24,file='wall_nodes.dat',status='unknown')		 
         open(25,file='convergence.dat',status='unknown')
         open(38,file='force.dat',status='unknown')		 
         open(12,file='generaloutput.dat',status='unknown')
         open(14,file='iw_velocity.dat',status='unknown')
         open(52,file='del_info.dat',status='unknown')
c         open(6,file='debuging.dat',status='unknown')  screen output
		 
	     write(12,*)' nest_incompres4.f'
	     write(12,*)' simple extrapolation in at i=nx'
	     write(12,*)' symmetry condition at j=1 and ny'
         c0=4.0/9.0
         c1=1.0/9.0
         c2=1.0/36.0
         c3=3.0
         c4=4.5
         c5=2.0/3.0        
	     c6=1.0/6.0
         c9=1.0/9.0
         c12=1.0/12.0
	     c36=1.0/36.0

	     cf(0)=4.0/9.0
	     do i=1,7,2
	       cf(i)=1.0/9.0
	       cf(i+1)=1.0/36.0
	     enddo
cc-------datain-------------------------------------------------
	     read(9,*) key,tau,Re,radius,nx,ny,imid,err_converge,nw
	     read(9,*) n_level
	     read(9,*) (n_march(k),k=1,n_level)
	     s2=1.63 
	     s3=1.14 
	     s5=1.92 
	     s7=1.92
	     read(9,*)s2,s3,s5,s7
	     s8=1.0/tau
	     s9=s8
cc--------------------------------------------------------------		 
	     do k=1,n_level
	        write(12,*) '  n_march(',k,') =',n_march(k)
	     enddo
	     nt_last=0
	     call initial_cond(key,nt_last)
         visc=(2.*tau-1.)/6.
	     jrd=jmid+radius+2
	     iw=imid+6*radius
	     write(14,*) 'iw=',iw
	     write(38,*)'variables ="n", "cdx-mrt", "cdy-mrt"'
	     write(38,*)'zone T="drag_mrt"'
	     write(12,*)'variables ="y", "ux-mrt", "del-rho-mrt"'
	     write(12,*)'zone T="velocity-mrt"'
	     write(6,*)'finishing initialization'
         write(25,*) 'variables =ii,err2,i_max,j_max,err_max'
         do i=1,nx
         do j=1,ny
	         un(i,j)=ux(i,j)
	     enddo
   	     enddo
c**************************************************************
c**********   START THE MAJOR LOOP IN LBM SIMULATION  *********
c**************************************************************
   !!!!START OF ITERATION LEVEL SETUP
       DO 2000 nn=1,n_level
	      ntimes=n_march(nn)
	      nxm=nx-1
	      nym=ny-1
	   !!!START OF TIME ITERATION
          DO 1000 ii=nt_last+1,nt_last+ntimes
             ttao = tau 
           !!collision(equilibrate)
	         call collision(ttao)
           !!force calculation
	         call force_mom(fi,fxm,fym)
	         cdx_m=-fxm/(radius*vt**2.)
	         cdy_m=-fym/(radius*vt**2.)
	         write(38,388) ii,cdx_m,cdy_m
           !!streaming:
	         call streaming(fi)
           !!calculate physical variables: rho, ux, uy
	         call macro_variables
           !!checking convergence
	         call difference(un,err2)
             write(25,885) ii,err2,i_max,j_max,err_max
             write(6,885)  ii,err2,i_max,j_max,err_max
	         write(14,655) ii,ux(nxm,jmid)/vt,uy(nxm,jmid)/vt			 
388	         FORMAT(3x,i7,2x,f19.9,2x, e15.6)
885 	     FORMAT(3x,i8,2x,e14.6,2x,2i4,2x,e15.6)
655	         FORMAT(3x,i7,2f18.8)
           !!output I:
            !printing some typical velocity and derivative profiles:
            !note: the i=imid line goes through the center of the cylinder
	         km=mod(ii,1000)
	         kw=0
             if(km.eq.0.and.ii.ge.200) kw=1
             if(ii.eq.ntimes+nt_last)  kw=1
             if(err2.lt.err_converge)  kw=1
             if(kw.eq.1) then
	            i=imid
c	            write(12,*)'ii=',ii,'err2=',err2
	            do i=2,nx-1
	               sum_flux(i)=0
	               do j=2,ny-2
	                  sum_flux(i)=sum_flux(i)+ux(i,j)
	               enddo
                enddo
	            h=ny-3
	            factor=visc*vt/radius
	            rho_ref=rho(2,jmid)
	            call stress_y(imid,tau_xx,tau_xy,tau_yy)
	            write(12,*)'j   y/r   u/U   Cp   shear'
	            do j=2,ny-1
	               shear=tau_xy(j)/factor
	               cp=(rho(imid,j)-rho_ref)/(0.5*vt**2)
	               if(walls(imid,j).eq.1) then
		             shear=0
		             cp=2
	               endif
	               write(12,40) j,float(j-jmid)/radius,ux(imid,j)/vt,
     &                          cp,shear
	            enddo
	            call stress_x(jmid,sig_xx,sig_xy,sig_yy)
               !output central line Velocity, Cp, Shear and Sum_flux 
	            write(12,*)'x-velocity on the y=0 axis'
	            do i=2,nx-1
	               xx=float(i-imid)/radius
	               shear_norm=sig_xx(i)/factor
	               cp=(rho(i,jmid)-rho_ref)/(0.5*vt**2)
	               if(walls(i,jmid).eq.1) then
		              shear=0
		              cp=3
	               endif
	               write(12,44) i,xx,ux(i,jmid)/vt,cp,
     &	                        1-sum_flux(i)/(vt*h),shear_norm
	            enddo
	            relative=1-sum_flux(nx-1)/sum_flux(2)
	            write(12,*)'ii=',ii, ' cd_x=',cdx_m,'  cd_y=',cdy_m
	            write(12,*)'relative error in ux_flux=',relative
             endif
40	         FORMAT(2x,i3,2x,f10.6,2x,3f13.8)
44	         FORMAT(2x,i3,2x,f10.6,2x,4f13.8)
           !!output II:
            !save the flow field every 1000 time steps:
            !note: results for ux,uy, rho at i=1, i=nx, j=1 and j=ny lines are 
	        !not physical.  Only some of the fi's are needed at thoselocations; 
	         km=mod(ii,nw)
	         kw=0
	         if(km.eq.0.and.ii.ge.200) kw=1
	         if(ii.eq.ntimes) kw=1
	         if(err2.lt.err_converge) kw=1
	         if(kw.eq.1) then
	            rewind(20)
	            rewind(30)
	            write(20,*) tau,vt,Re,radius,nx,ny,imid,jmid,ii
                write(20,*) 'variables =x,y,ux,uy,rho'
                write(20,*) 'zone I=', nx,'J=',ny,'f= point'
	            do i=1,nx
	            do j=1,ny
	               write(20,88) i,j,ux(i,j),uy(i,j),rho(i,j)
	            enddo
	            enddo		   
                write(20,89) 1,imid,jmid
                write(20,*) radius		   	 
	            write(30,*)tau,vt,Re,radius,nx,ny,imid,jmid,ii
	            do i=1,nx
	            do j=1,ny
	               write(30,98) i,j,fi(0,i,j),fi(1,i,j),fi(2,i,j),
     &		                    fi(3,i,j),fi(4,i,j),fi(5,i,j),
     &                          fi(6,i,j),fi(7,i,j),fi(8,i,j)
	            enddo
	            enddo
             endif
88	            FORMAT(x,2i5,x,3f16.10)
89              FORMAT("GEOMETRY ZN=" I5.0 "X= "I15.7 "Y=" I15.7 
     &                 "FC=CUST2 LT=0.2 CS=GRID T=CIRCLE C=BLACK")
98	            FORMAT(x,2i5,x,9f16.10)
			  
	         if(err2.lt.err_converge) goto 1500
1000	  CONTINUE
       !!!END OF TIME ITERATION
 
1500	  ks=0
	      n_used(nn)=ii-1   
	      if(nn.eq.n_level ) goto 2000
	      write(12,*)'nn=',nn,' this level done'
	      write(12,*)'switch to next grid level'
	      write(25,*)'switch to next grid level'
	      write(38,*)'switch to next grid level'
	      nt_last = ii-1
	      nx_c=nx
	      ny_c=ny
	      imid_c=imid
	      jmid_c=jmid
	      radius_c=radius
	      tau_c=tau
	      call c_2_f_transfer(un)
	      visc=(2*tau-1)/6	
2000   CONTINUE 
   !!!!END OF ITERATION LEVEL SETUP

       do nn=1,n_level
	      k=n_level-nn
	      if(nn.eq.1) jt=n_used(1)
	      if(nn.gt.1) jt=n_used(nn)-n_used(nn-1) 
	      ww=jt/4.0**k 
	      work= work+ ww 
	      write(12,226)nn,n_used(nn),jt,ww,work
       enddo
226    FORMAT(2x,i3,2x,'n_used=',i6,2x,'time=',i6,2x,'wk=',f9.2,3x,
     &          'cumu_work=',f10.2)
	  
       STOP
       END
	   
c*****************************************************************	   
c***************************subroutine collision****************
c*****************************************************************
      	 subroutine collision(ttao)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     implicit double precision (a-h,o-z)
	     common /solution/vx(nxd,nyd),vy(nxd,nyd),rho(nxd,nyd),
     &                    fi(0:8,nxd,nyd)
	     common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
         common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /cxy/cix(0:8),ciy(0:8),cf(0:8)
         common /con/ s2,s3,s5,s7,s8,s9
	     common /rev/kreverse(8)
	     common /constants/c0,c1,c2,c3,c4,c5,c6,c9,c12,c36

         s8=1./ttao
	     s9=s8
         omega=1.0/ttao
	     ug=0
	     vg=0
		 
cc--Determine the boundary conditions for fi's on wb(i,j)=1
         k_del=0
	     do i=1,nx-1
	     do 100 j=1,ny
	       if(wb(i,j).eq.0) goto 100
	       do 120 k=1,8
		     ii=cix(k)
		     jj=ciy(k)
		     ia=i+ii
		     ja =j+jj
		     if(walls(ia,ja).eq.1) goto 120
800 		 k_del=k_del+1
		     del=dels(k_del)
		     ko=kreverse(k)

 		     if(del.ge.0.5) then
		        omegai=omega*(2.*del-1.0)
		     else
		        omegai=omega*(2.*del-1.0)/(1.0-2.*omega)
		     endif
            !set boundary node velocity at (i,j) where wb(i,j)=1
		     ratio=(del-1)/del
		     vxb=vx(ia,ja)*ratio+ug/del
		     vyb=vy(ia,ja)*ratio+vg/del
	  	     if(del.lt.0.5) then
		        i2a=i+2*ii
		        j2a=j+2*jj
		        vxb=vx(i2a,j2a)
		        vyb=vy(i2a,j2a)
		     endif
            !evaluate the equilibrium function on the solid nodes: wb(i,j)=1
	      	 vxx=vx(ia,ja)
	      	 vyy=vy(ia,ja)
			 
 	 	     vmag=1.5*(vxx*vxx + vyy*vyy)
		     rh=rho(ia,ja)

		     gf=cix(ko)*vxx+ciy(ko)*vyy
		     gb=cix(ko)*vxb+ciy(ko)*vyb
		     f_eqb=cf(k)*(rh-vmag+c3*gb+c4*gf*gf)
		  
	         fi_eq_ko_ia_ja=cf(ko)*(rh-vmag+c3*gf+c4*gf*gf)	   
            !setting the b.c.'s for fi's
		     fi(k,i,j)=((1-omega)*fi(ko,i+ii,j+jj)
     &		      + omega*fi_eq_ko_ia_ja)*(1-omegai)
     &		      + omegai*f_eqb
120	       continue
100	    continue
        enddo
		
cc--obtain the inlet conditions for fi's at i=1
        i=1
        do 200 j=2,ny-1
	       do 220 k=1,8
	         ii=cix(k)
	         jj=ciy(k)
		 
		     ia=i+ii
		     ja=j+jj
		     if(ii.le.0) goto 220
	         ko=kreverse(k)
            !set boundary node velocity at i=1
		     vxx=vx(ia,ja)
		     vyy=vy(ia,ja)

		     vmag=1.5*(vxx*vxx + vyy*vyy) 
		     rh=rho(ia,ja) 
		     gf=cix(ko)*vxx+ciy(ko)*vyy
            !evaluate the equilibrium function for the inlet nodes                  
	         fi_eq_ko_ia_ja=cf(ko)*(rh-vmag+c3*gf+c4*gf*gf)
		     vmag_ij=1.5*vt**2
		     gf=cix(k)*vt
	         fi_eq_k_i_j=cf(k)*(0-vmag_ij+c3*gf+c4*gf*gf)
	         fi_neq_ko_ia_ja=fi(ko,ia,ja)-fi_eq_ko_ia_ja
            !set up the inlet fi's
	         fi(k,1,j)= fi_eq_k_i_j
c     &           	  + fi_neq_ko_ia_ja
c     &  	          + 6*vt*cf(k)
220	       continue
200	    continue

cc--obtain the BC conditions for fi's at j=1 & j=ny
        j=1
        do 400 i=1,nx-1
	      do 420 k=1,8
	         ii=cix(k)
	         jj=ciy(k)
		 
		     ia=i+ii
		     ja=j+jj
		     if(jj.le.0) goto 420
	         ko=kreverse(k)
            !set boundary node velocity at j=1
		     vxx=vx(ia,ja)
		     vyy=vy(ia,ja)

		     vmag=1.5*(vxx*vxx + vyy*vyy) 
		     rh=rho(ia,ja) 
		     gf=cix(ko)*vxx+ciy(ko)*vyy
                 
	         fi_eq_ko_ia_ja=cf(ko)*(rh-vmag+c3*gf+c4*gf*gf)
		     vmag_ij=1.5*vt**2
		     gf=cix(k)*vt
	         fi_eq_k_i_j=cf(k)*(0-vmag_ij+c3*gf+c4*gf*gf)
	         fi_neq_ko_ia_ja=fi(ko,ia,ja)-fi_eq_ko_ia_ja

	         fi(k,i,j)= fi_eq_k_i_j
c     &           	  + fi_neq_ko_ia_ja
c     &  	          + 6*vt*cf(k)
420	      continue
400	    continue
        j=ny
        do 600 i=1,nx-1
	      do 620 k=1,8
	         ii=cix(k)
	         jj=ciy(k)
		 
		     ia=i+ii
		     ja=j+jj
		     if(jj.ge.0) goto 620
	         ko=kreverse(k)
            !set boundary node velocity at j=ny
		     vxx=vx(ia,ja)
		     vyy=vy(ia,ja)

		     vmag=1.5*(vxx*vxx + vyy*vyy) 
		     rh=rho(ia,ja) 
		     gf=cix(ko)*vxx+ciy(ko)*vyy
c                 
	         fi_eq_ko_ia_ja=cf(ko)*(rh-vmag+c3*gf+c4*gf*gf)
		     vmag_ij=1.5*vt**2
		     gf=cix(k)*vt
	         fi_eq_k_i_j=cf(k)*(0-vmag_ij+c3*gf+c4*gf*gf)
	         fi_neq_ko_ia_ja=fi(ko,ia,ja)-fi_eq_ko_ia_ja

	         fi(k,i,j)= fi_eq_k_i_j
c     &           	  + fi_neq_ko_ia_ja
c     &  	          + 6*vt*cf(k)
620	      continue
600	    continue

cc--carry out relaxation step for the interior points
        do i=2,nx-1
	       do 500 j=2,ny-1
	          if(walls(i,j).eq.1) goto 500
	          f0=fi(0,i,j)
	          f1=fi(1,i,j)
	          f2=fi(2,i,j)
	          f3=fi(3,i,j)
	          f4=fi(4,i,j)
	          f5=fi(5,i,j)
	          f6=fi(6,i,j)
	          f7=fi(7,i,j)
	          f8=fi(8,i,j)
 
	          sum_axis     = f1+f3+f5+f7
	          sum_diagnoal = f2+f4+f6+f8
	          rh=rho(i,j)
	          uu=vx(i,j)
	          vv=vy(i,j)
	          u_square=uu*uu+vv*vv

	          peq1=  -2*rh +3*u_square
	          peq2=     rh -3*u_square
              peq4=  -uu
              peq6=  -vv
	          peq7=  uu*uu-vv*vv
	          peq8=  uu*vv

	          f04 =  4*f0

      	      p1=-f04-  sum_axis+2*sum_diagnoal - peq1
	          p2= f04-2*sum_axis+  sum_diagnoal - peq2
	          p4= 2*(f5-f1)+f2-f4-f6+f8 	- peq4
	          p6= 2*(f7-f3)+f2+f4-f6-f8  	- peq6
	          p7= f1-f3+f5-f7  			- peq7
	          p8= f2-f4+f6-f8  			- peq8

	          s2p1=s2*p1
	          s3p2=s3*p2
	          s5p4=s5*p4
	          s7p6=s7*p6

	          a1b2=(s2p1+2*s3p2)*c36
	          a2b1=(2*s2p1+s3p2)*c36

	          t4p6=(s5p4+s7p6)*c12
	          t4m6=(s5p4-s7p6)*c12
	          s546=s5p4*c6
	          s766=s7p6*c6
	          t7=0.25*s8*p7
	          t8=0.25*s9*p8

	          fi(0,i,j)=f0 + (s2p1-s3p2)*c9
	          fi(1,i,j)=f1 + a1b2 + s546 - t7
	          fi(2,i,j)=f2 - a2b1 - t4p6 - t8 
	          fi(3,i,j)=f3 + a1b2 + s766 + t7
	          fi(4,i,j)=f4 - a2b1 + t4m6 + t8
	          fi(5,i,j)=f5 + a1b2 - s546 - t7
	          fi(6,i,j)=f6 - a2b1 + t4p6 - t8
	          fi(7,i,j)=f7 + a1b2 - s766 + t7
	          fi(8,i,j)=f8 - a2b1 - t4m6 + t8

500	       continue	
        enddo
		
cc--simple extrapolation at the downstrem at x=nx.
        fi(4,nx,1:ny) = fi(4,nx-1,1:ny)
        fi(5,nx,1:ny) = fi(5,nx-1,1:ny)
        fi(6,nx,1:ny) = fi(6,nx-1,1:ny)

      RETURN
      END

c******************************************************************
c************************subroutine force_mom**********************
c******************************************************************
	     subroutine force_mom(fi,fx,fy)
	     implicit double precision (a-h,o-z)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension fi(0:8,nxd,nyd)
	     common /wall/wb(nxd,nyd),w(nxd,nyd),dels(ndel),k_del_max
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /cxy/cix(0:8),ciy(0:8),cf(0:8)

c w(i,j)=walls(i,j)
c walls(i,j)=0.0 for fluid and it is 1.0 for solid part including the
c	     interior of the solid.
c wb(i,j)= 1 for boundary nodes and 0 for the rest of points.

c The force (fx,fy) is actually the force of the solid on the fluid.
c One needs to multiply -1 to get the force on the solid body. 

 	     fx=0
	     fy=0
	     do i=2,nx-1
	        ip=i+1
	        im=i-1
	        do j=2,ny-1
	           jp=j+1
	           jm=j-1
	           if(wb(i,j).eq.1) then
		       fx= fx +(fi(5,i,j)+fi(1,im,j ))*cix(5)*(1-w(im,j ))
     &			      +(fi(1,i,j)+fi(5,ip,j ))*cix(1)*(1-w(ip,j ))
     &		 	      +(fi(6,i,j)+fi(2,im,jm))*cix(6)*(1-w(im,jm))
     &		 	      +(fi(2,i,j)+fi(6,ip,jp))*cix(2)*(1-w(ip,jp))
     &		          +(fi(8,i,j)+fi(4,ip,jm))*cix(8)*(1-w(ip,jm))
     &		          +(fi(4,i,j)+fi(8,im,jp))*cix(4)*(1-w(im,jp))
		       fy= fy +(fi(7,i,j)+fi(3,i,jm ))*ciy(7)*(1-w(i ,jm))
     &		 	      +(fi(3,i,j)+fi(7,i,jp ))*ciy(3)*(1-w(i ,jp))
     &		 	      +(fi(6,i,j)+fi(2,im,jm))*ciy(6)*(1-w(im,jm))
     &		 	      +(fi(2,i,j)+fi(6,ip,jp))*ciy(2)*(1-w(ip,jp))
     &		          +(fi(8,i,j)+fi(4,ip,jm))*ciy(8)*(1-w(ip,jm))
     &		          +(fi(4,i,j)+fi(8,im,jp))*ciy(4)*(1-w(im,jp))
	           endif
	        enddo
	     enddo
		 
	     RETURN
	     END
	
c**********************************************************************
c**********************subroutine finding_dels*************************
c**********************************************************************
	     subroutine finding_dels
	     parameter(nxd=1001,nyd=401,ndel=600)
	     implicit double precision (a-h,o-z)
	     common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /cxy/cix(0:8),ciy(0:8),cf(0:8)
	     common /rev/kreverse(8)

	     sqroot2=sqrt(2.0)
	     rsq=radius**2
	     k_del = 0
	     xcnt=imid
	     ycnt=jmid
		 
c  Determine the boundary conditions for fi's on wb(i,j)=1
	     DO i=2,nx-1
	        DO 100 j=2,ny-1
	           if(wb(i,j).eq.0) goto 100
	           xgrid=i
	           ygrid=j
	           DO 120 k=1,8
		          ii=cix(k)
		          jj=ciy(k)
		          ia=i+ii
		          ja=j+jj
		          if(walls(ia,ja).eq.1) goto 120
		          if(k.eq.3.or.k.eq.7) then
		             x=i
		             y1=ycnt+sqrt(rsq-(x-xcnt)**2)
		             y2=ycnt-sqrt(rsq-(x-xcnt)**2)
		             dist=min(abs(y1-ygrid),abs(y2-ygrid))
		             del=1.0-dist
		             goto 800
		          endif
		          if(k.eq.1.or.k.eq.5) then
		            y=j
		            x1=xcnt+sqrt(rsq-(y-ycnt)**2)
		            x2=xcnt-sqrt(rsq-(y-ycnt)**2)
		            dist=min(abs(x1-xgrid),abs(x2-xgrid))
		            del=1.0-dist
		            goto 800
		          endif
		          if(k.eq.2.or.k.eq.6) ys=ycnt+xgrid-ygrid
		          if(k.eq.4.or.k.eq.8) ys=xgrid+ygrid-ycnt
		          qsqroot=sqrt(2*rsq+2*xcnt*ys-(xcnt**2+ys**2))
		          x1=0.5*(xcnt+ys+qsqroot)
		          x2=0.5*(xcnt+ys-qsqroot)
		          y1=x1-xgrid+ygrid
		          y2=x2-xgrid+ygrid
		          dist1=sqrt((x1-xgrid)**2+(y1-ygrid)**2)
		          dist2=sqrt((x2-xgrid)**2+(y2-ygrid)**2)
		          del=1-min(dist1,dist2)/sqroot2
800 		      k_del=k_del+1
                  if(k_del.gt.ndel) then
                       write(6,*)' ndel is too small for dels(ndel)'
                       write(6,*)' Computation is terminated.'
                       write(6,*)' Please modify the program((k_del)'
		          endif
		          dels(k_del)=del
	              write(52,648) i-imid,j-jmid,k,del
				  
120	           CONTINUE
100	        CONTINUE
	     ENDDO
		 
	     k_del_max=k_del
	     write(12,*)'k_del_max=',k_del_max
	     do i=1,k_del_max
	        write(12,688) i,dels(i)
	     enddo
648	     FORMAT(3x,'iimid=',i3,3x,'jjmid=',i3,3x,
     &          'k=',i3,3x,'del=',f12.6)
688	     FORMAT(3x,i4,3x,f12.6)
	     write(12,*)'k_del_max=',k_del_max
         write(6,*)'done finding_dels'

	     RETURN
	     END
		 
c*****************************************************************		 
c*********************subroutine macro_variables******************
c*****************************************************************
	      subroutine macro_variables
	      implicit double precision (a-h,o-z)
	      parameter(nxd=1001,nyd=401,ndel=600)
	      common /solution/u(nxd,nyd),v(nxd,nyd),rho(nxd,nyd),
     &	                   fi(0:8,nxd,nyd)
	      common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
	      common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	      common /cxy/cix(0:8),ciy(0:8),cf(0:8)

c  calculate physical variables: rho, ux, uy:

	      do i=2,nx-1
	         do j=2,ny-1
	            if(walls(i,j).eq.1) then
	              u(i,j)=0.0
	              v(i,j)=0.0
	              rho(i,j)=2.0
	            else
	              f0=fi(0,i,j) 
	              f1=fi(1,i,j) 
	              f2=fi(2,i,j) 
	              f3=fi(3,i,j) 
	              f4=fi(4,i,j) 
	              f5=fi(5,i,j) 
	              f6=fi(6,i,j) 
	              f7=fi(7,i,j) 
	              f8=fi(8,i,j) 
	              rho(i,j)=f0+f1+f2+f3+f4+f5+f6+f7+f8
	              u(i,j)=f1+f2-f4-f5-f6+f8
	              v(i,j)=f2+f3+f4-f6-f7-f8
	            endif
	         enddo
          enddo
c    set symmetry condition for uy at the symmetry line:
	      v(1:nx,2)=0.0
	      v(1:nx,ny-1)=0.0
	      u(1,1:ny)=vt
	      v(1,1:ny)=0.0
		 
	      RETURN
	      END
	
c**************************************************************
c****************difference************************************	
c**************************************************************
	     subroutine difference(un,err2)
	     implicit double precision (a-h,o-z)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension un(nxd,nyd)
	     common /solution/u(nxd,nyd),v(nxd,nyd),rho(nxd,nyd),
     &	                   fi(0:8,nxd,nyd)
	     common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /cxy/cix(0:8),ciy(0:8),cf(0:8)
	     common /maxerr/err_max,i_max,j_max
	     sum=0
	     sux=0
	     err_max=0
	     do i=2,nx-1
	        do j=2,ny-1
	         if(walls(i,j).eq.0) then
	           uu=u(i,j)
	           diff=uu-un(i,j)
	           if(abs(diff).gt.err_max) then
	             err_max=diff
	   	         i_max=i
		         j_max=j
	           endif
	           sum=sum+diff*diff
	           sux=sux+uu*uu
	         endif
	        enddo
	     enddo
	     err2=sqrt(sum/sux)
	     do i=2,nx-1
	       do j=2,ny-1
	         un(i,j)=u(i,j)
	       enddo
	     enddo
		 
	     RETURN
	     END

c********************************************************************
c*********************wall_finding***********************************
c********************************************************************
	     subroutine wall_finding(radius,imid,jmid,nx,ny,walls,wb)
	     implicit double precision (a-h,o-z)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension walls(nxd,nyd),wb(nxd,nyd)

c   determin boundary and wall nodes	
         walls(1:nx,1:ny)=0
	     wb(1:nx,1:ny)=0.0
	     xcnt=imid
	     ycnt=jmid
c  identifing the solid nodes:
	     forall(i=1:nx,j=1:ny,sqrt((xcnt-i)**2+(ycnt-j)**2)
     &		    .le.radius) walls(i,j)=1
c  identifying the boundary nodes:
	     forall(i=1:nx,j=1:ny,walls(i,j).eq.1.and.
     &		(walls(i+1,j)+walls(i+1,j+1)
     &	 	+walls(i,j+1)+walls(i-1,j+1)
     &		+walls(i-1,j)+walls(i-1,j-1)
     &     	+walls(i,j-1)+walls(i+1,j-1)-8).lt.0) wb(i,j)=1
	     do i=1,nx
	        xx=i-imid
	        do j=1,ny
		       yy=j-jmid
		       if(wb(i,j).eq.1.0) write(24,*) 'i=',i,' j=',j
	        enddo
	     enddo
         write(6,*)'finishing wall boundary identification'

	     RETURN
	     END
		 
c**************************************************************************		 
c*********************subroutine initialize********************************
c**************************************************************************
	     subroutine initialize(tau_c)
	     implicit double precision (a-h,o-z)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension feq(0:8)
	     common /solution/vx(nxd,nyd),vy(nxd,nyd),rho(nxd,nyd),
     &			 fi(0:8,nxd,nyd)
	     common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /cxy/cix(0:8),ciy(0:8),cf(0:8)
	     common /rev/kreverse(8)
	     common /constants/c0,c1,c2,c3,c4,c5,c6,c9,c12,c36

c  evaluate the equilibrium distribution function at walls(i,j)=0

	     am=2
	     c_to_f=tau/am/tau_c

         do jx=2,nx-1,2
	        do 600 jy=2,ny-1,2
	           vxx=vx(jx,jy)
	           vyy=vy(jx,jy)

	           vmag=1.5*(vxx*vxx + vyy*vyy)
	           rh_vmag=rho(jx,jy)-vmag

	           g1=vxx 
	           g2=vxx+vyy 
	           g3=vyy 
	           g4=-vxx+vyy 
	           g5=-vxx 
	           g6=-g2 
	           g7=-vyy 
	           g8=-g4

	           feq(0)=c0*rh_vmag
	           feq(1)=c1*(rh_vmag+g1*(c3+c4*g1))
	           feq(2)=c2*(rh_vmag+g2*(c3+c4*g2))
	           feq(3)=c1*(rh_vmag+g3*(c3+c4*g3))
	           feq(4)=c2*(rh_vmag+g4*(c3+c4*g4))
	           feq(5)=c1*(rh_vmag+g5*(c3+c4*g5))
	           feq(6)=c2*(rh_vmag+g6*(c3+c4*g6))
	           feq(7)=c1*(rh_vmag+g7*(c3+c4*g7))
	           feq(8)=c2*(rh_vmag+g8*(c3+c4*g8))
	      
	           if(walls(i,j).eq.1) goto 600
	           do k=0,8
	             fi(k,jx,jy)=feq(k) + c_to_f*(fi(k,jx,jy)-feq(k))
	           enddo
	      
600	        continue
	     enddo

c now perform interpolation for fi
	     do i=3,nx-2,2
	        do 76 j=2,ny-1,2
	           if(walls(i,j).eq.1) goto 76
	           do k=0,8
		          fi(k,i,j) = 0.5*(fi(k,i+1,j)+fi(k,i-1,j))
	           enddo		 
76	        continue 	
	     enddo
		 
	     do i=2,nx-1,2
	        do 86 j=3,ny-2,2
	           if(walls(i,j).eq.1) goto 86
	           do k=0,8
		          fi(k,i,j) = 0.5*(fi(k,i,j+1)+fi(k,i,j-1))
	           enddo		 
86	        continue 	
	     enddo
		 
	     do i=3,nx-2,2
	        do 96 j=3,ny-2,2
	           if(walls(i,j).eq.1) goto 96
	           do k=0,8
		          fi(k,i,j) = 0.25*(fi(k,i+1,j+1)+fi(k,i+1,j-1)
     &				  +fi(k,i-1,j+1)+fi(k,i-1,j-1))
	           enddo		 
96	        continue 	
	     enddo		

c now handle the near wall part:
c  along vertica grids:
	     do i=2,nx-1,2
	        do 248 j=3,ny-2,2
	           if(walls(i,j).eq.1) goto 248
	           jp=j+1
	           jm=j-1
	           if(walls(i,jp).eq.1) then
		         fi(0:8,i,j)=2*fi(0:8,i,j-1)-fi(0:8,i,j-2)
	           endif
	           if(walls(i,jm).eq.1) then
		         fi(0:8,i,j)=2*fi(0:8,i,j+1)-fi(0:8,i,j+2)
	           endif
248	        continue
	     enddo
c  along horizontal grids:
	     do i=3,nx-2,2
	        ip=i+1
	        im=i-1
	        do 258 j=2,ny-1,2
	           if(walls(i,j).eq.1) goto 258
	           if(walls(ip,j).eq.1) then
		         fi(0:8,i,j)=2*fi(0:8,i-1,j)-fi(0:8,i-2,j)
	           endif
	           if(walls(im,j).eq.1) then
		         fi(0:8,i,j)=2*fi(0:8,i+1,j)-fi(0:8,i+2,j)
	           endif
258	        continue
	     enddo

	     do i=3,nx-2,2
	        ip=i+1
	        im=i-1
	        do 269 j=3,ny-2,2     
	           if(walls(i,j).eq.1) goto 268
	           jp=j+1
	           jm=j-1
	           if(i.eq.145.and.j.eq.139) then
		          write(41,255) fi(0,i-2,j),fi(0,i-1,j),
     &				fi(0,i,j-1),fi(0,i,j-2)
		          write(42,255) fi(1,i-2,j),fi(1,i-1,j),
     &				fi(1,i,j-1),fi(1,i,j-2)
	           endif
	           if(i.eq.145.and.j.eq.145) then
		          write(41,255) fi(0,i-2,j),fi(0,i-1,j),
     &				fi(0,i,j+1),fi(0,i,j+2)
		          write(42,255) fi(1,i-2,j),fi(1,i-1,j),
     &				fi(1,i,j+1),fi(1,i,j+2)
	           endif
255	           FORMAT(2x,4f12.7)
	           swv=99
	           swh=99
	           sw4=99
	           sw_diag=walls(im,jp)+walls(ip,jm)+walls(ip,jp)
     &                 +walls(im,jm)
	           if(sw_diag.eq.0) goto 268
	           sw4=walls(i,jp)+walls(i,jm)+walls(ip,j)+walls(im,j)
	           if(sw4.eq.0) then
		          fi(0:8,i,j)=0.25*(fi(0:8,i,jp)+fi(0:8,i,jm)
     &		 		 +fi(0:8,im,j)+fi(0:8,ip,j))
		          goto 268
	           endif
	           swh=walls(ip,j)+walls(im,j)
	           if(swh.eq.0) then
		         fi(0:8,i,j)=0.5*(fi(0:8,ip,j)+fi(0:8,im,j))
	             goto 268
	           endif
	           swv=walls(i,jp)+walls(i,jm)
	           if(swv.eq.0) then
		         fi(0:8,i,j)=0.5*(fi(0:8,i,jp)+fi(0:8,i,jm))
	             goto 268
	           endif
268	           kk=0
269	        continue
	     enddo
		 
	     do i=3,nx-2,2
	        ip=i+1
	        im=i-1
	        do 279 j=3,ny-2,2	     
	           if(walls(i,j).eq.1) goto 278
	           jp=j+1
	           jm=j-1
	           if(walls(ip,j).eq.1.and.walls(i,jm).eq.1) then
		         fi(0:8,i,j)=.5*(2*fi(0:8,im,j)-fi(0:8,i-2,j)+
     &			        2*fi(0:8,i,jp)-fi(0:8,i,j+2))	      
	             goto 278
	           endif
	           if(walls(ip,j).eq.1.and.walls(i,jp).eq.1) then
		         fi(0:8,i,j)=.5*(2*fi(0:8,im,j)-fi(0:8,i-2,j)+
     &			        2*fi(0:8,i,jm)-fi(0:8,i,j-2))	      
	             goto 278
	           endif
	           if(walls(im,j).eq.1.and.walls(i,jp).eq.1) then
		         fi(0:8,i,j)=.5*(2*fi(0:8,ip,j)-fi(0:8,i+2,j)+
     &			        2*fi(0:8,i,jm)-fi(0:8,i,j-2))	      
	             goto 278
	           endif
	           if(walls(im,j).eq.1.and.walls(i,jm).eq.1) then
		         fi(0:8,i,j)=.5*(2*fi(0:8,ip,j)-fi(0:8,i+2,j)+
     &			        2*fi(0:8,i,jp)-fi(0:8,i,j+2))	      
	             goto 278
	           endif
278	           kw=0
279	        continue
	     enddo
		 
235	     FORMAT(3x,2i4, 3x,4f7.2)
237	     FORMAT(2x,5f13.7)
238	     FORMAT(x,i4,5f13.7)
239	     FORMAT(17x,4f13.7)
98	     FORMAT(x,2i5,x,9f12.7)
         write(6,*)'done initialization'
		 
	     RETURN
	     END

c***********************************************************************
c******************subroutine initial_cond******************************
c***********************************************************************
cc initialize physical field for flow over a cylinder: rho, ux, uy
	     subroutine initial_cond(key,nt_last)
	     implicit double precision (a-h,o-z)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension af(0:8),un(nxd,nyd)
	     common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
         common /solution/ux(nxd,nyd),uy(nxd,nyd),rho(nxd,nyd),
     &			 fi(0:8,nxd,nyd)
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /bk_c/nx_c,ny_c,imid_c,jmid_c,tau_c,radius_c
	     common /cxy/cix,ciy,cf(0:8)
         common /con/s2,s3,s5,s7,s8,s9
	     common /maxerr/err_max,i_max,j_max
	     common /rev/kreverse
	     common /constants/c0,c1,c2,c3,c4,c5,c6,c9,c12,c36
	     INTEGER  kreverse(8)/5,6,7,8,1,2,3,4/
	     real*8 cix(0:8)/0.,1.,1.,0.,-1.,-1.,-1.,0.,1.0/
	     real*8 ciy(0:8)/0.,0.,1.,1.,1.,0.,-1.,-1.,-1.0/

cc-------using potential flow field------
	     if(key.eq.0) then  
	        visc=(2*tau-1)/6.
	        vt=Re*visc/(2*radius)
	   
	        jmid=(ny-1)/2+1
	        xcnt=imid
	        ycnt=jmid
	        write(12,20) tau,radius,visc,vt,nx,ny,Re,imid,jmid
	        den_0 = 1.0
           
	        call wall_finding(radius,imid,jmid,nx,ny,walls,wb)
	        call finding_dels
	
            do i=1,nx
	           xx=i-imid
	           do j=1,ny
	             if(walls(i,j).eq.0) then
		         yy=j-jmid
		         rrr=sqrt((xcnt-i)**2+(ycnt-j)**2)
		         sinq=yy/rrr
		         cosq=xx/rrr
		         rrr=rrr/radius
		         ur=vt*(1-1./rrr**2)*cosq
		         uq=-vt*(1+1./rrr**2)*sinq

		         ux(i,j)=ur*cosq-uq*sinq
		         uy(i,j)=ur*sinq+uq*cosq
 	             rho(i,j)=0-1.5*(vt**2-ur*ur-uq*uq)
	             endif 
	           enddo
	        enddo
	        ux(1,:)=vt
	        uy(1,:)=0.0
	        write(6,*)'finishing vel specify from potential solution'

            do jx=1,nx
	        do jy=1,ny
	           vxx=ux(jx,jy)
	           vyy=uy(jx,jy)

c              vmag=1.0-1.5*(vxx*vxx + vyy*vyy)
               vmag=-1.5*(vxx*vxx + vyy*vyy)
	           rh=rho(jx,jy)

	           g0i0=c0*(rh+vmag)
	           g1=cix(1)*vxx+ciy(1)*vyy
	           g2=cix(2)*vxx+ciy(2)*vyy
	           g3=cix(3)*vxx+ciy(3)*vyy
	           g4=cix(4)*vxx+ciy(4)*vyy
	           g5=cix(5)*vxx+ciy(5)*vyy
	           g6=cix(6)*vxx+ciy(6)*vyy
	           g7=cix(7)*vxx+ciy(7)*vyy
	           g8=cix(8)*vxx+ciy(8)*vyy

	           g0i1=c1*(rh+vmag+c3*g1+c4*g1*g1)
	           g0i2=c2*(rh+vmag+c3*g2+c4*g2*g2)
	           g0i3=c1*(rh+vmag+c3*g3+c4*g3*g3)
	           g0i4=c2*(rh+vmag+c3*g4+c4*g4*g4)
	           g0i5=c1*(rh+vmag+c3*g5+c4*g5*g5)
	           g0i6=c2*(rh+vmag+c3*g6+c4*g6*g6)
	           g0i7=c1*(rh+vmag+c3*g7+c4*g7*g7)
	           g0i8=c2*(rh+vmag+c3*g8+c4*g8*g8)
	
	           fi(0,jx,jy)=g0i0
	           fi(1,jx,jy)=g0i1
	           fi(2,jx,jy)=g0i2
	           fi(3,jx,jy)=g0i3
	           fi(4,jx,jy)=g0i4
	           fi(5,jx,jy)=g0i5
	           fi(6,jx,jy)=g0i6
	           fi(7,jx,jy)=g0i7
	           fi(8,jx,jy)=g0i8
	        enddo
	        enddo
	        write(6,*)'finishing 1st collision for key=0'
			
	        call streaming(fi)
	        call macro_variables
			
	        do i=1,nx
	        do j=1,ny
	           un(i,j)=ux(i,j)
            enddo
            enddo
	        write(6,*)'finishing 1st macro_var for key=0'
	     endif
		 
cc-------start from saved information---
	     if(key.eq.1) then
	        read(22,*)tau,vt,Re,radius,nx,ny,imid,jmid,nt_last
	        s8=1./tau
	        s9=s8
	        jmid=(ny-1)/2+1
      	    visc=(2*tau-1)/6.
	        write(12,20) tau,radius,visc,vt,nx,ny,Re,imid,jmid
	        write(6,*)'reading saved velocity, density and fi fields'
	        ntotal=ny*nx
	        do ir=1,ntotal
	           read(22,*) i,j,ux(i,j),uy(i,j),rho(i,j)
	        enddo
	        read(32,*)tau, vt,Re,radius,nx,ny,imid,jmid,nt_last
	        do ir=1,ntotal
	           read(32,*) i,j,fi(0,i,j),fi(1,i,j),fi(2,i,j),
     &		             fi(3,i,j),fi(4,i,j),fi(5,i,j),
     &			         fi(6,i,j),fi(7,i,j),fi(8,i,j)
	        enddo
	        call wall_finding(radius,imid,jmid,nx,ny,walls,wb)
	        call finding_dels
	     endif 
		 
cc-------perform nested iteration right after reading the saved field.----
	     if(key.eq.2)then 
	        read(22,*)tau_c,vt,Re,radius_c,nx_c,ny_c,imid_c,jmid_c,nt_last
   	        am=2
       	    tau=0.5 + am*(tau_c - 0.5)
	        s8=1./tau
	        s9=s8
	        radius=2*radius_c
            nx=(nx_c-3)*am+3
            ny=(ny_c-3)*am+3
            imid=imid_c*am-2
            jmid=(ny-1)/2+1
      	    visc=(2*tau-1)/6.
	        write(12,20) tau,radius,visc,vt,nx,ny,Re,imid,jmid

   	        rho(1:nx,1:ny)=0
            ux(1:nx,1:ny)=0
            ux(1,1:ny)=vt
            uy(1,1:ny)=0
            fi(0:8,1:nx,1:ny)=0

	        call wall_finding(radius,imid,jmid,nx,ny,walls,wb)
	        call finding_dels

	        ntotal=nx_c*ny_c
            do ir=1,ntotal
               read(22,*) i,j,au,av,ar
	           ix=2*i-2
               jy=2*j-2
	           if(ix.eq.0) ix=1
	           if(jy.eq.0) jy=1
	           ux(ix,jy)=au
	           uy(ix,jy)=av
	           rho(ix,jy)=ar
            enddo
            read(32,*) tau_c,vt,Re,del_c,nx_c,ny_c,imid_c,jmid_c,nt_last
            do ir=1,ntotal
               read(32,*) i,j,(af(k),k=0,8)
	           ix=2*i-2
               jy=2*j-2
	           if(ix.eq.0) ix=1
	           if(jy.eq.0) jy=1
	           fi(0:8,ix,jy)=af(0:8)
   	        enddo

   	        write(6,*)'tau_c=',tau_c
c	        call initialize(tau_c) 
c	        write(6,*)'done initialization for key=2'
      	    call macro_variables
            write(6,*)'done  macro for key=2'
	     endif
cc-------finishing initialization-----
	     write(12,*)'finishing initialization; here is the ux(imid,y)'
 	     i_front=imid-radius
	     do j=2,ny-1
	        write(12,406) j,float(j-jmid)/radius,ux(imid,j)/vt,
     &			          uy(imid,j),rho(imid,j)
c	        write(15,406) j,float(j-jmid)/radius,ux(i_front,j)/vt,
c     &			          uy(i_front,j),rho(i_front,j)
	     enddo

		 
20	     FORMAT(3x,'tau=',f10.4,2x,'radius=',f6.2,2x,'visc=',f10.6,2x,
     &          'vt=',e13.5/3x,'nx=',i3,3x,'ny=',i3/3x,'Re=',f12.4,2x,
     &          'imid=',i3,2x,'jmid=',i3/)
406	     FORMAT(2x,i3,4f12.7)

	     RETURN
	     END

c***************************************************************************		 
c********************subroutine c_2_f_transfer******************************		 
c***************************************************************************	
	     subroutine c_2_f_transfer(un)
	     implicit double precision (a-h,o-z)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension un(nxd,nyd),temp(nxd,nyd)
	     common /wall/wb(nxd,nyd),walls(nxd,nyd),dels(ndel),k_del_max
	     common /solution/ux(nxd,nyd),uy(nxd,nyd),rho(nxd,nyd),
     &			          fi(0:8,nxd,nyd)
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /bk_c/nx_c,ny_c,imid_c,jmid_c,tau_c,radius_c

   	     am=2
         tau=0.5 + am*(tau_c - 0.5)
	     radius=2*radius_c
         nx=(nx_c-3)*am+3
         ny=(ny_c-3)*am+3
         imid=imid_c*am-2
         jmid=(ny-1)/2+1

         visc=(2*tau-1)/6.

	     call wall_finding(radius,imid,jmid,nx,ny,walls,wb)
	     call finding_dels

	     do i=1,nx_c
	        do j=1,ny_c
               temp(i,j)=ux(i,j)
	           ux(i,j)=0
	        enddo
	     enddo
	     do i=1,nx_c
	        do j=1,ny_c
	           ix=2*i-2
               jy=2*j-2
	           if(ix.eq.0) ix=1
	           if(jy.eq.0) jy=1
	           ux(ix,jy)=temp(i,j)
	        enddo
	     enddo
	     do i=1,nx_c
	        do j=1,ny_c
               temp(i,j)=uy(i,j)
	           uy(i,j)=0
	        enddo
	     enddo
	     do i=1,nx_c
	        do j=1,ny_c
	           ix=2*i-2
               jy=2*j-2
	           if(ix.eq.0) ix=1
	           if(jy.eq.0) jy=1
	           uy(ix,jy)=temp(i,j)
	        enddo
	     enddo	
	     do i=1,nx_c
	        do j=1,ny_c
               temp(i,j)=rho(i,j)
	           rho(i,j)=0
	        enddo
	     enddo
	     do i=1,nx_c
	        do j=1,ny_c
	           ix=2*i-2
               jy=2*j-2
	           if(ix.eq.0) ix=1
	           if(jy.eq.0) jy=1
	           rho(ix,jy)=temp(i,j)
	        enddo
	     enddo	
		 
	     do k=0,8
	        do i=1,nx_c
	           do j=1,ny_c
                 temp(i,j)=fi(k,i,j)
		         fi(k,i,j)=0
	           enddo
   	        enddo
	        do i=1,nx_c
	           do j=1,ny_c
	              ix=2*i-2
                  jy=2*j-2
	              if(ix.eq.0) ix=1
	              if(jy.eq.0) jy=1
	              fi(k,ix,jy)=temp(i,j)
	           enddo
   	        enddo	
	     enddo	
		         
	      call initialize(tau_c) 
	      write(6,*)'done initialization for key=3'
          call macro_variables
	      write(6,*)'done  macro for key=3'
   	      do i=1,nx
	      do j=1,ny
	         un(i,j)=ux(i,j)
	      enddo
   	      enddo
		  
	     RETURN
	     END
	
c*************************************************************	
c******************subroutine streaming***********************
c*************************************************************
	     subroutine streaming(fi)
	     implicit double precision (a-h,o-z)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     dimension fi(0:8,nxd,nyd)
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid

 	      do i=2,nx-1
             do j=2,ny-1
		        fi(5,i,j)=fi(5,i+1,j)
	         enddo
	      enddo
 	      do i=nx-1,2,-1
             do j=2,ny-1
		        fi(1,i,j)=fi(1,i-1,j)
	         enddo
	      enddo
 	      do i=2,nx-1
		     do j=2,ny-1
		        fi(7,i,j)=fi(7,i,j+1)
	         enddo
	      enddo
 	      do i=2,nx-1
		     do j=ny-1,2,-1
		        fi(3,i,j)=fi(3,i,j-1)
	         enddo
	      enddo
 	      do i=2,nx-1
		     do j=2,ny-1
		        fi(6,i,j)=fi(6,i+1,j+1)
	         enddo
	      enddo
 	      do i=nx-1,2,-1
		     do j=ny-1,2,-1
		        fi(2,i,j)=fi(2,i-1,j-1)
	         enddo
	      enddo
 	      do i=2,nx-1
		     do j=ny-1,2,-1
		        fi(4,i,j)=fi(4,i+1,j-1)
	         enddo
	      enddo
 	      do i=nx-1,2,-1
		     do j=2,ny-1
		        fi(8,i,j)=fi(8,i-1,j+1)
	         enddo
	      enddo
		  
	      RETURN
	      END
		  
c***************************************************************************
c********************subroutine stress_y************************************  
c***************************************************************************
	     subroutine stress_y(i_stress,sig_xx,sig_xy,sig_yy)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     implicit double precision (a-h,o-z)
	     dimension sig_xx(nxd),sig_xy(nxd),sig_yy(nxd),
     ^		cc(0:8),f_neq(0:8),feq(0:8)
	     common /solution/vx(nxd,nyd),vy(nxd,nyd),rho(nxd,nyd),
     &			 fi(0:8,nxd,nyd)
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /cxy/cix(0:8),ciy(0:8),cf(0:8)
	     common /con/ s2,s3,s5,s7,s8,s9
	     common /constants/c0,c1,c2,c3,c4,c5,c6,c9,c12,c36
	     common /rev/kreverse(8)
	
         conv=1.0-0.5/tau
	     do k=1,8
	        cc(k)=cix(k)**2+ciy(k)**2
	     enddo
	
	     jx=i_stress
	     nxm=nx-1

	     do jy=2,ny-1
	      vxx=vx(jx,jy)
	      vyy=vy(jx,jy)

          vmag=-1.5*(vxx*vxx + vyy*vyy)
	      rh=rho(jx,jy)

	      g0i0=c0*(rh+vmag)
	      g1=cix(1)*vxx+ciy(1)*vyy
	      g2=cix(2)*vxx+ciy(2)*vyy
	      g3=cix(3)*vxx+ciy(3)*vyy
	      g4=cix(4)*vxx+ciy(4)*vyy
	      g5=cix(5)*vxx+ciy(5)*vyy
	      g6=cix(6)*vxx+ciy(6)*vyy
	      g7=cix(7)*vxx+ciy(7)*vyy
	      g8=cix(8)*vxx+ciy(8)*vyy

	      feq(1)=c1*(rh+vmag+c3*g1+c4*g1*g1)
	      feq(2)=c2*(rh+vmag+c3*g2+c4*g2*g2)
	      feq(3)=c1*(rh+vmag+c3*g3+c4*g3*g3)
	      feq(4)=c2*(rh+vmag+c3*g4+c4*g4*g4)
	      feq(5)=c1*(rh+vmag+c3*g5+c4*g5*g5)
	      feq(6)=c2*(rh+vmag+c3*g6+c4*g6*g6)
	      feq(7)=c1*(rh+vmag+c3*g7+c4*g7*g7)
	      feq(8)=c2*(rh+vmag+c3*g8+c4*g8*g8)
	
	      s_xx=0.0
	      s_xy=0.0
	      s_yy=0.0

	      do k=1,8
	         f_neq(k)=fi(k,jx,jy)-feq(k)
		     s_xx=s_xx+ f_neq(k)*(cix(k)*cix(k)-0.5*cc(k))
		     s_xy=s_xy+ f_neq(k)*cix(k)*ciy(k)
		     s_yy=s_yy+ f_neq(k)*(ciy(k)*ciy(k)-0.5*cc(k))
	      enddo
	      sig_xx(jy)=-s_xx*conv
	      sig_xy(jy)=-s_xy*conv
	      sig_yy(jy)=-s_yy*conv

	     enddo
		 
	     RETURN
	     END
		 
c***************************************************************************
c********************subroutine stress_x************************************  
c***************************************************************************
	     subroutine stress_x(j_stress,sig_xx,sig_xy,sig_yy)
	     parameter(nxd=1001,nyd=401,ndel=600)
	     implicit double precision (a-h,o-z)
	     dimension sig_xx(nxd),sig_xy(nxd),sig_yy(nxd),
     ^		cc(0:8),f_neq(0:8),feq(0:8)
	     common /solution/vx(nxd,nyd),vy(nxd,nyd),rho(nxd,nyd),
     &			 fi(0:8,nxd,nyd)
	     common /bk1/nx,ny,tau,vt,Re,radius,imid,jmid
	     common /cxy/cix(0:8),ciy(0:8),cf(0:8)
	     common /con/ s2,s3,s5,s7,s8,s9
	     common /constants/c0,c1,c2,c3,c4,c5,c6,c9,c12,c36
	     common /rev/kreverse(8)
	
         conv=1.0-0.5/tau
	     do k=1,8
	        cc(k)=cix(k)**2+ciy(k)**2
	     enddo
	
	     jy=j_stress
	     nxm=nx-1

	     do jx=2,nx-1
	        vxx=vx(jx,jy)
	        vyy=vy(jx,jy)

            vmag=-1.5*(vxx*vxx + vyy*vyy)
	        rh=rho(jx,jy)

	        g0i0=c0*(rh+vmag)
	        g1=cix(1)*vxx+ciy(1)*vyy
	        g2=cix(2)*vxx+ciy(2)*vyy
	        g3=cix(3)*vxx+ciy(3)*vyy
	        g4=cix(4)*vxx+ciy(4)*vyy
	        g5=cix(5)*vxx+ciy(5)*vyy
	        g6=cix(6)*vxx+ciy(6)*vyy
	        g7=cix(7)*vxx+ciy(7)*vyy
	        g8=cix(8)*vxx+ciy(8)*vyy

	        feq(1)=c1*(rh+vmag+c3*g1+c4*g1*g1)
	        feq(2)=c2*(rh+vmag+c3*g2+c4*g2*g2)
	        feq(3)=c1*(rh+vmag+c3*g3+c4*g3*g3)
	        feq(4)=c2*(rh+vmag+c3*g4+c4*g4*g4)
	        feq(5)=c1*(rh+vmag+c3*g5+c4*g5*g5)
	        feq(6)=c2*(rh+vmag+c3*g6+c4*g6*g6)
	        feq(7)=c1*(rh+vmag+c3*g7+c4*g7*g7)
	        feq(8)=c2*(rh+vmag+c3*g8+c4*g8*g8)
	
	        s_xx=0.0
	        s_xy=0.0
	        s_yy=0.0

	        do k=1,8
		       f_neq(k)=fi(k,jx,jy)-feq(k)
		       s_xx=s_xx+ f_neq(k)*(cix(k)*cix(k)-0.5*cc(k))
		       s_xy=s_xy+ f_neq(k)*cix(k)*ciy(k)
		       s_yy=s_yy+ f_neq(k)*(ciy(k)*ciy(k)-0.5*cc(k))
	        enddo

	        if (jx.ge.130.and.jx.le.181) then
	  	       write(16,28) jx,(f_neq(k),k=1,5)
		       write(26,58) jx,s_xx,f_neq(6),f_neq(7),f_neq(8)
	        endif
28  	 	FORMAT(x,i3,2x,5f13.7)
58 		    FORMAT(3x,i3,3x,4f14.8)
	        sig_xx(jx)=-s_xx*conv
	        sig_xy(jx)=-s_xy*conv
	        sig_yy(jx)=-s_yy*conv

	     enddo
		 
	     RETURN
	     END

