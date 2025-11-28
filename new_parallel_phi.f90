module my_subs2
implicit none

contains

	FUNCTION cross(a, b)
	  Real(8), DIMENSION(3) :: cross
	  Real(8), DIMENSION(3), INTENT(IN) :: a, b

	  cross(1) = a(2) * b(3) - a(3) * b(2)
	  cross(2) = a(3) * b(1) - a(1) * b(3)
	  cross(3) = a(1) * b(2) - a(2) * b(1)
	END FUNCTION cross

	Subroutine routateHodd(angle_theta,angle_phi,TR,ham_r,Ndim,nrpts,ham,ham_odd,zpos1)			
		Integer,Intent(In)::Ndim,nrpts
		Real(8), Intent(In)::angle_phi(2),angle_theta(2),zpos1(Ndim)
		complex*16, INTENT(IN) ::ham_r(Ndim,Ndim,nrpts),TR(Ndim,Ndim)
		complex*16, INTENT(InOut)::ham(Ndim,Ndim,nrpts),ham_odd(Ndim,Ndim,nrpts)
		complex*16::Rot(Ndim,Ndim),TRH(Ndim,Ndim),TRHTR(Ndim,Ndim),TROdd(Ndim,Ndim),TREven(Ndim,Ndim),RotTOdd(Ndim,Ndim),RotTOddRot(Ndim,Ndim)
		Integer::irpt,ii,jj
		Real(8)::Pi
		PI=4.D0*DATAN(1.D0)
		Rot=0d0


                do ii=1, Ndim/2
                        if (zpos1(ii)<=0.5d0) then
                                Rot(ii  ,ii)=Dcmplx(Cos(angle_theta(1)/2d0)*Cos(angle_phi(1)/2d0),-Cos(angle_theta(1)/2d0)*Sin(angle_phi(1)/2d0))
                                Rot(ii,ii+Ndim/2)=Dcmplx(-Sin(angle_theta(1)/2d0)*Cos(angle_phi(1)/2d0),Sin(angle_theta(1)/2d0)*Sin(angle_phi(1)/2d0))
                                Rot(ii+Ndim/2,ii)=Dcmplx(Sin(angle_theta(1)/2d0)*Cos(angle_phi(1)/2d0),Sin(angle_theta(1)/2d0)*Sin(angle_phi(1)/2d0))
                                Rot(ii+Ndim/2,ii+Ndim/2)=Dcmplx(Cos(angle_theta(1)/2d0)*Cos(angle_phi(1)/2d0),Cos(angle_theta(1)/2d0)*Sin(angle_phi(1)/2d0))
                        else
                                Rot(ii,ii)=Dcmplx(Cos(angle_theta(2)/2d0)*Cos(angle_phi(2)/2d0),-Cos(angle_theta(2)/2d0)*Sin(angle_phi(2)/2d0))
                                Rot(ii,ii+Ndim/2)=Dcmplx(-Sin(angle_theta(2)/2d0)*Cos(angle_phi(2)/2d0),Sin(angle_theta(2)/2d0)*Sin(angle_phi(2)/2d0))
                                Rot(ii+Ndim/2,ii)=Dcmplx(Sin(angle_theta(2)/2d0)*Cos(angle_phi(2)/2d0),Sin(angle_theta(2)/2d0)*Sin(angle_phi(2)/2d0))
                                Rot(ii+Ndim/2,ii+Ndim/2)=Dcmplx(Cos(angle_theta(2)/2d0)*Cos(angle_phi(2)/2d0),Cos(angle_theta(2)/2d0)*Sin(angle_phi(2)/2d0))
                                !Rot(ii  ,ii)  =Dcmplx(1d0,0d0)
                               !Rot(ii,ii+Ndim/2)=Dcmplx(0d0,0d0)
                               !Rot(ii+Ndim/2,ii)  =Dcmplx(0d0,0d0)
                                !Rot(ii+Ndim/2,ii+Ndim/2)=Dcmplx(1d0,0d0)
                        endif
                enddo



		do irpt=1,nrpts
			!~~~Get time-reversal odd/even parts from the ham_r, odd part would be magnetization~~~~~
			TRH=0d0
			TRHTR=0d0
			TROdd=0d0
			TREven=0d0
			CALL zgemm ('T', 'N', Ndim, Ndim, Ndim, Dcmplx(1d0,0d0),TR, Ndim, ham_r(:,:,irpt), Ndim, Dcmplx(0d0,0d0),TRH, Ndim)
			CALL zgemm ('N', 'N', Ndim, Ndim, Ndim, Dcmplx(1d0,0d0),TRH, Ndim, TR, Ndim, Dcmplx(0d0,0d0),TRHTR, Ndim)		
			do ii=1,Ndim
			do jj=1,Ndim
				TROdd(ii,jj)=0.5d0*(ham_r(ii,jj,irpt)-Conjg(TRHTR(ii,jj)))
				TREven(ii,jj)=0.5d0*(ham_r(ii,jj,irpt)+Conjg(TRHTR(ii,jj)))
			!	if (irpt==332) write(19,"(5i4,20f22.14)") irvec(:,irpt),ii,jj,TRodd(ii,jj),TReven(ii,jj),ham_r(ii,jj,irpt)
			enddo
			enddo		
			CALL zgemm ('N', 'N', Ndim, Ndim, Ndim, Dcmplx(1d0,0d0),Rot, Ndim, TROdd, Ndim, Dcmplx(0d0,0d0),RotTOdd, Ndim)
			CALL zgemm ('N', 'C', Ndim, Ndim, Ndim, Dcmplx(1d0,0d0),RotTOdd, Ndim, Rot, Ndim, Dcmplx(0d0,0d0),RotTOddRot, Ndim)
				
			do ii=1,Ndim
			do jj=1,Ndim
				ham(ii,jj,irpt)=TREven(ii,jj)+RotTOddRot(ii,jj)
				ham_odd(ii,jj,irpt)=RotTOddRot(ii,jj)
			enddo
			enddo								
		enddo
	End Subroutine routateHodd
end module my_subs2

Program MoTe2_fromWF
use my_subs2
implicit none
! Load MPI definitions   
include 'mpif.h'
integer:: size,status(MPI_STATUS_SIZE),& 
&iat,jat,counteri,counterj,ichain,jchain,counter,nstate,Nchain,Ndim,N_mu,N_smr,i_smr,i_mu,nelectron,Nphi,Nversion, ntyp,iflag,rank, Njob,jobcounter,ierror, flag,nadp,tag,Nk,Nkx,Nky,Nkz,flag_n,istep,i,j,irpt,file_unit,num_wann,nrpts,ierr,ir,jr,LWORK,INFO,Nstep,Nstep_theta,Nstep_phi,ii,jj,&
&N_phi,num_wannier_plot,num_atoms,nsp,num_species,loop,loop2,loop_b,ngs(3),&
&ngx,ngy,ngz,nat,nx,ny,nz,istat,errio,ln,loop_bb,i_center,x_center,&
&y_center,z_center,nx_center,ny_center,nz_center,whichband,whichk,ix,iy,iz,kk,nband,kk2,ixi,ixf,iyi,iyf,izi,izf,i_phi,i_theta,kcounter,theta_0,theta_n,phi_0,phi_n

character (len=33) :: seedname,header
character (len=9)  :: cdate,ctime
complex*16, allocatable :: Umirror(:,:),ham_r(:,:,:),ham_r_up(:,:,:),ham_r_dn(:,:,:),WORK(:),ham_k(:,:),hso(:,:),ham(:,:,:),ham_odd(:,:,:),ham_even(:,:,:),&
&Vy(:,:),Vx(:,:),TR(:,:),TRH(:,:),TRHTR(:,:),&
&VxU(:,:),UVxU(:,:),VyU(:,:),UVyU(:,:),&
&Sz(:,:),Sx(:,:),Sy(:,:),TREven(:,:),TROdd(:,:),&
&Lx(:,:),Ly(:,:),Lz(:,:),&
&Gammax(:,:,:),Gammay(:,:,:),Gammaz(:,:,:),&
&GammaxU(:,:,:),GammayU(:,:,:),GammazU(:,:,:),&
&UGammaxU(:,:,:),UGammayU(:,:,:),UGammazU(:,:,:),&
&SxU(:,:),SyU(:,:),SzU(:,:),&
&USxU(:,:),USyU(:,:),USzU(:,:),&
&TROddSx(:,:),TROddSy(:,:),TROddSz(:,:),SxTROdd(:,:),SyTROdd(:,:),SzTROdd(:,:),&
&Sxatom(:,:,:),Syatom(:,:,:),Szatom(:,:,:),Rot(:,:)
integer,   allocatable :: irvec(:,:),atoms_species_num(:), ndegen(:),norbs(:),natoms(:),orbs(:)
Double Precision, allocatable:: alpha(:),xpos(:),ypos(:),zpos(:),xpos1(:),ypos1(:),zpos1(:),xval(:),kpath(:,:),Uconv(:,:,:),RWORK(:),atom_pos(:,:),eg(:,:),Eigenvalue(:),&
	&fac(:,:),Fermifac(:),fac_odd(:),tetotal(:,:),density(:,:),tgammax(:,:,:),tgammay(:,:,:),tgammaz(:,:,:),tgammaoddx(:,:,:),tgammaoddy(:,:,:),tgammaoddz(:,:,:),&
	&sum_gammax(:),sum_gammay(:),sum_gammaz(:),sum_gammaoddx(:),sum_gammaoddy(:),sum_gammaoddz(:)
Real(8):: magx(2),magy(2),magz(2), k(3),a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vol,aB,kweight,ght,phasexy,tVp2,tVp0,tVz0,tVp,tVz,tVpE,tVzE,zcurrent,smr_0,smr_odd,smr,smr_step,mu_step,ef_0,Vran,minuskdotR,kxadpmesh,kyadpmesh,QSEr,etotal,alpha_Co,alpha_Pb,alpha_Te,tol,b1x,b2x,b3x,b1y,b1z,b2y,b2z,b3y,b3z,a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,&
		&kxmesh,kymesh,kzmesh,kx,ky,kz,kdotR,real_lattice(3,3),mulow,muup,mu,densitynew,&
		&r_x,r_y,r_z,phase,Pi,start,finish,Efield,ef,&
		&Temp,deltae,kxi,kxf,kyi,kyf,kzi,kzf,delta_theta,delta_phi,angle_phi(2),angle_theta(2),&
		&tGammar,tGammaEr,abstorque,newabstorque,deltatorque,threshold,kQpx,kQpy,kQpz,&
		&sum_etotal,sum_density	
character(len=256) :: char_read	
CHARACTER(LEN=1024) ::filename
Parameter(aB=0.529177210903d0)

	
	Pi=4.D0*DATAN(1.D0)
    ! Initialize MPI
        call MPI_INIT(ierror)
        ! Get the number of processes 
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
        ! Get my process number (rank) 
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        if (rank==0) write(*,*) "total process is", size

!~~~~~~~System dependent~~~~~~~~~~~~~~
        a1=[ 4.3613760d0,0d0,0d0]
        a2=[- 4.3613760d0/2d0, 4.3613760d0*sqrt(3d0)/2d0,0d0]
        a3=[0d0,0d0,64.5226130d0]
	vol=Dot_product(a1,cross(a2,a3))
	b1=2d0*Pi*cross(a2,a3)/vol
	b2=2d0*Pi*cross(a3,a1)/vol
	b3=2d0*Pi*cross(a1,a2)/vol

!~~~~~~~System dependent~~~~~~~~~~~~~~
	ntyp=3 !~~~3 types atoms Fe Ge Te
	allocate (natoms(ntyp),orbs(ntyp),stat=ierr)
	natoms=[2,4,8] !~~~4 each
	nat=sum(natoms)	
	orbs=[5,3,3] !~~~atomic projections are spd;sp;sp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    allocate(norbs(nat),stat=ierr)
	i=0
	do ii=1,ntyp
	do jj=1,natoms(ii)
		i=i+1
		norbs(i)=orbs(ii) !We have 2 spins
	enddo
	enddo
	tol=1d-5
	!ef=6.1605d0
	Temp=0.01d0!~~~About 100K
	!write(*,*) 'whichband', 'whichk' 




    call cpu_time(start)	
!=============================================================!
!!  Read the Hamiltonian in the WF basis from seedname_hr.dat   
!=============================================================!
    file_unit=18 
    
    seedname= 'wannier90_symmed'

    open(file_unit,file=trim(seedname)//'_trim_hr.dat',form='formatted',&
         status='unknown')
    read(file_unit,*) header ! Date and time
    read(file_unit,*) num_wann
    read(file_unit,*) nrpts
    allocate(irvec(3,nrpts),stat=ierr)
    if (ierr/=0) then
	write(*,*) 'Error in allocating irvec in hamiltonian_Asetup'
    endif
    allocate(ndegen(nrpts),stat=ierr)
    if (ierr/=0) then
	write(*,*) 'Error in allocating ndegen in hamiltonian_Asetup'
    endif
    allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
    if (ierr/=0) then
	write(*,*) 'Error in allocating ham_r in hamiltonian_Asetup'
    endif

    read(file_unit,'(15I5)') (ndegen(i),i=1,nrpts)
    do irpt=1,nrpts
       do i=1,num_wann
          do j=1,num_wann
		    read(file_unit,'(5I5,2F22.16)') irvec(:,irpt), jr, ir,ham_r(j,i,irpt)
          end do
       end do
    end do
    close(file_unit)

	!~~~Read atom positions~~~
	open(file_unit,file='MnBi2Te4_pos.dat',form='formatted',status='unknown')
    read(file_unit,*) nat !how many atoms 
    allocate(atom_pos(nat,3),stat=ierr)
    if (ierr/=0) write(*,*) 'Error in allocating atom_pos'
    do ii=1,nat
            read(file_unit,*) char_read,atom_pos(ii,:)
    end do
	close(file_unit)
 
	
    
	Ndim=num_wann !consider spin and Wannier90 has 2 spins already
        LWORK=5*Ndim

	mu_step=0.01d0
	smr_step=0.01d0
	N_mu=400
	N_smr=10
        smr_0=5d-3
        ef_0=-0.42653847
	allocate(tetotal(N_mu,N_smr),density(N_mu,N_smr),Rot(Ndim,Ndim),&
			&tgammax(N_mu,N_smr,2),tgammay(N_mu,N_smr,2),tgammaz(N_mu,N_smr,2),&
			&tgammaoddx(N_mu,N_smr,2),tgammaoddy(N_mu,N_smr,2),tgammaoddz(N_mu,N_smr,2),&			
			!&tgammayx(N_mu,N_smr),tgammayy(N_mu,N_smr),tgammayz(N_mu,N_smr),&
			!&tgammazx(N_mu,N_smr),tgammazy(N_mu,N_smr),tgammazz(N_mu,N_smr),&
			!&tToddxx(N_mu,N_smr),tToddxy(N_mu,N_smr),tToddxz(N_mu,N_smr),&
			!&tToddyx(N_mu,N_smr),tToddyy(N_mu,N_smr),tToddyz(N_mu,N_smr),&
			!&tToddzx(N_mu,N_smr),tToddzy(N_mu,N_smr),tToddzz(N_mu,N_smr),&
			!&tQSxx(N_mu,N_smr),tQSxy(N_mu,N_smr),tQSxz(N_mu,N_smr),&
			!&tQSyx(N_mu,N_smr),tQSyy(N_mu,N_smr),tQSyz(N_mu,N_smr),&
			!&tQSzx(N_mu,N_smr),tQSzy(N_mu,N_smr),tQSzz(N_mu,N_smr),&
			!&thsoxx(N_mu,N_smr),thsoxy(N_mu,N_smr),thsoxz(N_mu,N_smr),&
			!&thsoyx(N_mu,N_smr),thsoyy(N_mu,N_smr),thsoyz(N_mu,N_smr),&
			!&thsozx(N_mu,N_smr),thsozy(N_mu,N_smr),thsozz(N_mu,N_smr),&			
			&Uconv(Ndim,Ndim,3),ham_k(Ndim,Ndim),ham(Ndim,Ndim,nrpts),ham_odd(Ndim,Ndim,nrpts),ham_even(Ndim,Ndim,nrpts),&
			&TR(Ndim,Ndim),TROdd(Ndim,Ndim),TREven(Ndim,Ndim),TRH(Ndim,Ndim),TRHTR(Ndim,Ndim),&
			&Vx(Ndim,Ndim),Vy(Ndim,Ndim),&!Vz(Ndim,Ndim),&
			&VxU(Ndim,Ndim),UVxU(Ndim,Ndim),VyU(Ndim,Ndim),UVyU(Ndim,Ndim),&!VzU(Ndim,Ndim),UVzU(Ndim,Ndim),&
			&fac(Ndim,Ndim),Fermifac(Ndim),fac_odd(Ndim),&
			&xpos(Ndim/2),ypos(Ndim/2),zpos(Ndim/2),xpos1(Ndim),ypos1(Ndim),zpos1(Ndim),&
			&Lx(Ndim,Ndim),Ly(Ndim,Ndim),Lz(Ndim,Ndim),Sx(Ndim,Ndim),Sy(Ndim,Ndim),Sz(Ndim,Ndim),&
			&hso(Ndim,Ndim),&
			!&hsoSx(Ndim,Ndim),hsoSy(Ndim,Ndim),hsoSz(Ndim,Ndim),&
			!&Sxhso(Ndim,Ndim),Syhso(Ndim,Ndim),Szhso(Ndim,Ndim),&
			!&hsox(Ndim,Ndim),hsoy(Ndim,Ndim),hsoz(Ndim,Ndim),&			
			&Gammax(Ndim,Ndim,2),Gammay(Ndim,Ndim,2),Gammaz(Ndim,Ndim,2),&
			&SxTROdd(Ndim,Ndim),SyTROdd(Ndim,Ndim),SzTROdd(Ndim,Ndim),TROddSx(Ndim,Ndim),TROddSy(Ndim,Ndim),TROddSz(Ndim,Ndim),&
			!&Toddx(Ndim,Ndim),Toddy(Ndim,Ndim),Toddz(Ndim,Ndim),&
			!&QSx(Ndim,Ndim),QSy(Ndim,Ndim),QSz(Ndim,Ndim),&			
			!&hsoxU(Ndim,Ndim),hsoyU(Ndim,Ndim),hsozU(Ndim,Ndim),UhsoxU(Ndim,Ndim),UhsoyU(Ndim,Ndim),UhsozU(Ndim,Ndim),&
			!&ToddxU(Ndim,Ndim),ToddyU(Ndim,Ndim),ToddzU(Ndim,Ndim),UToddxU(Ndim,Ndim),UToddyU(Ndim,Ndim),UToddzU(Ndim,Ndim),&
			&GammaxU(Ndim,Ndim,2),GammayU(Ndim,Ndim,2),GammazU(Ndim,Ndim,2),UGammaxU(Ndim,Ndim,2),UGammayU(Ndim,Ndim,2),UGammazU(Ndim,Ndim,2),&
			&SxU(Ndim,Ndim),SyU(Ndim,Ndim),SzU(Ndim,Ndim),USxU(Ndim,Ndim),USyU(Ndim,Ndim),USzU(Ndim,Ndim),&
			!&QSxU(Ndim,Ndim),QSyU(Ndim,Ndim),QSzU(Ndim,Ndim),UQSxU(Ndim,Ndim),UQSyU(Ndim,Ndim),UQSzU(Ndim,Ndim),&
			&RWORK(3*Ndim-2),WORK(LWORK),Eigenvalue(Ndim),stat=ierr)

    if (ierr/=0) then
		write(*,*) 'Error in allocating ham_k in hamiltonian_Asetup'
    endif
   	call cpu_time(start)
	
!~~~Define spin matrices, time-reversal operator matrix, and mirror transformation matrix~~~
	Sz=0d0
	Sx=0d0
	Sy=0d0
	TR=0d0
	do ii=1, Ndim/2
		Sz(ii,ii)=0.5d0
		Sz(ii+Ndim/2,ii+Ndim/2)=-0.5d0
		Sx(ii,ii+Ndim/2)=0.5d0
		Sx(ii+Ndim/2,ii)=0.5d0
		Sy(ii,ii+Ndim/2)=DCMPLX(0d0,-0.5d0)
		Sy(ii+Ndim/2,ii)=DCMPLX(0d0,0.5d0)
		TR(ii,ii+Ndim/2)=-1d0
		TR(ii+Ndim/2,ii)=1d0		
	enddo

	counter=0
	do ii=1,nat
	do jj=1,norbs(ii)
		counter=counter+1
		xpos(counter)=atom_pos(ii,1)
		ypos(counter)=atom_pos(ii,2)
		zpos(counter)=atom_pos(ii,3)	
	enddo
	enddo
        xpos1=[xpos,xpos]
        ypos1=[ypos,ypos]
        zpos1=[zpos,zpos] 
        write(*,*) "positon=",zpos1
	
	Uconv=0d0
	do ii=1,Ndim
	do jj=1,Ndim
		Uconv(ii,jj,1)=xpos1(jj)-xpos1(ii)
		Uconv(ii,jj,2)=ypos1(jj)-ypos1(ii)
		Uconv(ii,jj,3)=zpos1(jj)-zpos1(ii)	
              !  write(*,*) "positon=", Uconv(ii,jj,1),Uconv(ii,jj,2),Uconv(ii,jj,3)
	enddo
       
	enddo		
       
	call cpu_time(finish)
	!if (rank==0) 
        write(*,*) "finish setup, used time=", finish-start

!~~~Rotate time-reversal odd part in spin space with \theta, \phi~~~
		!angle_phi=delta_phi*i_phi
		!write(*,*) "input angles"
		!read(*,*) angle_theta,angle_phi
                !angle_theta(1)=0.0d0*pi
                !angle_phi(1)=0.0d0*pi
                !angle_theta(2)=0.0d0*pi
                !angle_phi(2)=0.0d0*pi
		

                Nstep_theta=80
                Nstep_phi=160
                delta_phi=2.0d0*Pi/(Nstep_phi)
                delta_theta=Pi/(Nstep_theta)
                theta_0=0
                theta_n=20
                phi_n=160
                phi_0=0
                !read(*,*) angle_theta,angle_phi
                !angle_theta=angle_theta*pi
                !angle_phi=angle_phi*pi
                write(*,*) "(theta,phi)=", angle_theta,angle_phi




                Sxatom=0d0
                Syatom=0d0
                Szatom=0d0
                tgammaoddx=0d0
                tgammaoddy=0d0
                tgammaoddz=0d0
                tgammax=0d0
                tgammay=0d0
                tgammaz=0d0

                magx=0d0
                magy=0d0
                magz=0d0

       allocate(Sxatom(Ndim,Ndim,2),Syatom(Ndim,Ndim,2),Szatom(Ndim,Ndim,2),stat=ierr)
                if (ierr/=0) then
                        write(*,*) 'Error in allocating Satom in setup'
                        stop
                endif
                do ii=1,Ndim/2,1
                        if (zpos1(ii)<=0.5d0) then
                                Szatom(ii,ii,1)=0.5d0
                                Szatom(ii+Ndim/2,ii+Ndim/2,1)=-0.5d0
                                Sxatom(ii,ii+Ndim/2,1)=0.5d0
                                Sxatom(ii+Ndim/2,ii,1)=0.5d0
                                Syatom(ii,ii+Ndim/2,1)=DCMPLX(0d0,-0.5d0)
                                Syatom(ii+Ndim/2,ii,1)=DCMPLX(0d0,0.5d0)
                        else
                                Szatom(ii,ii,2)=0.5d0
                                Szatom(ii+Ndim/2,ii+Ndim/2,2)=-0.5d0
                                Sxatom(ii,ii+Ndim/2,2)=0.5d0
                                Sxatom(ii+Ndim/2,ii,2)=0.5d0
                                Syatom(ii,ii+Ndim/2,2)=DCMPLX(0d0,-0.5d0)
                                Syatom(ii+Ndim/2,ii,2)=DCMPLX(0d0,0.5d0)
                       endif
                enddo






         open(unit=1000,file="Torque_data_theta_0_20_.dat",form='formatted',status='unknown')

        do i_theta=theta_0,theta_n
                angle_theta(1)=delta_theta*i_theta !+Pi/2d0 Pi/Nstep/2d0+
                angle_theta(2)=angle_theta(1)
        do i_phi=phi_0,phi_n
                angle_phi(1)=delta_phi*i_phi
                angle_phi(2)=pi-angle_phi(1) !delta_phi*i_phi
                write(*,*) "my phi angle is", angle_phi

		call routateHodd(angle_theta,angle_phi,TR,ham_r,Ndim,nrpts,ham,ham_odd,zpos1)	


                tgammaoddx=0d0
                tgammaoddy=0d0
                tgammaoddz=0d0
                tgammax=0d0
                tgammay=0d0
                tgammaz=0d0



		tetotal=0d0		
		flag=0	
		kcounter=0
		density=0d0	
		jobcounter=0
		Nkx=200
		Nky=200
                kweight=1d0/(Nkx*Nky)!factor is 1/(N)

!	open(unit=1000,file="Torque_data_phi.dat",form='formatted',status='unknown')		
	!open(unit=2000,file="spinSx_polarization_bottom.dat",form='formatted',status='unknown')		
!	open(unit=3000,file="torque_z.dat",form='formatted',status='unknown')
!        open(unit=4000,file="velocity.dat",form='formatted',status='unknown')





        !allocate(kpath(100,3))
 
  
              do ix=0,Nkx-1
              do iy=0, Nky-1
              jobcounter=jobcounter+1
              if (mod(jobcounter-1,size)==rank) then
                kx=-b1(1)/2d0-b2(1)/2d0+ix*b1(1)/Nkx+iy*b2(1)/Nky
                ky=-b1(2)/2d0-b2(2)/2d0+ix*b1(2)/Nkx+iy*b2(2)/Nky
                kz=0d0
                !kx=0
                !ky=0
                !kz=0
                !istep=1
                !if (rank.eq.0) then
                !        read(*,*) kpath(istep,1),kpath(istep,2),kpath(istep,3)
                !endif
         
                !kx=kpath(istep,1)*b1(1)+kpath(istep,2)*b2(1)+kpath(istep,3)*b3(1)
                !ky=kpath(istep,1)*b1(2)+kpath(istep,2)*b2(2)+kpath(istep,3)*b3(2)
                !kz=kpath(istep,1)*b1(3)+kpath(istep,2)*b2(3)+kpath(istep,3)*b3(3)



                if (rank.eq.0) then
                           write(*,*) "input k in cartesian"
                           write(*,*) "k=",kx ,ky ,kz
                endif
		!write(*,*) "k=",kx,ky
	        ham_k=cmplx(0d0)
				!Vx=cmplx(0d0)
		Vx=cmplx(0d0)
		TROdd=cmplx(0d0)
                
		do irpt=1,nrpts
					do ii=1,Ndim
						do jj=1,Ndim
							r_x=(irvec(1,irpt)+Uconv(jj,ii,1))*a1(1)+(irvec(2,irpt)+Uconv(jj,ii,2))*a2(1)+(irvec(3,irpt)+Uconv(jj,ii,3))*a3(1)
							r_y=(irvec(1,irpt)+Uconv(jj,ii,1))*a1(2)+(irvec(2,irpt)+Uconv(jj,ii,2))*a2(2)+(irvec(3,irpt)+Uconv(jj,ii,3))*a3(2)
							r_z=(irvec(1,irpt)+Uconv(jj,ii,1))*a1(3)+(irvec(2,irpt)+Uconv(jj,ii,2))*a2(3)+(irvec(3,irpt)+Uconv(jj,ii,3))*a3(3)
							kdotR=kx*r_x+ky*r_y+kz*r_z
							ham_k(jj,ii)=ham_k(jj,ii)+ham(jj,ii,irpt)*DCMPLX(Cos(kdotR),Sin(kdotR))/real(ndegen(irpt))
                                                        TROdd(jj,ii)=TROdd(jj,ii)+ham_odd(jj,ii,irpt)*DCMPLX(Cos(kdotR),Sin(kdotR))/real(ndegen(irpt))
                                                        Vx(jj,ii)=Vx(jj,ii)+DCMPLX(0d0,r_x)*ham(jj,ii,irpt)*DCMPLX(Cos(kdotR),Sin(kdotR))/real(ndegen(irpt))
						end do
					end do
				end do
				!TReven=ham_k-TRodd
				!ham_k=ham_k+hso
				CALL    ZHEEV('V', 'U', Ndim, ham_k, Ndim, Eigenvalue, WORK, LWORK, RWORK, INFO)
				if (INFO .NE. 0) then
					Write(*,*) 'The algorithm failed to compute eigenvalues of ham_k. k=', kx,ky,kz
					STOP
				endif	
                                if (rank.eq.0) then
                                        write(*,*) Eigenvalue
                                endif
                                !write(*,*) "Eg",Eigenvalue
                                !~~~Torque
                                !operator = -i
                                ![S,H_odd]
                                Gammax=0d0
                                Gammay=0d0
                                Gammaz=0d0
                                GammaxU=0d0
                                GammayU=0d0
                                GammazU=0d0
                                UGammaxU=0d0
                                UGammayU=0d0
                                UGammazU=0d0

                                do i=1,2!two layers
                                        !Here we look at time-reversal odd part
                                        !torque~~~
                                        TRoddSx=0d0
                                        TRoddSy=0d0
                                        TRoddSz=0d0
                                        SxTRodd=0d0
                                        SyTRodd=0d0
                                        SzTRodd=0d0
                                        CALL zgemm ('N', 'N',Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),TROdd,Ndim, Sxatom(:,:,i),Ndim, Dcmplx(0d0,0d0),TRoddSx,Ndim)
                                        CALL zgemm ('N', 'N',Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),Sxatom(:,:,i),Ndim, TRodd,Ndim, Dcmplx(0d0,0d0),SxTRodd,Ndim)
                                        CALL zgemm ('N', 'N',Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),TRodd,Ndim, Syatom(:,:,i),Ndim, Dcmplx(0d0,0d0),TRoddSy,Ndim)
                                        CALL zgemm ('N', 'N',Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),Syatom(:,:,i),Ndim, TRodd,Ndim, Dcmplx(0d0,0d0),SyTRodd,Ndim)
                                        CALL zgemm ('N', 'N',Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),TRodd,Ndim, Szatom(:,:,i),Ndim, Dcmplx(0d0,0d0),TRoddSz,Ndim)
                                        CALL zgemm ('N', 'N',Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),Szatom(:,:,i),Ndim, TRodd,Ndim, Dcmplx(0d0,0d0),SzTRodd,Ndim)
                                        do ii=1,Ndim
                                        do jj=1,Ndim
                                                if ((zpos1(ii)-0.5d0)*(zpos1(jj)-0.5d0)>=0d0) then!same layer local torque
                                                        Gammax(ii,jj,i)=Dcmplx(0,1d0)*(TRoddSx(ii,jj)-SxTRodd(ii,jj))
                                                        Gammay(ii,jj,i)=Dcmplx(0,1d0)*(TRoddSy(ii,jj)-SyTRodd(ii,jj))
                                                        Gammaz(ii,jj,i)=Dcmplx(0,1d0)*(TRoddSz(ii,jj)-SzTRodd(ii,jj))
                                                endif
                                        enddo
                                        enddo

                                        CALL zhemm ('L', 'U', Ndim,Ndim,Dcmplx(1d0,0d0),Gammax(:,:,i), Ndim, ham_k, Ndim,Dcmplx(0d0,0d0),GammaxU(:,:,i),Ndim)
                                        CALL zhemm ('L', 'U', Ndim,Ndim,Dcmplx(1d0,0d0),Gammay(:,:,i), Ndim, ham_k, Ndim,Dcmplx(0d0,0d0),GammayU(:,:,i),Ndim)
                                        CALL zhemm ('L', 'U', Ndim,Ndim,Dcmplx(1d0,0d0),Gammaz(:,:,i), Ndim, ham_k, Ndim,Dcmplx(0d0,0d0),GammazU(:,:,i),Ndim)
                                        CALL zgemm ('C', 'N', Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),ham_k,Ndim, GammaxU(:,:,i),Ndim,Dcmplx(0d0,0d0),UGammaxU(:,:,i),Ndim)
                                        CALL zgemm ('C', 'N', Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),ham_k,Ndim, GammayU(:,:,i),Ndim,Dcmplx(0d0,0d0),UGammayU(:,:,i),Ndim)
                                        CALL zgemm ('C', 'N', Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),ham_k,Ndim, GammazU(:,:,i),Ndim,Dcmplx(0d0,0d0),UGammazU(:,:,i),Ndim)

                                enddo
                                !write(3000,"(200f32.16)") (UGammazU) 
                                !~~~We use Kubo formula to calculate torque,
                                !even(odd) part is the imaginary(real) part of
                                !<i|T|j><j|v|i>
                                !VxU=0d0
                                !UVxU=0d0
                                VxU=0d0
                                UVxU=0d0
                                CALL zhemm ('L', 'U', Ndim,Ndim,Dcmplx(1d0,0d0),Vx, Ndim, ham_k, Ndim,Dcmplx(0d0,0d0), VxU,Ndim)
                                CALL zgemm ('C', 'N', Ndim,Ndim,Ndim,Dcmplx(1d0,0d0),ham_k,Ndim, VxU,Ndim, Dcmplx(0d0,0d0),UVxU,Ndim)

                                do i_smr=0,N_smr-1
                                        smr=smr_0+i_smr*smr_step
                                        do i_mu=0,N_mu-1
                                                ef=ef_0-2.0d0+i_mu*mu_step
                                                !if (ef>-0.3d0 .and. ef<0.3d0) ef=ef+0.6d0                                             

                                                do whichband=1,Ndim
                                                        Fermifac(whichband)=1d0/(EXP((Eigenvalue(whichband)-ef)/Temp)+1d0)
                                                !~~~Time-reversal odd has the
                                                !-df/de=0.5*beta/(1+cosh(beta(e-\mu)))
                                                        fac_odd(whichband)=0.5d0/(1d0+cosh((Eigenvalue(whichband)-ef)/Temp))
                                                enddo
                                                do ii=1, Ndim
                                                do jj=1, Ndim
                                                        if (jj==ii) then
                                                                fac(ii,jj)=0d0
                                                        else
                                                                fac(ii,jj)=1d0/(((Eigenvalue(ii)-Eigenvalue(jj))**2+smr**2))
                                                        endif
                                                enddo
                                                enddo
                                                do i=1,2!nat
                                                do ii=1, Ndim
                                                        do jj=1, Ndim           
                                                              tgammax(i_mu+1,i_smr+1,i)=tgammax(i_mu+1,i_smr+1,i)+fac(ii,jj)*(-2d0*aimag(UGammaxU(ii,jj,i)*UVxU(jj,ii)))*Fermifac(ii)
                                                              tgammay(i_mu+1,i_smr+1,i)=tgammay(i_mu+1,i_smr+1,i)+fac(ii,jj)*(-2d0*aimag(UGammayU(ii,jj,i)*UVxU(jj,ii)))*Fermifac(ii)
                                                              tgammaz(i_mu+1,i_smr+1,i)=tgammaz(i_mu+1,i_smr+1,i)+fac(ii,jj)*(-2d0*aimag(UGammazU(ii,jj,i)*UVxU(jj,ii)))*Fermifac(ii)                                                         
                                                        enddo
                                                        tgammaoddx(i_mu+1,i_smr+1,i)=tgammaoddx(i_mu+1,i_smr+1,i)+fac_odd(ii)*real(UGammaxU(ii,ii,i)*UVxU(ii,ii))
                                                        tgammaoddy(i_mu+1,i_smr+1,i)=tgammaoddy(i_mu+1,i_smr+1,i)+fac_odd(ii)*real(UGammayU(ii,ii,i)*UVxU(ii,ii))
                                                        tgammaoddz(i_mu+1,i_smr+1,i)=tgammaoddz(i_mu+1,i_smr+1,i)+fac_odd(ii)*real(UGammazU(ii,ii,i)*UVxU(ii,ii))
                                                enddo
                                                enddo

                                                do kk=1,Ndim
                                                    tetotal(i_mu+1,i_smr+1)=tetotal(i_mu+1,i_smr+1)+1d0/(EXP((Eigenvalue(kk)-ef)/Temp)+1d0)*Eigenvalue(kk)
                                                    density(i_mu+1,i_smr+1)=density(i_mu+1,i_smr+1)+1d0/(EXP((Eigenvalue(kk)-ef)/Temp)+1d0) 
                                                enddo           
                                        enddo
                                enddo
                                kcounter=kcounter+1
                                !if (rank==0 .and.mod(kcounter,Nkx*Nky/size/10)==0) then !Nkx*Nky*Nkz/size/100
                                 !       call cpu_time(finish)
                                  !      write(*,'("node=",i4,"process=",f12.4,"using time=",3f12.4)') rank,real(kcounter)*size/(Nkx*Nky),finish-start        
                                !endif
                                endif

                          enddo
                          enddo

        call cpu_time(finish)
        if (rank==0) write(*,'("Solve h_k Time = ",f10.3," seconds.")') finish-start
        allocate(sum_gammax(2),sum_gammay(2),sum_gammaz(2),sum_gammaoddx(2),sum_gammaoddy(2),sum_gammaoddz(2),stat=ierr)


        do i_smr=0,N_smr-1
                smr=smr_0+i_smr*smr_step
                do i_mu=0,N_mu-1
                        ef=ef_0-2.0d0+i_mu*mu_step
                        !if (ef>-0.3d0 .and. ef<0.3d0) ef=ef+0.6d0                                               
                        
                        sum_etotal=0d0
                        sum_density=0d0
                        sum_gammax=0d0
                        sum_gammay=0d0
                        sum_gammaz=0d0                  
                        CALL MPI_REDUCE(tetotal(i_mu+1,i_smr+1), sum_etotal, 1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
                        CALL MPI_REDUCE(density(i_mu+1,i_smr+1), sum_density, 1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)                   
                        do i=1,2!nat
                                CALL MPI_REDUCE(tgammax(i_mu+1,i_smr+1,i),sum_gammax(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
                                CALL MPI_REDUCE(tgammay(i_mu+1,i_smr+1,i),sum_gammay(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
                                CALL MPI_REDUCE(tgammaz(i_mu+1,i_smr+1,i),sum_gammaz(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)                                 
                                CALL MPI_REDUCE(tgammaoddx(i_mu+1,i_smr+1,i),sum_gammaoddx(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
                                CALL MPI_REDUCE(tgammaoddy(i_mu+1,i_smr+1,i),sum_gammaoddy(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
                                CALL MPI_REDUCE(tgammaoddz(i_mu+1,i_smr+1,i),sum_gammaoddz(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)                                      
                        enddo
                        if(rank .eq. 0) then
                                !write(*,*) "reduce sum",sum_etotal
                                sum_etotal=sum_etotal*kweight
                                sum_density=sum_density*kweight
                                do i=1,2!nat
                                        sum_gammax(i)=sum_gammax(i)*kweight
                                        sum_gammay(i)=sum_gammay(i)*kweight
                                        sum_gammaz(i)=sum_gammaz(i)*kweight
                                        sum_gammaoddx(i)=sum_gammaoddx(i)*kweight/Temp/smr
                                        sum_gammaoddy(i)=sum_gammaoddy(i)*kweight/Temp/smr
                                        sum_gammaoddz(i)=sum_gammaoddz(i)*kweight/Temp/smr              
                                enddo
                                write(1000,"(200f32.16)") angle_theta(1),angle_theta(2),angle_phi(1),angle_phi(2),smr,ef,sum_etotal,sum_density,&
                                &sum_gammax,sum_gammay,sum_gammaz,sum_gammaoddx,sum_gammaoddy,sum_gammaoddz!,sum_gammayx,sum_gammayy,sum_gammayz,sum_gammazx,sum_gammazy,sum_gammazz
        
                        endif                           
       if(rank .eq. 0) then
                  write(*,*) "Torque_A=",sum_gammax(1),sum_gammay(1),sum_gammaz(1),sum_gammaoddx(1),sum_gammaoddy(1),sum_gammaoddz(1)
                  write(*,*) "Torque_B=",sum_gammax(2),sum_gammay(2),sum_gammaz(2),sum_gammaoddx(2),sum_gammaoddy(2),sum_gammaoddz(2)
        endif
!        enddo

                enddo
        enddo
        enddo
        enddo
        ! Finalize 
    call MPI_FINALIZE(ierror)   
END Program MoTe2_fromWF
