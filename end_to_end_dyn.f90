!*************************************************************************************!

!Precision_control 
!*************************************************************************************!    
      module Precision_Control
      implicit none
      integer, parameter :: single = selected_real_kind(p=6,r=37) !single precision
      integer, parameter :: double = selected_real_kind(p=15,r=37) !double precision
      integer, parameter :: quad= selected_real_kind(p=33,r=4931)  !quadruple precision
      save
      end module


!*************************************************************************************!

!Parameter definition
!*************************************************************************************! 

      module Param_Control
      Use Precision_Control
      implicit none
      integer:: step_num, Nstep,Nstepskip, Npolymer, Nresi, Nbead_per_resi, Ndim
      real(kind=double) ::rho

      CHARACTER(LEN=100),PARAMETER ::FMT1="(3f8.3)"
      save
      end module

!*************************************************************************************!

!Variable definition
!*************************************************************************************! 
      module Var_Control
      Use Precision_Control
      Use Param_Control
      implicit none
      integer ::istep
      real(kind=double), allocatable, dimension(:) :: x, y, z
      real(kind=double)::boxlen, binsize
      real(kind=double), allocatable, dimension(:,:) :: com_lys,com_glu,com
      real(kind=double), allocatable, dimension(:,:) ::com_x,com_y,com_z
      real(kind=double), allocatable, dimension(:,:) ::com_droplet_x,com_droplet_y, com_droplet_z
      integer::Nbin
      real(kind=double),allocatable,dimension(:):: unnorm_summed_numbers
      real(kind=double)::av_rg_droplet, av_inter_droplet,av_intra_droplet
      real(kind=double):: av_intra_droplet_hyd, av_inter_droplet_hyd
      real(kind=double):: av_rg_free, av_inter_free, av_intra_free
      real(kind=double):: sq_rg_droplet_av, sq_inter_droplet_av,sq_intra_droplet_av
      real(kind=double):: sq_rg_free_av, sq_inter_free_av,sq_intra_free_av
      real(kind=double)::sq_inter_droplet_hyd_av, sq_intra_droplet_hyd_av
      integer, allocatable,dimension(:) ::  num_neighbor
      integer, allocatable,dimension(:,:)::neighborlist
      integer, allocatable,dimension(:) ::cluster_head, cluster_tail, cluster_size, cluster_list, cluster_checked
      integer,allocatable, dimension(:) :: state
      integer ::clustersum, av_droplet
      integer,allocatable,dimension(:,:)::listed_id,id_droplet,id_free_single
      integer,allocatable,dimension(:,:) ::id_particle_in_cluster
      integer:: count_number_droplet, count_number_free
      character(len=1), allocatable,dimension(:)::myseq,seq_char
      real(kind=double) ::energy_electrostatic_per_polymer,energy_hydrophobic_per_polymer
      real(kind=double)::sq_elec, sq_hyd
      integer::largest_size_of_cluster
      real(kind=double), allocatable, dimension(:,:)::distx,disty,distz,d_matrix
      save 
      end module

!*************************************************************************************!

!CORE
!*************************************************************************************! 
      program main
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer::  iresi,ibin,ip,xx
      integer ::ijunk
      real(kind=double) :: constant, vol, pi =3.14, num_ideal
      call INITIALIZE
      !do iresi = 1, Nbead_per_resi
      !   read(2,'(a1)') myseq(iresi)
      !end do
      !do ip =1, Npolymer
      !  do iresi=1, Nbead_per_resi
      !     xx=((ip-1)*Nbead_per_resi) + iresi
      !     seq_char(xx) = myseq(iresi)
      !  end do
      !end do
      
      read(1,*)
      read(1,*)
      do iresi = 1 , Nresi
         read(1,*) 
      end do
      
      do ijunk = 1, Nstepskip*4002
         read(1,*)
      end do
      
      do istep =1 , (Nstep-Nstepskip)
         
         call READ_FILE
         call FINDCOM
         call RDF
         call GETCLUSTER
         call GETRG
      end do
      
      constant = (1.33333333333)* pi
      do ibin = 1, Nbin
         vol = constant*(((ibin*binsize)**3)-((ibin-1)*binsize)**3)
         num_ideal = vol*rho
         
         unnorm_summed_numbers(ibin)=unnorm_summed_numbers(ibin)/(num_ideal)
      end do
      do ibin = 1, Nbin
         !write(200,*) unnorm_summed_numbers(ibin)
         !write(100, *)ibin*binsize,unnorm_summed_numbers(ibin)/(Npolymer*(Nstep-Nstepskip))
      end do
      !write(1500,*) real(av_droplet)/(Nstep-Nstepskip)

      !do ip =1, Npolymer
      !   write(2000,*) id_droplet(ip,Nstep-Nstepskip)
      !end do
      !call GETDIFFUSION
      call EXCHANGE_NBRS
      end program main

!*************************************************************************************!

!INITIALIZE
!*************************************************************************************! 
      subroutine INITIALIZE
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      Nstep  = 15500
      Nstepskip =5500
      !open (unit = 2 ,file &
      ! ="Nstepinfo", &
      !                                           action="read")
      !read (2,*) Nstep
      !Nstep =Nstep-100
      !Nstepskip=Nstep-10000
      Npolymer = 100
      Nbead_per_resi = 40
      Nresi = Npolymer*Nbead_per_resi
      Ndim           =3
      boxlen = 302.38
      Nbin = 500
      binsize = boxlen/(Nbin)
      rho     = Nresi/(boxlen**3)
      allocate(x(Nresi), y(Nresi), z(Nresi))
      allocate(com_lys(Npolymer, Ndim), com_glu(Npolymer,Ndim),com(Npolymer, Ndim))
      allocate(com_x(Npolymer,(Nstep-Nstepskip)))
      allocate(com_y(Npolymer,(Nstep-Nstepskip)))
      allocate(com_z(Npolymer,(Nstep-Nstepskip)))
      allocate(unnorm_summed_numbers(Nbin+1))

      unnorm_summed_numbers = 0
      av_droplet= 0
      

      av_rg_droplet=0
      av_inter_droplet =0
      av_intra_droplet =0
      av_rg_free = 0
      av_inter_free=0
      av_intra_free =0
      av_intra_droplet_hyd=0 
      av_inter_droplet_hyd=0
      count_number_droplet =0
      count_number_free    =0
      sq_rg_droplet_av=0
      sq_inter_droplet_av=0
      sq_intra_droplet_av=0
      sq_rg_free_av=0
      sq_inter_free_av=0
      sq_intra_free_av=0
      sq_inter_droplet_hyd_av=0 
      sq_intra_droplet_hyd_av=0
      energy_electrostatic_per_polymer=0
      sq_elec =0
      energy_hydrophobic_per_polymer=0
      sq_hyd =0

      allocate(num_neighbor(Npolymer),neighborlist(Npolymer,(Npolymer-1)))
      allocate(state(Npolymer))
      allocate(cluster_head(Npolymer),cluster_tail(Npolymer),cluster_size(Npolymer), cluster_checked(Npolymer))
      allocate (cluster_list(Npolymer))
      allocate(listed_id(Npolymer,Npolymer), id_droplet(Npolymer,(Nstep-Nstepskip)))
      allocate(id_free_single(Npolymer, (Nstep-Nstepskip)))
      allocate(com_droplet_x(Npolymer,(Nstep-Nstepskip)))
      allocate(com_droplet_y(Npolymer, (Nstep-Nstepskip)))
      allocate(com_droplet_z(Npolymer, (Nstep-Nstepskip)))
      allocate(id_particle_in_cluster(Npolymer,(Nstep-Nstepskip)))
      allocate(myseq(Nbead_per_resi))
      allocate(seq_char(Nresi))
      allocate(distx(Npolymer, (Nstep-Nstepskip)))
      allocate(disty(Npolymer, (Nstep-Nstepskip)))
      allocate(distz(Npolymer, (Nstep-Nstepskip)))
      allocate(d_matrix(Npolymer, (Nstep-Nstepskip)))
      open (unit = 1 ,file &
       ="temp", &
                                                 action="read")
      !open (unit = 2 ,file="myseq.dat", &
      !                                           action="read")
      open (unit = 5000 ,file="energy_elec.dat", &
                                                 action="write")
      open (unit = 6000 ,file="energy_hyd.dat", &
                                                 action="write")
      end subroutine INITIALIZE

!*************************************************************************************!

!READ TRAJECTORY
!*************************************************************************************! 

      subroutine READ_FILE
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer::  ip,iresi,xx
         read(1,*)
         read(1,"(i12)") step_num
         do iresi = 1, Nresi
            read(1,FMT1) x(iresi),y(iresi),z(iresi)
         end do

         
         
      !print*, x(Nresi),y (Nresi),z(Nresi)
      end subroutine READ_FILE


!*************************************************************************************!

!FINDING COM FOR POSITIVE AND NEGATIVE RESIDUES OF EACH CHAIN
!*************************************************************************************! 

      subroutine FINDCOM
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer::ii, ip,ll, iresi
      do ip =1 , Npolymer
         
         do ii = 1,3
            com_lys(ip,ii) =0
            com_glu(ip,ii) =0
            com(ip, ii)    =0
         end do
      end do
      do ip = 1, Npolymer
         ll = (ip-1)*Nbead_per_resi 
         do iresi = ll+1 , ll+10 
            com_lys(ip,1)= com_lys(ip,1)+x(iresi) 
            com_lys(ip,2)= com_lys(ip,2)+y(iresi)
            com_lys(ip,3)= com_lys(ip,3)+z(iresi)
         end do
         do iresi = ll +11, ll+20
            com_glu(ip,1) = com_glu(ip,1)+x(iresi)
            com_glu(ip,2) = com_glu(ip,2)+y(iresi)
            com_glu(ip,3) = com_glu(ip,3)+z(iresi)
         end do
      end do
      do ip = 1, Npolymer
         ll = (ip-1)*Nbead_per_resi
         do iresi = ll+1 , ll+Nbead_per_resi
            com(ip,1)= com(ip,1)+x(iresi)
            com(ip,2)= com(ip,2)+y(iresi)
            com(ip,3)= com(ip,3)+z(iresi)
         end do
      end do
      do ip = 1, Npolymer
         do  ii = 1,3
             com(ip,ii) =com(ip,ii)/Nbead_per_resi
             com_lys(ip,ii)=(2*com_lys(ip,ii)/Nbead_per_resi)
             com_glu(ip,ii) = (2*com_glu(ip,ii)/Nbead_per_resi)
         end do
      end do
      do ip =1, Npolymer
         com_x(ip,istep)=com(ip,1)
         com_y(ip, istep) = com(ip,2)
         com_z(ip, istep) = com(ip,3)
      end do
      end subroutine FINDCOM

!*************************************************************************************!

!BEAD to BEAD RDF
!*************************************************************************************! 

      subroutine RDF
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      real(kind=double) ::d(3),mod_d, pi=3.14, constant
      integer::iresi,jresi, whichbin
      integer:: ip,jp,ll,mm
      constant = (4/3)* pi
      do ip = 1, Npolymer-1
       do jp = ip+1, Npolymer
          if (ip .eq. jp) cycle
                d(1) = com(jp,1) -com(ip,1)
                d(2) = com(jp,2) -com(ip,2)
                d(3) = com(jp,3) -com(ip,3)
                mod_d = sqrt(d(1)**2 +d(2)**2 +d(3)**2)
                if (mod_d .gt. 3.7 .and. mod_d .lt. boxlen) then
                   !if(mod_d .lt. 3.8) then
                   !  write(500,*) iresi, jresi, mod_d
                   !end if
                   whichbin = int(mod_d/binsize)
                   unnorm_summed_numbers(whichbin+1) = unnorm_summed_numbers(whichbin+1)+2
                end if

       end do
      end do 
      end subroutine RDF
      
!*************************************************************************************!

!CLUSTERING TO SEE NEXT
!*************************************************************************************!

      subroutine GETCLUSTER
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer :: max_neighbor 
      integer :: ip, jp, pp, kk, mm, nbrs, nn,zz,qq,ll,rr
      integer :: ic,jc,ii
      real(kind=double) ::comij_x,comij_y, comij_z, modrij 
      real(kind=double):: fraction_droplet
      logical:: finding_not_completed
      integer :: temp, largest_cluster, droplet_polymer
      !Initializing at every time step
      clustersum = 0
      cluster_size = 0
      cluster_head  =0
      cluster_tail =0
      cluster_list =0
      cluster_checked =0 
      droplet_polymer=0
      do ip = 1, Npolymer
         do jp = 1, Npolymer-1
            neighborlist(ip, jp) = 0
         end do
      end do
      do ip = 1, Npolymer
         state(ip)        = 0
         num_neighbor(ip) = 0
      end do
      do ip =1, Npolymer
       do jp =1, Npolymer
          listed_id(ip,jp) =0
       end do
      end do
      
      !way to find the clusters
      !create the neighbor list 
      !go along the neighbor list to find all the neighbors of neighbors
      !which have not been counted to any cluster till now
      max_neighbor= Npolymer

      !search the neighbors
      do ip = 1, Npolymer -1
         do jp = ip+1, Npolymer
            comij_x = com(jp,1) -com(ip,1)
            comij_y = com(jp,2) -com(ip,2)
            comij_z = com(jp,3) -com(ip,3)
            modrij = sqrt (comij_x**2 + comij_y**2 + comij_z**2)
            if (modrij .le. 30) then 
               num_neighbor(ip) = num_neighbor(ip)+1
               num_neighbor(jp) = num_neighbor(jp)+1
               neighborlist(ip, num_neighbor(ip)) = jp
               neighborlist(jp, num_neighbor(jp)) = ip
            end if 
         end do
      end do
      
       
      
      !you have made your neighbor list 
      !now go on to find the clusters
      do ip = 1, Npolymer
         if (state(ip) .eq. 0) then
            clustersum =clustersum+1
            pp = ip
            state(pp) = clustersum
            cluster_head(clustersum) = pp
            cluster_tail(clustersum) = pp
            cluster_size(clustersum) = 1
            listed_id(clustersum, cluster_size(clustersum)) =pp
            id_particle_in_cluster(ip,istep) =clustersum
            !defined a logical variable in order to understand when to
            !end and when to start
            finding_not_completed  = .TRUE.
            do while (finding_not_completed)
               nbrs = num_neighbor(pp)
               kk = 1
               do while (kk .le. nbrs)
                  mm = neighborlist(pp, kk)
                  if (state(mm) .eq. 0) then 
                     cluster_list(cluster_tail(clustersum))= mm
                     cluster_tail(clustersum)              = mm
                     cluster_size(clustersum)              =cluster_size(clustersum)+1
                     state(mm) = clustersum
                     listed_id(clustersum, cluster_size(clustersum)) =mm
                     id_particle_in_cluster(mm,istep) = clustersum
                     !print*, listed_id(clustersum,cluster_size(clustersum))
                  end if
                  kk =kk+1
               end do
               cluster_checked(clustersum) = pp
               if (pp .eq.cluster_tail(clustersum)) then
                  finding_not_completed = .FALSE.
                  cluster_list(pp) = 0
               end if
               pp = cluster_list(pp)
            end do
         end if
      end do
      !print*, com(100,1)
      !do ic = 1, clustersum
      !   write(1100, *) cluster_size(ic)
      !end do
      do ic = 1, clustersum
         if(cluster_size(ic) .gt. 3) then
           droplet_polymer = droplet_polymer +cluster_size(ic)
         end if
      end do
      av_droplet= av_droplet+droplet_polymer
      fraction_droplet = droplet_polymer/Npolymer 
      largest_cluster= maxval(cluster_size)
      largest_size_of_cluster=largest_cluster 
       do ic = 1, clustersum
         if( cluster_size(ic) .eq. largest_cluster) then
           zz = cluster_size(ic)
           !print*, zz
           do nn =1, zz
              id_droplet(listed_id(ic,nn),istep) =1
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)
              
           end do
         else 
          zz=cluster_size(ic)
          do nn =1, zz
              id_droplet(listed_id(ic,nn),istep) =0
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)
           end do
         end if
       end do

       do ic = 1, clustersum
         if( cluster_size(ic) .lt. 2) then
           zz = cluster_size(ic)
           !print*, zz
           do nn =1, zz
              id_free_single(listed_id(ic,nn),istep) =1
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)

           end do
         else 
           zz = cluster_size(ic)
           do nn =1, zz
              id_free_single(listed_id(ic,nn),istep) =0
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)

           end do
         end if
       end do





       !calculate the COM of droplet
       do ic = 1, clustersum

          do nn=1,cluster_size(ic)
             ll=listed_id(ic,nn)
             !qq=(ll-1)*Nbead_per_resi
             !do rr=1, Nbead_per_resi

              !  com_droplet_x(ic,istep)=com_droplet_x(ic,istep)+x(qq+rr)
              !  com_droplet_y(ic,istep)=com_droplet_y(ic,istep)+y(qq+rr)
              !  com_droplet_z(ic,istep)=com_droplet_z(ic,istep)+z(qq+rr)
               com_droplet_x(ic,istep)=com_droplet_x(ic,istep)+com_x(listed_id(ic,nn),istep)
               com_droplet_y(ic,istep)=com_droplet_y(ic,istep)+com_y(listed_id(ic,nn),istep)
               com_droplet_z(ic,istep)=com_droplet_z(ic,istep)+com_z(listed_id(ic,nn),istep)
             !end do
          end do
          !com_droplet_x(ic,istep)=com_droplet_x(ic,istep)/(Nbead_per_resi*cluster_size(ic))
          !com_droplet_y(ic,istep)=com_droplet_y(ic,istep)/(Nbead_per_resi*cluster_size(ic))
          !com_droplet_z(ic,istep)=com_droplet_z(ic,istep)/(Nbead_per_resi*cluster_size(ic))
          com_droplet_x(ic,istep)=com_droplet_x(ic,istep)/(cluster_size(ic))
          com_droplet_y(ic,istep)=com_droplet_y(ic,istep)/(cluster_size(ic))
          com_droplet_z(ic,istep)=com_droplet_z(ic,istep)/(cluster_size(ic))

       end do
 



       !Now find out the largest cluster

     


     !  do ic = 1, clustersum
     !    do jc = ic, clustersum 
     !       if(cluster_size(ic) .lt. cluster_size(jc)) then
     !         temp = cluster_size(jc)
     !         cluster_size(jc) = cluster_size(ic)
     !         cluster_size(ic) = temp
     !       end if
     !    end do
     ! end do
     ! largest_cluster = cluster_size(1)        
      !write(700,*) largest_cluster, clustersum
      !write(700,*) "stepdone=", istep 
      !write(800,*) step_num,real(Npolymer-droplet_polymer)/Npolymer,real(droplet_polymer)/Npolymer
      !write(900,*) step_num,largest_cluster

             
           
         

      end subroutine GETCLUSTER

!*************************************************************************************!

!Radius of gyration of monomers in droplet
!*************************************************************************************!
      subroutine GETRG
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      logical::status_go
      real(kind=double)::d(3), d_sq, sum_dis_sq
      real(kind=double):: d_from_droplet_com(3), dijbead(3),dij,dijintra(3), dijin
      real(kind=double)::rg_monomer
      integer ::ip,jp,ll,mm,iresi,jresi,id_poly,qiresi,qjresi
      real(kind=double)::Ecol,eps_solvent=80.0,salt_conc=0.02,kappa,tot_E_col, tot_E_hyd, Ehyd
      real(kind=double)::r12,r10,epsil=0.2, sigma=7.0
      integer :: count_free_pts, test_var
      real(kind=double) ::av_cn_free
      logical:: state_poly
      do ip = 1, Npolymer
       if (id_droplet(ip,istep).eq. 1 ) then
         ll = (ip-1)*Nbead_per_resi
         
           
            iresi = ll+1

            jresi = ll+Nbead_per_resi
               
            dijintra(1) = x(jresi)-x(iresi)
            dijintra(2) = y(jresi)-y(iresi)
            dijintra(3) = z(jresi)-z(iresi)
            dijin= sqrt(dijintra(1)**2+ dijintra(2)**2+dijintra(3)**2)
            distx(ip, istep) =dijintra(1)
            disty(ip, istep) =dijintra(2)
            distz(ip, istep) =dijintra(3)
            d_matrix(ip,istep)=dijin
            write(500,*) dijin   
       end if
      end do
      end subroutine GETRG
!*************************************************************************************!

!time correlation
!*************************************************************************************!
      subroutine EXCHANGE_NBRS
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      real(kind=double),allocatable,dimension(:) :: stay_connect, d_av
      integer::ip,jp,iframe,jframe,kframe,corr_frame,count_such_event,count_ip
      real(kind=double):: norm,tt,rr,zz
      logical::status_jump, status_go
      corr_frame=int((Nstep-Nstepskip)/5)
      allocate(stay_connect(corr_frame))
      allocate (d_av(Npolymer))
      print*,"entered"
      !print*, dist(10,1)
      do ip =1, Npolymer
         count_ip = 0
         do iframe = 1, (Nstep-Nstepskip)
            if (id_droplet(ip,iframe).eq. 1 ) then
                d_av(ip) = d_av(ip)+d_matrix(ip,iframe)
                count_ip = count_ip+1
            end if
         end do
         d_av(ip) = d_av(ip)/count_ip
      end do

      norm =0
      do iframe =1, corr_frame

         stay_connect(iframe) =0
      end do
      do jframe =1, (Nstep-Nstepskip)
         do ip =1, Npolymer
              
               if(id_droplet(ip,jframe).eq. 1) then
                  !norm = norm + distx(ip, jframe)**2 + disty(ip, jframe)**2 + distz(ip, jframe)**2
                   norm = norm +(d_matrix(ip,jframe)-d_av(ip))*(d_matrix(ip,jframe)-d_av(ip))
               end if
               
         end do
      end do
      print*, norm
      do iframe =1, corr_frame
         print*, "entered frame=", iframe
         do jframe =1, (Nstep-Nstepskip)
            if (jframe+iframe .gt. (Nstep-Nstepskip)) cycle
            do ip =1, Npolymer
               if(id_droplet(ip,jframe).eq. 1 .and. id_droplet(ip,jframe+iframe).eq. 1) then
                 status_go=.true.
                 do kframe = jframe+1, jframe+iframe-1
                    if (id_droplet(ip,kframe) .eq. 0) then
                     status_go = .false.
                     exit
                    end if
                 end do
                 if( status_go .eqv. .false.) then
                   cycle
                 else
                   tt = distx(ip, jframe)* distx(ip,jframe+iframe) 
                   rr = disty(ip, jframe)* disty(ip,jframe+iframe)
                   zz=distz(ip, jframe)* distz(ip,jframe+iframe)
                   !stay_connect(iframe) = stay_connect(iframe)+tt+rr+zz
                   stay_connect(iframe) =stay_connect(iframe)+(d_matrix(ip,jframe)-d_av(ip))*(d_matrix(ip,(jframe+iframe))-d_av(ip))
                 end if
               end if
            end do
         end do
         !print*, "done with=", iframe
      end do
      do iframe = 1, corr_frame
         write(2000,*) iframe, stay_connect(iframe)/norm
      end do
      end subroutine EXCHANGE_NBRS
